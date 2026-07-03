/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <config/config.hpp>

#ifdef ESPRESSO_DIPOLES

#include "magnetostatics/dipolar_direct_sum.hpp"

#include "BoxGeometry.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include "magnetostatics/dipolar_direct_sum_kernels.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>
#include <Kokkos_ScatterView.hpp>

#include <utils/Vector.hpp>
#include <utils/mpi/iall_gatherv.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

/**
 * @brief Position and dipole moment of one particle.
 */
struct PosMom {
  Utils::Vector3d pos;
  Utils::Vector3d m;

  template <class Archive> void serialize(Archive &ar, long int) {
    ar & pos & m;
  }
};

static auto gather_particle_data(BoxGeometry const &box_geo,
                                 ParticleRange const &particles) {
  auto const &comm = ::comm_cart;
  std::vector<Particle *> local_particles;
  std::vector<PosMom> local_posmom;
  std::vector<PosMom> all_posmom;
  std::vector<boost::mpi::request> reqs;

  local_particles.reserve(particles.size());
  local_posmom.reserve(particles.size());

  for (auto &p : particles) {
    if (p.dipm() != 0.0) {
      local_particles.emplace_back(&p);
      local_posmom.emplace_back(
          PosMom{box_geo.folded_position(p.pos()), p.calc_dip()});
    }
  }

  auto const local_size = static_cast<int>(local_posmom.size());
  std::vector<int> all_sizes;
  boost::mpi::all_gather(comm, local_size, all_sizes);

  auto const offset =
      std::accumulate(all_sizes.begin(), all_sizes.begin() + comm.rank(), 0);
  auto const total_size =
      std::accumulate(all_sizes.begin() + comm.rank(), all_sizes.end(), offset);

  if (comm.size() > 1) {
    all_posmom.resize(total_size);
    reqs = Utils::Mpi::iall_gatherv(comm, local_posmom.data(), local_size,
                                    all_posmom.data(), all_sizes.data());
  } else {
    std::swap(all_posmom, local_posmom);
  }

  return std::make_tuple(std::move(local_particles), std::move(all_posmom),
                         std::move(reqs), offset);
}

static auto get_n_cut(BoxGeometry const &box_geo, int n_replicas) {
  return n_replicas * Utils::Vector3i{static_cast<int>(box_geo.periodic(0)),
                                      static_cast<int>(box_geo.periodic(1)),
                                      static_cast<int>(box_geo.periodic(2))};
}

struct PosMomViews {
  Kokkos::View<double *[3], Kokkos::LayoutLeft, Kokkos::HostSpace> pos;
  Kokkos::View<double *[3], Kokkos::LayoutLeft, Kokkos::HostSpace> m;
};

static PosMomViews make_posmom_views(std::size_t n_total) {
  return {decltype(PosMomViews::pos)("dds_pos", n_total),
          decltype(PosMomViews::m)("dds_m", n_total)};
}

static void fill_posmom_views(PosMomViews &views,
                              std::vector<PosMom> const &all_posmom,
                              std::size_t begin, std::size_t end) {
  for (auto i = begin; i < end; ++i) {
    for (int c = 0; c < 3; ++c) {
      views.pos(i, c) = all_posmom[i].pos[c];
      views.m(i, c) = all_posmom[i].m[c];
    }
  }
}

/** Broadcast a scalar 3-vector into a 3-vector of simd registers. */
static Utils::Vector<simd_double, 3> broadcast_simd(Utils::Vector3d const &v) {
  return {simd_double(v[0]), simd_double(v[1]), simd_double(v[2])};
}

/** Load the dipole moments of lanes j..j+width-1 into a simd 3-vector.
 *  LayoutLeft makes component @c c contiguous across particles, so a simd
 *  register can be constructed directly from @c &views.m(j, c). */
static Utils::Vector<simd_double, 3> load_simd_moment(PosMomViews const &views,
                                                      std::size_t j) {
  Utils::Vector<simd_double, 3> m;
  for (int c = 0; c < 3; ++c)
    m[c] = simd_double(&views.m(j, c), Kokkos::Experimental::simd_flag_default);
  return m;
}

/** Load the primary distance vectors of lanes j..j+width-1 relative to
 *  @c pos_i. With replicas: full SIMD (broadcast pos_i, load pos_j, subtract).
 *  Without replicas: per-lane minimum image (preserves Lees-Edwards
 *  semantics), scalarizing only the distance computation. */
static Utils::Vector<simd_double, 3>
primary_distance_simd(Utils::Vector3d const &pos_i, PosMomViews const &views,
                      std::size_t j, bool with_replicas,
                      BoxGeometry const &box_geo) {
  constexpr std::size_t w = simd_double::size();
  Utils::Vector<simd_double, 3> d;
  if (with_replicas) {
    for (int c = 0; c < 3; ++c) {
      simd_double const pj(&views.pos(j, c),
                           Kokkos::Experimental::simd_flag_default);
      d[c] = simd_double(pos_i[c]) - pj;
    }
  } else if (box_geo.type() == BoxType::CUBOID) {
    /* Vectorized, branchless minimum image, bit-exact with
     * detail::get_mi_coord: wrap = dx - round(dx * box_l_inv) * box_l is
     * selected per lane only where |dx| > box_l_half (strictly greater, as in
     * the scalar reference), leaving dx untouched otherwise. Only applied to
     * periodic dimensions. */
    auto const &box_l = box_geo.length();
    auto const &box_l_inv = box_geo.length_inv();
    auto const &box_l_half = box_geo.length_half();
    for (int c = 0; c < 3; ++c) {
      simd_double const pj(&views.pos(j, c),
                           Kokkos::Experimental::simd_flag_default);
      simd_double dx = simd_double(pos_i[c]) - pj;
      if (box_geo.periodic(static_cast<unsigned>(c))) {
        simd_double const wrapped =
            dx - Kokkos::round(dx * simd_double(box_l_inv[c])) *
                     simd_double(box_l[c]);
        auto const m = Kokkos::abs(dx) > simd_double(box_l_half[c]);
        /* select per lane: wrapped where |dx| > box_l_half, else dx */
        dx = Kokkos::Experimental::condition(m, wrapped, dx);
      }
      d[c] = dx;
    }
  } else {
    /* Lees-Edwards (or any non-cuboid): per-lane scalar fallback preserves the
     * shear min-image semantics of box_geo.get_mi_vector. */
    double buf[3][w];
    for (std::size_t l = 0; l < w; ++l) {
      Utils::Vector3d const pos_j{views.pos(j + l, 0), views.pos(j + l, 1),
                                  views.pos(j + l, 2)};
      auto const mi = box_geo.get_mi_vector(pos_i, pos_j);
      for (int c = 0; c < 3; ++c)
        buf[c][l] = mi[c];
    }
    for (int c = 0; c < 3; ++c)
      d[c] = simd_double(&buf[c][0], Kokkos::Experimental::simd_flag_default);
  }
  return d;
}

/** Real-space image shifts n .* box_l inside the |ncut| sphere; index 0 is the
 *  primary (zero) shift so self-interaction loops start at index 1. */
static std::vector<Utils::Vector3d>
make_image_shifts(Utils::Vector3i const &ncut, Utils::Vector3d const &box_l) {
  auto const ncut2 = ncut.norm2();
  std::vector<Utils::Vector3d> shifts;
  shifts.push_back({0., 0., 0.});
  for (int nx = -ncut[0]; nx <= ncut[0]; ++nx)
    for (int ny = -ncut[1]; ny <= ncut[1]; ++ny)
      for (int nz = -ncut[2]; nz <= ncut[2]; ++nz) {
        if (nx == 0 && ny == 0 && nz == 0)
          continue;
        if (nx * nx + ny * ny + nz * nz <= ncut2)
          shifts.push_back({nx * box_l[0], ny * box_l[1], nz * box_l[2]});
      }
  return shifts;
}

/**
 * @brief Calculate and add the interaction forces/torques to the particles.
 *
 * This employs a parallel N-square loop over all particle pairs.
 * The computation the partitioned into several steps so that the
 * communication latency can be hidden behind some local computation:
 *
 * 1. The local particle positions and momenta are packed into
 *    one array.
 * 2. The asynchronous distribution of the local arrays to all
 *    ranks is started.
 * 3. The interaction for the local pairs is started, here every
 *    pair is visited only once, and the force is added to both particles.
 * 4. Wait for the data from the other nodes.
 * 5. Calculate the interaction with the rest of the particles. Here
 *    every pair is visited twice (not necessarily on the same rank)
 *    so that no reduction of the forces is needed.
 *
 * Logically this is equivalent to the potential calculation
 * in @ref DipolarDirectSum::long_range_energy, which calculates
 * a naive N-square sum, but has better performance and scaling.
 */
void DipolarDirectSum::add_long_range_forces_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &box_l = box_geo.length();
  auto const particles = system.cell_structure->local_particles();
  auto [local_particles, all_posmom, reqs, offset_signed] =
      gather_particle_data(box_geo, particles);

  /* Number of image boxes considered */
  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_l);

  auto const offset = static_cast<std::size_t>(offset_signed);
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();
  constexpr std::size_t w = simd_double::size();

  auto const prefactor_local = prefactor;

  /* SoA views; fill the local slice only so the remote gather can overlap. */
  auto views = make_posmom_views(n_total);
  fill_posmom_views(views, all_posmom, offset, offset + n_local);

  using execution_space = Kokkos::DefaultExecutionSpace;
  using ForceView =
      Kokkos::View<double *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
  using ScatterForce =
      Kokkos::Experimental::ScatterView<double *[3], Kokkos::LayoutRight>;
  ForceView local_force("dds_force", n_local);
  ForceView local_torque("dds_torque", n_local);
  ScatterForce scatter_force(local_force);
  ScatterForce scatter_torque(local_torque);

  /* Raw pointers so the Kokkos lambdas do not capture std::vector by value. */
  auto *local_particles_ptr = local_particles.data();
  auto const *shifts_ptr = shifts.data();
  auto const n_shifts = shifts.size();

  /* Phase A: local pairs. Each i owns its own force/torque accumulation
   * (written directly, unique owner, no race); the Newton's-third-law
   * partner-j contributions go through the ScatterView with a per-lane
   * scatter. */
  Kokkos::parallel_for(
      "dds_local_pairs",
      Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(
          std::size_t{0}, n_local)
          .set_chunk_size(64),
      [=](std::size_t const i) {
        auto const gi = offset + i;
        Utils::Vector3d const pos_i{views.pos(gi, 0), views.pos(gi, 1),
                                    views.pos(gi, 2)};
        Utils::Vector3d const m_i{views.m(gi, 0), views.m(gi, 1),
                                  views.m(gi, 2)};
        PairForce<double> fi{};

        /* (a) self-images (shifts[1..], primary excluded) */
        for (std::size_t s = 1; s < n_shifts; ++s)
          fi += pair_force<double>(shifts_ptr[s], m_i, m_i);

        auto force_access = scatter_force.access();
        auto torque_access = scatter_torque.access();
        auto const m_i_s = broadcast_simd(m_i);

        /* (b) SIMD body over j in (gi, offset + n_local) */
        auto j = gi + 1;
        for (; j + w <= offset + n_local; j += w) {
          auto const m_j = load_simd_moment(views, j);
          auto const d0 =
              primary_distance_simd(pos_i, views, j, with_replicas, box_geo);
          PairForce<simd_double> fij{}, fji{};
          for (std::size_t s = 0; s < n_shifts; ++s) {
            auto const &shift = shifts_ptr[s];
            Utils::Vector<simd_double, 3> const rn{
                d0[0] + shift[0], d0[1] + shift[1], d0[2] + shift[2]};
            auto const pf = pair_force<simd_double>(rn, m_i_s, m_j);
            fij += pf;
            fji.f -= pf.f;
            /* Conservation of angular momentum mandates that
             * 0 = t_i + r_ij x F_ij + t_j */
            fji.torque += vector_product(pf.f, rn) - pf.torque;
          }
          /* i-side: horizontal reduce into fi */
          for (int c = 0; c < 3; ++c) {
            fi.f[c] += Kokkos::Experimental::reduce(fij.f[c], std::plus<>{});
            fi.torque[c] +=
                Kokkos::Experimental::reduce(fij.torque[c], std::plus<>{});
          }
          /* j-side: per-lane scatter */
          for (std::size_t l = 0; l < w; ++l) {
            auto const jl = (j + l) - offset;
            for (int c = 0; c < 3; ++c) {
              force_access(jl, c) += fji.f[c][l];
              torque_access(jl, c) += fji.torque[c][l];
            }
          }
        }
        /* (c) scalar tail over remaining j in (.., offset + n_local) */
        for (; j < offset + n_local; ++j) {
          Utils::Vector3d const pos_j{views.pos(j, 0), views.pos(j, 1),
                                      views.pos(j, 2)};
          Utils::Vector3d const m_j{views.m(j, 0), views.m(j, 1),
                                    views.m(j, 2)};
          auto const d0 = with_replicas ? (pos_i - pos_j)
                                        : box_geo.get_mi_vector(pos_i, pos_j);
          auto const jl = j - offset;
          for (std::size_t s = 0; s < n_shifts; ++s) {
            auto const rn = d0 + shifts_ptr[s];
            auto const pf = pair_force<double>(rn, m_i, m_j);
            fi.f += pf.f;
            fi.torque += pf.torque;
            for (int c = 0; c < 3; ++c) {
              force_access(jl, c) -= pf.f[c];
              torque_access(jl, c) += (vector_product(pf.f, rn) - pf.torque)[c];
            }
          }
        }
        /* (d) write i's own total directly (unique owner, no race) */
        local_particles_ptr[i]->force() += prefactor_local * fi.f;
        local_particles_ptr[i]->torque() += prefactor_local * fi.torque;
      });
  Kokkos::fence();

  /* Wait for remote data, fill remote slices. */
  boost::mpi::wait_all(reqs.begin(), reqs.end());
  fill_posmom_views(views, all_posmom, 0, offset);
  fill_posmom_views(views, all_posmom, offset + n_local, n_total);

  /* Phase B: remote pairs (red [0, offset) + black [offset + n_local,
   * n_total)), visit-twice, no scatter — accumulate only i. */
  Kokkos::parallel_for(
      "dds_remote_pairs",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i) {
        auto const gi = offset + i;
        Utils::Vector3d const pos_i{views.pos(gi, 0), views.pos(gi, 1),
                                    views.pos(gi, 2)};
        Utils::Vector3d const m_i{views.m(gi, 0), views.m(gi, 1),
                                  views.m(gi, 2)};
        auto const m_i_s = broadcast_simd(m_i);
        PairForce<double> fi{};

        /* Two remote ranges: red [0, offset) and black [offset + n_local,
         * n_total). Mirror the (b) SIMD body + (c) scalar tail of Phase A,
         * but without any scatter (each remote pair is visited twice, once
         * per owning rank, so only i accumulates). */
        std::size_t const ranges[2][2] = {{std::size_t{0}, offset},
                                          {offset + n_local, n_total}};
        for (auto const &range : ranges) {
          auto const range_begin = range[0];
          auto const range_end = range[1];
          auto j = range_begin;
          /* SIMD body */
          for (; j + w <= range_end; j += w) {
            auto const m_j = load_simd_moment(views, j);
            auto const d0 =
                primary_distance_simd(pos_i, views, j, with_replicas, box_geo);
            PairForce<simd_double> fij{};
            for (std::size_t s = 0; s < n_shifts; ++s) {
              auto const &shift = shifts_ptr[s];
              Utils::Vector<simd_double, 3> const rn{
                  d0[0] + shift[0], d0[1] + shift[1], d0[2] + shift[2]};
              fij += pair_force<simd_double>(rn, m_i_s, m_j);
            }
            for (int c = 0; c < 3; ++c) {
              fi.f[c] += Kokkos::Experimental::reduce(fij.f[c], std::plus<>{});
              fi.torque[c] +=
                  Kokkos::Experimental::reduce(fij.torque[c], std::plus<>{});
            }
          }
          /* scalar tail */
          for (; j < range_end; ++j) {
            Utils::Vector3d const pos_j{views.pos(j, 0), views.pos(j, 1),
                                        views.pos(j, 2)};
            Utils::Vector3d const m_j{views.m(j, 0), views.m(j, 1),
                                      views.m(j, 2)};
            auto const d0 = with_replicas ? (pos_i - pos_j)
                                          : box_geo.get_mi_vector(pos_i, pos_j);
            for (std::size_t s = 0; s < n_shifts; ++s) {
              auto const rn = d0 + shifts_ptr[s];
              auto const pf = pair_force<double>(rn, m_i, m_j);
              fi.f += pf.f;
              fi.torque += pf.torque;
            }
          }
        }
        local_particles_ptr[i]->force() += prefactor_local * fi.f;
        local_particles_ptr[i]->torque() += prefactor_local * fi.torque;
      });
  Kokkos::fence();

  /* Reduce the Newton's-third-law contributions and add to particles. */
  Kokkos::Experimental::contribute(local_force, scatter_force);
  Kokkos::Experimental::contribute(local_torque, scatter_torque);
  Kokkos::parallel_for(
      "dds_reduction",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i) {
        local_particles_ptr[i]->force() +=
            prefactor_local * Utils::Vector3d{local_force(i, 0),
                                              local_force(i, 1),
                                              local_force(i, 2)};
        local_particles_ptr[i]->torque() +=
            prefactor_local * Utils::Vector3d{local_torque(i, 0),
                                              local_torque(i, 1),
                                              local_torque(i, 2)};
      });
  Kokkos::fence();

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (not m_is_gpu) {
    dipole_field_at_part_cpu();
  }
#endif
}

/**
 * @brief Calculate the interaction potential.
 *
 * This employs a parallel N-square loop over all particle pairs.
 */
double DipolarDirectSum::long_range_energy_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &box_l = box_geo.length();
  auto const particles = system.cell_structure->local_particles();
  auto [local_particles, all_posmom, reqs, offset_signed] =
      gather_particle_data(box_geo, particles);

  /* Number of image boxes considered */
  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_l);

  auto const offset = static_cast<std::size_t>(offset_signed);
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();
  constexpr std::size_t w = simd_double::size();

  /* SoA views; fill the local slice only so the remote gather can overlap. */
  auto views = make_posmom_views(n_total);
  fill_posmom_views(views, all_posmom, offset, offset + n_local);

  using execution_space = Kokkos::DefaultExecutionSpace;

  /* Raw pointers so the Kokkos lambdas do not capture std::vector by value. */
  auto const *shifts_ptr = shifts.data();
  auto const n_shifts = shifts.size();

  /* Phase A: local-upper triangular sum over j in [gi, offset + n_local),
   * i.e. the self-image energy (shifts[1..], primary excluded) plus the pairs
   * with j in (gi, offset + n_local). Computed from the local SoA slice while
   * the remote data is still in flight. */
  double uA = 0.;
  Kokkos::parallel_reduce(
      "dds_energy_local",
      Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(
          std::size_t{0}, n_local)
          .set_chunk_size(64),
      [=](std::size_t const i, double &u_local) {
        auto const gi = offset + i;
        Utils::Vector3d const pos_i{views.pos(gi, 0), views.pos(gi, 1),
                                    views.pos(gi, 2)};
        Utils::Vector3d const m_i{views.m(gi, 0), views.m(gi, 1),
                                  views.m(gi, 2)};

        /* (a) self-images (shifts[1..], primary excluded) */
        for (std::size_t s = 1; s < n_shifts; ++s)
          u_local += pair_potential<double>(shifts_ptr[s], m_i, m_i);

        auto const m_i_s = broadcast_simd(m_i);

        /* (b) SIMD body over j in (gi, offset + n_local) */
        auto j = gi + 1;
        for (; j + w <= offset + n_local; j += w) {
          auto const m_j = load_simd_moment(views, j);
          auto const d0 =
              primary_distance_simd(pos_i, views, j, with_replicas, box_geo);
          simd_double acc{0.};
          for (std::size_t s = 0; s < n_shifts; ++s) {
            auto const &shift = shifts_ptr[s];
            Utils::Vector<simd_double, 3> const rn{
                d0[0] + shift[0], d0[1] + shift[1], d0[2] + shift[2]};
            acc += pair_potential<simd_double>(rn, m_i_s, m_j);
          }
          u_local += Kokkos::Experimental::reduce(acc, std::plus<>{});
        }
        /* (c) scalar tail over remaining j in (.., offset + n_local) */
        for (; j < offset + n_local; ++j) {
          Utils::Vector3d const pos_j{views.pos(j, 0), views.pos(j, 1),
                                      views.pos(j, 2)};
          Utils::Vector3d const m_j{views.m(j, 0), views.m(j, 1),
                                    views.m(j, 2)};
          auto const d0 = with_replicas ? (pos_i - pos_j)
                                        : box_geo.get_mi_vector(pos_i, pos_j);
          for (std::size_t s = 0; s < n_shifts; ++s)
            u_local += pair_potential<double>(d0 + shifts_ptr[s], m_i, m_j);
        }
      },
      uA);

  /* Wait for remote data, fill the black slice. The red range [0, offset) is
   * never summed by the energy kernel (each pair is counted once on the rank
   * owning its lower index). */
  boost::mpi::wait_all(reqs.begin(), reqs.end());
  fill_posmom_views(views, all_posmom, offset + n_local, n_total);

  /* Phase B: remote-black sum over j in [offset + n_local, n_total). No self
   * term and no primary exclusion — the range is entirely remote. */
  double uB = 0.;
  Kokkos::parallel_reduce(
      "dds_energy_remote",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i, double &u_local) {
        auto const gi = offset + i;
        Utils::Vector3d const pos_i{views.pos(gi, 0), views.pos(gi, 1),
                                    views.pos(gi, 2)};
        Utils::Vector3d const m_i{views.m(gi, 0), views.m(gi, 1),
                                  views.m(gi, 2)};
        auto const m_i_s = broadcast_simd(m_i);

        /* SIMD body over j in [offset + n_local, n_total) */
        auto j = offset + n_local;
        for (; j + w <= n_total; j += w) {
          auto const m_j = load_simd_moment(views, j);
          auto const d0 =
              primary_distance_simd(pos_i, views, j, with_replicas, box_geo);
          simd_double acc{0.};
          for (std::size_t s = 0; s < n_shifts; ++s) {
            auto const &shift = shifts_ptr[s];
            Utils::Vector<simd_double, 3> const rn{
                d0[0] + shift[0], d0[1] + shift[1], d0[2] + shift[2]};
            acc += pair_potential<simd_double>(rn, m_i_s, m_j);
          }
          u_local += Kokkos::Experimental::reduce(acc, std::plus<>{});
        }
        /* scalar tail */
        for (; j < n_total; ++j) {
          Utils::Vector3d const pos_j{views.pos(j, 0), views.pos(j, 1),
                                      views.pos(j, 2)};
          Utils::Vector3d const m_j{views.m(j, 0), views.m(j, 1),
                                    views.m(j, 2)};
          auto const d0 = with_replicas ? (pos_i - pos_j)
                                        : box_geo.get_mi_vector(pos_i, pos_j);
          for (std::size_t s = 0; s < n_shifts; ++s)
            u_local += pair_potential<double>(d0 + shifts_ptr[s], m_i, m_j);
        }
      },
      uB);

  return prefactor * (uA + uB);
}

/**
 * @brief Calculate total dipole field of each particle.
 *
 * This employs a parallel N-square loop over all particles.
 * Logically this is equivalent to the potential calculation
 * in @ref DipolarDirectSum::long_range_energy, which calculates
 * a naive N-square sum. The difference is summation range,
 * and the kernel calculates the dipole field rather than the energy.
 */
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
void DipolarDirectSum::dipole_field_at_part_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &box_l = box_geo.length();
  auto const particles = system.cell_structure->local_particles();
  /* collect particle data */
  auto [local_particles, all_posmom, reqs, offset_signed] =
      gather_particle_data(box_geo, particles);

  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_l);

  auto const offset = static_cast<std::size_t>(offset_signed);
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();
  constexpr std::size_t w = simd_double::size();

  auto const prefactor_local = prefactor;

  /* The field sweeps over all j, so every view slice is needed. Unlike the
   * force/energy kernels there is no local-only computation to overlap with,
   * so wait for the remote data first, then fill all slices. */
  boost::mpi::wait_all(reqs.begin(), reqs.end());
  auto views = make_posmom_views(n_total);
  fill_posmom_views(views, all_posmom, 0, n_total);

  using execution_space = Kokkos::DefaultExecutionSpace;

  /* Raw pointers so the Kokkos lambdas do not capture std::vector by value. */
  auto *local_particles_ptr = local_particles.data();
  auto const *shifts_ptr = shifts.data();
  auto const n_shifts = shifts.size();

  Kokkos::parallel_for(
      "dds_dipole_field",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i) {
        auto const gi = offset + i;
        Utils::Vector3d const pos_i{views.pos(gi, 0), views.pos(gi, 1),
                                    views.pos(gi, 2)};
        Utils::Vector3d const m_i{views.m(gi, 0), views.m(gi, 1),
                                  views.m(gi, 2)};
        Utils::Vector3d u{};

        /* (a) self-image term over shifts[1..] (primary excluded) */
        for (std::size_t s = 1; s < n_shifts; ++s)
          u += dipole_field<double>(shifts_ptr[s], m_i);

        /* Sweep over all j in [0, n_total), self-primary excluded by splitting
         * the range into [0, gi) and [gi + 1, n_total). Each sub-range is a
         * SIMD body plus a scalar tail. */
        std::size_t const ranges[2][2] = {{std::size_t{0}, gi},
                                          {gi + 1, n_total}};
        for (auto const &range : ranges) {
          auto const range_begin = range[0];
          auto const range_end = range[1];
          auto j = range_begin;
          /* SIMD body */
          for (; j + w <= range_end; j += w) {
            auto const m_j = load_simd_moment(views, j);
            auto const d0 =
                primary_distance_simd(pos_i, views, j, with_replicas, box_geo);
            Utils::Vector<simd_double, 3> acc{};
            for (std::size_t s = 0; s < n_shifts; ++s) {
              auto const &shift = shifts_ptr[s];
              Utils::Vector<simd_double, 3> const rn{
                  d0[0] + shift[0], d0[1] + shift[1], d0[2] + shift[2]};
              acc += dipole_field<simd_double>(rn, m_j);
            }
            for (int c = 0; c < 3; ++c)
              u[c] += Kokkos::Experimental::reduce(acc[c], std::plus<>{});
          }
          /* scalar tail */
          for (; j < range_end; ++j) {
            Utils::Vector3d const pos_j{views.pos(j, 0), views.pos(j, 1),
                                        views.pos(j, 2)};
            Utils::Vector3d const m_j{views.m(j, 0), views.m(j, 1),
                                      views.m(j, 2)};
            auto const d0 = with_replicas ? (pos_i - pos_j)
                                          : box_geo.get_mi_vector(pos_i, pos_j);
            for (std::size_t s = 0; s < n_shifts; ++s)
              u += dipole_field<double>(d0 + shifts_ptr[s], m_j);
          }
        }
        local_particles_ptr[i]->dip_fld() = prefactor_local * u;
      });
  Kokkos::fence();
}
#endif

DipolarDirectSum::DipolarDirectSum(double prefactor, int n_replicas, bool gpu) {
  set_prefactor(prefactor);
  m_is_gpu = gpu;
  this->n_replicas = n_replicas;
  if (n_replicas < 0) {
    throw std::domain_error("Parameter 'n_replicas' must be >= 0");
  }
}

#endif // ESPRESSO_DIPOLES
