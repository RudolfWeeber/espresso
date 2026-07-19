/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_KOKKOS_FFT

#include "p3m/P3MFFTBackend.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <KokkosFFT.hpp>
#include <Kokkos_Core.hpp>

#include <cassert>
#include <complex>
#include <memory>

/**
 * @brief Single-MPI-rank P3M FFT backend built on kokkos-fft.
 *
 * kokkos-fft wraps the host FFT library (FFTW) behind a @c Kokkos::View API
 * and performs purely local transforms, so it is only selected when the
 * simulation runs on a single MPI rank; multi-rank runs keep the heFFTe
 * backend (@ref P3MFFTHeffte). It reproduces heFFTe's row-major r2c layout
 * (reduced last axis, size @c mesh[2]/2+1) and its unscaled convention
 * (@c scale::none in both directions; P3M folds the @c 1/N into the influence
 * function), so forces and energies agree with the heFFTe path to
 * floating-point round-off. The transform runs on the host regardless of the
 * Kokkos default device, matching the CPU heFFTe backend it replaces.
 */
template <typename FloatType, class FFTConfig>
struct P3MFFTKokkos final : public P3MFFTBackend<FloatType, FFTConfig> {
  static_assert(FFTConfig::use_r2c,
                "kokkos-fft P3M backend implements the r2c transform only");
  static_assert(FFTConfig::r2c_dir == 2u,
                "kokkos-fft P3M backend reduces the contiguous last axis");
  static_assert(FFTConfig::r_space_order == Utils::MemoryOrder::ROW_MAJOR and
                    FFTConfig::k_space_order == Utils::MemoryOrder::ROW_MAJOR,
                "kokkos-fft P3M backend supports row-major layout only");

  using Base = P3MFFTBackend<FloatType, FFTConfig>;
  using ComplexType = typename Base::ComplexType;
  using RSpaceScalar = typename Base::RSpaceScalar;
  using KComplex = Kokkos::complex<FloatType>;
  using ExecSpace = Kokkos::DefaultHostExecutionSpace;

  using RealView =
      Kokkos::View<FloatType ***, Kokkos::LayoutRight, Kokkos::HostSpace>;
  using CplxView =
      Kokkos::View<KComplex ***, Kokkos::LayoutRight, Kokkos::HostSpace>;
  using RealViewU =
      Kokkos::View<FloatType ***, Kokkos::LayoutRight, Kokkos::HostSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using CplxViewU =
      Kokkos::View<KComplex ***, Kokkos::LayoutRight, Kokkos::HostSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using ForwardPlan = KokkosFFT::Plan<ExecSpace, RealView, CplxView, 3>;
  using BackwardPlan = KokkosFFT::Plan<ExecSpace, CplxView, RealView, 3>;

  P3MFFTKokkos(boost::mpi::communicator comm,
               Utils::Vector3i const &global_mesh,
               Utils::Vector3i const &rs_local_ld_index,
               Utils::Vector3i const &rs_local_ur_index,
               Utils::Vector3i const & /* node_grid */)
      : m_mesh{global_mesh},
        m_ks_size{global_mesh[0], global_mesh[1], global_mesh[2] / 2 + 1} {
    // kokkos-fft is local-only; single rank means the no-halo local mesh spans
    // the whole global mesh, so there is no domain decomposition to mirror.
    assert(comm.size() == 1);
    assert((rs_local_ur_index - rs_local_ld_index) == global_mesh);
    static_cast<void>(comm);
    static_cast<void>(rs_local_ld_index);
    static_cast<void>(rs_local_ur_index);

    m_real = RealView(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "P3MFFTKokkos::real"),
        m_mesh[0], m_mesh[1], m_mesh[2]);
    m_cplx = CplxView(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "P3MFFTKokkos::cplx"),
        m_ks_size[0], m_ks_size[1], m_ks_size[2]);
    auto const axes = KokkosFFT::axis_type<3>({0, 1, 2});
    m_forward = std::make_unique<ForwardPlan>(
        ExecSpace{}, m_real, m_cplx, KokkosFFT::Direction::forward, axes);
    m_backward = std::make_unique<BackwardPlan>(
        ExecSpace{}, m_cplx, m_real, KokkosFFT::Direction::backward, axes);
  }

  Utils::Vector3i ks_local_ld_index() const override { return {0, 0, 0}; }
  Utils::Vector3i ks_local_ur_index() const override { return m_ks_size; }
  Utils::Vector3i ks_local_size() const override { return m_ks_size; }
  Utils::Vector3i rs_local_size() const override { return m_mesh; }

  void forward(RSpaceScalar const *in, ComplexType *out) override {
    RealViewU in_view(const_cast<FloatType *>(in), m_mesh[0], m_mesh[1],
                      m_mesh[2]);
    Kokkos::deep_copy(m_real, in_view);
    KokkosFFT::execute(*m_forward, m_real, m_cplx,
                       KokkosFFT::Normalization::none);
    CplxViewU out_view(reinterpret_cast<KComplex *>(out), m_ks_size[0],
                       m_ks_size[1], m_ks_size[2]);
    Kokkos::deep_copy(out_view, m_cplx);
  }

  void backward(ComplexType const *in, RSpaceScalar *out) override {
    CplxViewU in_view(
        reinterpret_cast<KComplex *>(const_cast<ComplexType *>(in)),
        m_ks_size[0], m_ks_size[1], m_ks_size[2]);
    Kokkos::deep_copy(m_cplx, in_view);
    KokkosFFT::execute(*m_backward, m_cplx, m_real,
                       KokkosFFT::Normalization::none);
    RealViewU out_view(out, m_mesh[0], m_mesh[1], m_mesh[2]);
    Kokkos::deep_copy(out_view, m_real);
  }

private:
  Utils::Vector3i m_mesh;
  Utils::Vector3i m_ks_size;
  RealView m_real;
  CplxView m_cplx;
  std::unique_ptr<ForwardPlan> m_forward;
  std::unique_ptr<BackwardPlan> m_backward;
};

#endif // ESPRESSO_KOKKOS_FFT
