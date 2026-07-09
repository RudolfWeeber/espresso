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
/** \file
 *  Ghost particles and particle exchange.
 *
 *  For more information on ghosts,
 *  see \ref ghosts.hpp "ghosts.hpp"
 *
 * Note on variable naming:
 * - a "GhostCommunicator" is always named "gcr",
 * - a "GhostCommunication" is always named "ghost_comm".
 */
#include "ghosts.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "particle_store/ParticleStore.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>
#include <utils/serialization/memcpy_archive.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/range/numeric.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <span>
#include <vector>

/** Tag for ghosts communications. */
#define REQ_GHOST_SEND 100

/**
 * Class that stores marshalled data for ghost communications.
 * To store and retrieve data, use the adapter classes below.
 */
class CommBuf {
public:
  /** Returns a pointer to the non-bond storage.
   */
  char *data() { return buf.data(); }
  const char *data() const { return buf.data(); }

  /** Returns the number of elements in the non-bond storage.
   */
  std::size_t size() const { return buf.size(); }

  /** Resizes the underlying storage s.t. the object is capable
   * of holding "new_size" chars.
   * @param new_size new size
   */
  void resize(std::size_t new_size) { buf.resize(new_size); }

  /** Returns a reference to the bond storage.
   */
  auto &bonds() { return bondbuf; }
  const auto &bonds() const { return bondbuf; }

  auto make_span() { return std::span(buf.data(), buf.size()); }

private:
  std::vector<char> buf;     ///< Buffer for everything but bonds
  std::vector<char> bondbuf; ///< Buffer for bond lists
};

/** @brief Pseudo-archive to calculate the size of the serialization buffer. */
class SerializationSizeCalculator {
  std::size_t m_size = 0;

public:
  auto size() const { return m_size; }

  template <class T> auto &operator<<(T &) {
    m_size += sizeof(T);
    return *this;
  }

  template <class T> auto &operator&(T &t) { return *this << t; }
};

/**
 * @brief Raw column pointers + strides for the mixed-parts POSITION fast path.
 *
 * Captured once per communication from the (clean) ParticleStore's host views
 * so the generic per-particle POSITION serialization can read/write the four
 * position columns directly by row instead of through per-field accessor
 * proxies. See make_position_row_context() below for how it is built and for
 * the applicability rationale. A row-r component-c element of a column lives at
 * base + r*row_stride + c*comp_stride (comp_stride==1 and one row contiguous
 * under the store's LayoutRight).
 */
struct PositionRowContext {
  double *pos; ///< base of the position column (row 0, component 0)
  int *img;    ///< base of the image-box column
  std::size_t pos_row_stride;
  std::size_t pos_comp_stride;
  std::size_t img_row_stride;
  std::size_t img_comp_stride;
#ifdef ESPRESSO_ROTATION
  double *quat;
  std::size_t quat_row_stride;
  std::size_t quat_comp_stride;
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  double *plast;
  std::size_t plast_row_stride;
  std::size_t plast_comp_stride;
#endif
};

/** @brief Raw-pointer view of the force/torque columns, hoisted once per
 *  communication so the per-part_list bulk loops pay no view-handle
 *  (refcount) or DualView bookkeeping cost per list. */
struct ForceRowContext {
  double *force;
  std::size_t force_row_stride;
  std::size_t force_comp_stride;
#ifdef ESPRESSO_ROTATION
  double *torque;
  std::size_t torque_row_stride;
  std::size_t torque_comp_stride;
#endif
};

/**
 * @brief Raw column pointers + strides for the mixed-parts MOMENTUM fast path.
 *
 * Mirror of PositionRowContext for the velocity/angular-velocity columns.
 * Captured once per communication from the (clean) ParticleStore's host views
 * so the generic per-particle MOMENTUM serialization can read/write the
 * velocity (and, under ROTATION, angular-velocity) columns directly by row
 * instead of through the per-field accessor proxies. See
 * make_momentum_row_context() for how it is built and for the applicability
 * rationale. A row-r component-c element lives at
 * base + r*row_stride + c*comp_stride (comp_stride==1 and one row contiguous
 * under the store's LayoutRight).
 */
struct MomentumRowContext {
  double *vel; ///< base of the velocity column (row 0, component 0)
  std::size_t vel_row_stride;
  std::size_t vel_comp_stride;
#ifdef ESPRESSO_ROTATION
  double *omega; ///< base of the angular-velocity column
  std::size_t omega_row_stride;
  std::size_t omega_comp_stride;
#endif
};

/**
 * @brief Raw column pointers + strides (and a store handle for the cold PODs)
 *        for the mixed-parts PROPERTIES value path.
 *
 * Mirror of MomentumRowContext, but for the parameter group. Captured once per
 * communication from the (clean) ParticleStore so the generic per-particle
 * PROPERTIES serialization can read/write the hot parameter columns directly by
 * row instead of through the per-field accessor proxies, and read/write the
 * three cold parameter PODs through the store's host sidecars by row.
 *
 * The scalar columns (id/mol_id/type/propagation, and the ifdef-gated uint8
 * bitfields rotation/ext_flag and doubles mass/q/dipm) are 1-D DualViews: a
 * row-r element lives at base + r (stride 1). The Vector3d columns
 * (rinertia/mu_E/gamma/gamma_rot/ext_force/ext_torque) are 2-D LayoutRight
 * columns, so a row-r component-c element lives at base + r*row_stride +
 * c*comp_stride (comp_stride==1 and one row contiguous). The gamma/gamma_rot
 * columns are scalar (double) when ESPRESSO_PARTICLE_ANISOTROPY is off and 2-D
 * Vector3d columns when it is on; the pointer + optional strides follow the
 * GammaColumn typedef switch, exactly like the accessor.
 *
 * The cold PODs (swim/magnetodynamics/vs_relative) are NOT Kokkos columns; they
 * live in host std::vector sidecars indexed by store row. The context therefore
 * carries the ParticleStore* and reads/writes them via the by-row sidecar
 * accessors (store->swimming(row) etc.).
 *
 * Wired like MomentumRowContext: built per communication (only when the store
 * is present and clean, WITHOUT the data_parts==POSITION restriction of
 * columnar_eligible), passed to serialize_and_reduce as an optional pointer
 * (nullptr = accessor fallback), and consulted only when p.store_row() >= 0.
 */
struct ParameterRowContext {
  ParticleStore *store; ///< sidecar access for the cold PODs, by row
  int *id;
  int *mol_id;
  int *type;
  int *propagation;
#ifdef ESPRESSO_ROTATION
  std::uint8_t *rotation;
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  double *rinertia;
  std::size_t rinertia_row_stride;
  std::size_t rinertia_comp_stride;
#endif
#endif
#ifdef ESPRESSO_MASS
  double *mass;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  double *q;
#endif
#ifdef ESPRESSO_DIPOLES
  double *dipm;
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  double *mu_E;
  std::size_t mu_E_row_stride;
  std::size_t mu_E_comp_stride;
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  double *gamma;
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  std::size_t gamma_row_stride;
  std::size_t gamma_comp_stride;
#endif
#ifdef ESPRESSO_ROTATION
  double *gamma_rot;
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  std::size_t gamma_rot_row_stride;
  std::size_t gamma_rot_comp_stride;
#endif
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  std::uint8_t *ext_flag;
  double *ext_force;
  std::size_t ext_force_row_stride;
  std::size_t ext_force_comp_stride;
#ifdef ESPRESSO_ROTATION
  double *ext_torque;
  std::size_t ext_torque_row_stride;
  std::size_t ext_torque_comp_stride;
#endif
#endif
};

/* Forward declaration: used by debug assertions inside serialize_and_reduce
 * (defined below with the columnar machinery). */
static ParticleStore *active_particle_store();

/** @brief Type of reduction to carry out during serialization. */
enum class ReductionPolicy {
  /** @brief Reduction for domain-to-domain particle communication. */
  MOVE,
  /** @brief Reduction for cell-to-cell particle update. */
  UPDATE,
};

/** @brief Whether to save the state to or load the state from the archive. */
enum class SerializationDirection { SAVE, LOAD };

/**
 * @brief Serialize particle data, possibly with reduction.
 * The reduction can take place during the save stage, e.g. to apply
 * a ghost shift to the particle position, or during the load stage,
 * e.g. to transfer momentum between particles in two local cells.
 */
template <class Archive>
static void
serialize_and_reduce(Archive &ar, Particle &p, unsigned int data_parts,
                     ReductionPolicy policy, SerializationDirection direction,
                     BoxGeometry const &box_geo,
                     Utils::Vector3d const *ghost_shift,
                     PositionRowContext const *pos_ctx = nullptr,
                     MomentumRowContext const *mom_ctx = nullptr,
                     ParameterRowContext const *param_ctx = nullptr) {
  if (data_parts & GHOSTTRANS_PROPRTS) {
    /* Properties have no reduction policy: like POSITION/MOMENTUM this group is
     * MOVE semantics for BOTH ReductionPolicy values. Enumerating all four
     * (policy, direction) combinations of the previous policy-blind
     * `ar & p.id() & ...`:
     *   (MOVE,  SAVE): write each field to the archive from the particle.
     *   (UPDATE,SAVE): identical -- policy is not consulted.
     *   (MOVE,  LOAD): read each field FROM the archive, overwriting.
     *   (UPDATE,LOAD): identical -- overwrite, NO accumulation (no `+=`).
     * So on LOAD we always write the received value INTO the particle
     * regardless of policy; on SAVE we always read the value FROM the particle.
     * Never bind field-storage memory to the archive directly (it may be a
     * proxy, which cannot be serialized): local copies are made for/from the
     * archive. The wire layout is byte-packed (fixed field order, fixed ifdef
     * structure, tightly packed field sizes: the three PODs are
     * bitwise-serializable, so each is packed whole as sizeof(T)).
     *
     * When a ParameterRowContext is present AND the particle is attached
     * (store_row() >= 0), read/write the hot parameter columns and the cold POD
     * sidecars directly by row rather than through the per-field accessor
     * proxies -- this is the mixed-parts hot path. An unattached particle
     * (resort-step fresh ghost) or an absent context falls back to the
     * accessor. */
    auto const row = p.store_row();
    bool const use_ctx = param_ctx != nullptr and row >= 0;
    auto const urow = static_cast<std::size_t>(row);
#ifndef NDEBUG
    // When use_ctx is true the particle must belong to the active store: a
    // stale param_ctx from a different store generation or a different store
    // object would silently corrupt the wrong memory region.
    if (use_ctx) {
      assert(p.store() == active_particle_store() and
             "serialize_and_reduce: use_ctx row used but p.store() != active "
             "store");
    }
#endif
    if (direction == SerializationDirection::LOAD) {
      int id = 0, mol_id = 0, type = 0, propagation = 0;
      ar & id & mol_id & type & propagation;
      if (use_ctx) {
        param_ctx->id[urow] = id;
        param_ctx->mol_id[urow] = mol_id;
        param_ctx->type[urow] = type;
        param_ctx->propagation[urow] = propagation;
      } else {
        p.id() = id;
        p.mol_id() = mol_id;
        p.type() = type;
        p.propagation() = propagation;
      }
    } else {
      int id = 0, mol_id = 0, type = 0, propagation = 0;
      if (use_ctx) {
        id = param_ctx->id[urow];
        mol_id = param_ctx->mol_id[urow];
        type = param_ctx->type[urow];
        propagation = param_ctx->propagation[urow];
      } else {
        id = p.id();
        mol_id = p.mol_id();
        type = p.type();
        propagation = p.propagation();
      }
      ar & id & mol_id & type & propagation;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      std::uint8_t rotation = 0u;
      ar & rotation;
      if (use_ctx) {
        param_ctx->rotation[urow] = rotation;
      } else {
        p.rotation() = rotation;
      }
    } else {
      std::uint8_t rotation =
          use_ctx ? param_ctx->rotation[urow] : p.rotation();
      ar & rotation;
    }
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d rinertia{};
      ar & rinertia;
      if (use_ctx) {
        auto *r = param_ctx->rinertia + urow * param_ctx->rinertia_row_stride;
        r[0] = rinertia[0];
        r[param_ctx->rinertia_comp_stride] = rinertia[1];
        r[2u * param_ctx->rinertia_comp_stride] = rinertia[2];
      } else {
        p.rinertia() = rinertia;
      }
    } else {
      Utils::Vector3d rinertia{};
      if (use_ctx) {
        auto const *r =
            param_ctx->rinertia + urow * param_ctx->rinertia_row_stride;
        rinertia = {r[0], r[param_ctx->rinertia_comp_stride],
                    r[2u * param_ctx->rinertia_comp_stride]};
      } else {
        rinertia = p.rinertia();
      }
      ar & rinertia;
    }
#endif
#endif
#ifdef ESPRESSO_MASS
    if (direction == SerializationDirection::LOAD) {
      double mass = 0.;
      ar & mass;
      if (use_ctx) {
        param_ctx->mass[urow] = mass;
      } else {
        p.mass() = mass;
      }
    } else {
      double mass = use_ctx ? param_ctx->mass[urow] : p.mass();
      ar & mass;
    }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    if (direction == SerializationDirection::LOAD) {
      double q = 0.;
      ar & q;
      if (use_ctx) {
        param_ctx->q[urow] = q;
      } else {
        p.q() = q;
      }
    } else {
      double q = use_ctx ? param_ctx->q[urow] : p.q();
      ar & q;
    }
#endif
#ifdef ESPRESSO_DIPOLES
    if (direction == SerializationDirection::LOAD) {
      double dipm = 0.;
      ar & dipm;
      if (use_ctx) {
        param_ctx->dipm[urow] = dipm;
      } else {
        p.dipm() = dipm;
      }
    } else {
      double dipm = use_ctx ? param_ctx->dipm[urow] : p.dipm();
      ar & dipm;
    }
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d mu_E;
      ar & mu_E;
      if (use_ctx) {
        auto *m = param_ctx->mu_E + urow * param_ctx->mu_E_row_stride;
        m[0] = mu_E[0];
        m[param_ctx->mu_E_comp_stride] = mu_E[1];
        m[2u * param_ctx->mu_E_comp_stride] = mu_E[2];
      } else {
        p.mu_E() = mu_E;
      }
    } else {
      Utils::Vector3d mu_E;
      if (use_ctx) {
        auto const *m = param_ctx->mu_E + urow * param_ctx->mu_E_row_stride;
        mu_E = {m[0], m[param_ctx->mu_E_comp_stride],
                m[2u * param_ctx->mu_E_comp_stride]};
      } else {
        mu_E = p.mu_E();
      }
      ar & mu_E;
    }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    if (direction == SerializationDirection::LOAD) {
      VirtualSitesRelativeParameters vs_relative;
      ar & vs_relative;
      if (use_ctx) {
        param_ctx->store->vs_relative(row) = vs_relative;
      } else {
        p.vs_relative() = vs_relative;
      }
    } else {
      VirtualSitesRelativeParameters vs_relative =
          use_ctx ? param_ctx->store->vs_relative(row) : p.vs_relative();
      ar & vs_relative;
    }
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    if (direction == SerializationDirection::LOAD) {
      ParticleStore::GammaValue gamma;
      ar & gamma;
      if (use_ctx) {
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
        auto *g = param_ctx->gamma + urow * param_ctx->gamma_row_stride;
        g[0] = gamma[0];
        g[param_ctx->gamma_comp_stride] = gamma[1];
        g[2u * param_ctx->gamma_comp_stride] = gamma[2];
#else
        param_ctx->gamma[urow] = gamma;
#endif
      } else {
        p.gamma() = gamma;
      }
    } else {
      ParticleStore::GammaValue gamma;
      if (use_ctx) {
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
        auto const *g = param_ctx->gamma + urow * param_ctx->gamma_row_stride;
        gamma = {g[0], g[param_ctx->gamma_comp_stride],
                 g[2u * param_ctx->gamma_comp_stride]};
#else
        gamma = param_ctx->gamma[urow];
#endif
      } else {
        gamma = p.gamma();
      }
      ar & gamma;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      ParticleStore::GammaValue gamma_rot;
      ar & gamma_rot;
      if (use_ctx) {
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
        auto *g = param_ctx->gamma_rot + urow * param_ctx->gamma_rot_row_stride;
        g[0] = gamma_rot[0];
        g[param_ctx->gamma_rot_comp_stride] = gamma_rot[1];
        g[2u * param_ctx->gamma_rot_comp_stride] = gamma_rot[2];
#else
        param_ctx->gamma_rot[urow] = gamma_rot;
#endif
      } else {
        p.gamma_rot() = gamma_rot;
      }
    } else {
      ParticleStore::GammaValue gamma_rot;
      if (use_ctx) {
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
        auto const *g =
            param_ctx->gamma_rot + urow * param_ctx->gamma_rot_row_stride;
        gamma_rot = {g[0], g[param_ctx->gamma_rot_comp_stride],
                     g[2u * param_ctx->gamma_rot_comp_stride]};
#else
        gamma_rot = param_ctx->gamma_rot[urow];
#endif
      } else {
        gamma_rot = p.gamma_rot();
      }
      ar & gamma_rot;
    }
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    if (direction == SerializationDirection::LOAD) {
      std::uint8_t ext_flag = 0u;
      ar & ext_flag;
      if (use_ctx) {
        param_ctx->ext_flag[urow] = ext_flag;
      } else {
        p.fixed() = ext_flag;
      }
    } else {
      std::uint8_t ext_flag = use_ctx ? param_ctx->ext_flag[urow] : p.fixed();
      ar & ext_flag;
    }
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d ext_force{};
      ar & ext_force;
      if (use_ctx) {
        auto *f = param_ctx->ext_force + urow * param_ctx->ext_force_row_stride;
        f[0] = ext_force[0];
        f[param_ctx->ext_force_comp_stride] = ext_force[1];
        f[2u * param_ctx->ext_force_comp_stride] = ext_force[2];
      } else {
        p.ext_force() = ext_force;
      }
    } else {
      Utils::Vector3d ext_force{};
      if (use_ctx) {
        auto const *f =
            param_ctx->ext_force + urow * param_ctx->ext_force_row_stride;
        ext_force = {f[0], f[param_ctx->ext_force_comp_stride],
                     f[2u * param_ctx->ext_force_comp_stride]};
      } else {
        ext_force = p.ext_force();
      }
      ar & ext_force;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d ext_torque{};
      ar & ext_torque;
      if (use_ctx) {
        auto *t =
            param_ctx->ext_torque + urow * param_ctx->ext_torque_row_stride;
        t[0] = ext_torque[0];
        t[param_ctx->ext_torque_comp_stride] = ext_torque[1];
        t[2u * param_ctx->ext_torque_comp_stride] = ext_torque[2];
      } else {
        p.ext_torque() = ext_torque;
      }
    } else {
      Utils::Vector3d ext_torque{};
      if (use_ctx) {
        auto const *t =
            param_ctx->ext_torque + urow * param_ctx->ext_torque_row_stride;
        ext_torque = {t[0], t[param_ctx->ext_torque_comp_stride],
                      t[2u * param_ctx->ext_torque_comp_stride]};
      } else {
        ext_torque = p.ext_torque();
      }
      ar & ext_torque;
    }
#endif
#endif
#ifdef ESPRESSO_ENGINE
    if (direction == SerializationDirection::LOAD) {
      ParticleParametersSwimming swimming;
      ar & swimming;
      if (use_ctx) {
        param_ctx->store->swimming(row) = swimming;
      } else {
        p.swimming() = swimming;
      }
    } else {
      ParticleParametersSwimming swimming =
          use_ctx ? param_ctx->store->swimming(row) : p.swimming();
      ar & swimming;
    }
#endif
  }
  if (data_parts & GHOSTTRANS_POSITION) {
    /* Position has no reduction policy: it is MOVE semantics for both
     * ReductionPolicy values. Branch on direction first: on LOAD we always
     * write the received value INTO the particle; on SAVE we read the value
     * FROM the particle, optionally with the ghost shift + fold applied. Never
     * bind field-storage memory to the archive directly (it may be a proxy,
     * which cannot be serialized). Wire layout is byte-packed.
     *
     * When a PositionRowContext is present AND the particle is attached
     * (store_row() >= 0), read/write pos/image (and quat, pos_last_time_step)
     * directly through the raw column pointers rather than the per-field
     * accessor proxies -- this is the mixed-parts (POSITION|PROPERTIES) hot
     * path. An unattached particle (resort-step fresh ghost) or an absent
     * context falls back to the accessor. Local copies are made for the archive
     * (SAVE) / from the archive (LOAD), so the shift+fold and the assignment
     * semantics are preserved. */
    auto const row = p.store_row();
    bool const use_ctx = pos_ctx != nullptr and row >= 0;
#ifndef NDEBUG
    // When use_ctx is true the particle must belong to the active store: a
    // stale pos_ctx from a different store generation or a different store
    // object would silently corrupt the wrong memory region.
    if (use_ctx) {
      assert(p.store() == active_particle_store() and
             "serialize_and_reduce: use_ctx row used but p.store() != active "
             "store");
    }
#endif
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d position{};
      Utils::Vector3i image_box;
      ar & position;
      ar & image_box;
      if (use_ctx) {
        auto *pos = pos_ctx->pos +
                    static_cast<std::size_t>(row) * pos_ctx->pos_row_stride;
        auto *img = pos_ctx->img +
                    static_cast<std::size_t>(row) * pos_ctx->img_row_stride;
        pos[0 * pos_ctx->pos_comp_stride] = position[0];
        pos[1 * pos_ctx->pos_comp_stride] = position[1];
        pos[2 * pos_ctx->pos_comp_stride] = position[2];
        img[0 * pos_ctx->img_comp_stride] = image_box[0];
        img[1 * pos_ctx->img_comp_stride] = image_box[1];
        img[2 * pos_ctx->img_comp_stride] = image_box[2];
      } else {
        p.pos() = position;
        p.image_box() = image_box;
      }
    } else {
      Utils::Vector3d position{};
      Utils::Vector3i image_box;
      if (use_ctx) {
        auto const *pos = pos_ctx->pos + static_cast<std::size_t>(row) *
                                             pos_ctx->pos_row_stride;
        auto const *img = pos_ctx->img + static_cast<std::size_t>(row) *
                                             pos_ctx->img_row_stride;
        position = {pos[0 * pos_ctx->pos_comp_stride],
                    pos[1 * pos_ctx->pos_comp_stride],
                    pos[2 * pos_ctx->pos_comp_stride]};
        image_box = {img[0 * pos_ctx->img_comp_stride],
                     img[1 * pos_ctx->img_comp_stride],
                     img[2 * pos_ctx->img_comp_stride]};
      } else {
        position = p.pos();
        image_box = p.image_box();
      }
      if (ghost_shift != nullptr) {
        /* ok, this is not nice, but perhaps fast */
        position += *ghost_shift;
        box_geo.fold_position(position, image_box);
      }
      ar & position;
      ar & image_box;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      Utils::Quaternion<double> quat{};
      ar & quat;
      if (use_ctx) {
        auto *q = pos_ctx->quat +
                  static_cast<std::size_t>(row) * pos_ctx->quat_row_stride;
        q[0 * pos_ctx->quat_comp_stride] = quat[0];
        q[1 * pos_ctx->quat_comp_stride] = quat[1];
        q[2 * pos_ctx->quat_comp_stride] = quat[2];
        q[3 * pos_ctx->quat_comp_stride] = quat[3];
      } else {
        p.quat() = quat;
      }
    } else {
      Utils::Quaternion<double> quat{};
      if (use_ctx) {
        auto const *q = pos_ctx->quat + static_cast<std::size_t>(row) *
                                            pos_ctx->quat_row_stride;
        quat[0] = q[0 * pos_ctx->quat_comp_stride];
        quat[1] = q[1 * pos_ctx->quat_comp_stride];
        quat[2] = q[2 * pos_ctx->quat_comp_stride];
        quat[3] = q[3 * pos_ctx->quat_comp_stride];
      } else {
        quat = p.quat();
      }
      ar & quat;
    }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d pos_last_time_step{};
      ar & pos_last_time_step;
      if (use_ctx) {
        auto *pl = pos_ctx->plast +
                   static_cast<std::size_t>(row) * pos_ctx->plast_row_stride;
        pl[0 * pos_ctx->plast_comp_stride] = pos_last_time_step[0];
        pl[1 * pos_ctx->plast_comp_stride] = pos_last_time_step[1];
        pl[2 * pos_ctx->plast_comp_stride] = pos_last_time_step[2];
      } else {
        p.pos_last_time_step() = pos_last_time_step;
      }
    } else {
      Utils::Vector3d pos_last_time_step{};
      if (use_ctx) {
        auto const *pl = pos_ctx->plast + static_cast<std::size_t>(row) *
                                              pos_ctx->plast_row_stride;
        pos_last_time_step = {pl[0 * pos_ctx->plast_comp_stride],
                              pl[1 * pos_ctx->plast_comp_stride],
                              pl[2 * pos_ctx->plast_comp_stride]};
      } else {
        pos_last_time_step = p.pos_last_time_step();
      }
      ar & pos_last_time_step;
    }
#endif
  }
  if (data_parts & GHOSTTRANS_MOMENTUM) {
    /* Momentum (velocity, and angular velocity under ROTATION) has no reduction
     * policy: like POSITION it is MOVE semantics for BOTH ReductionPolicy
     * values. Enumerating all four (policy, direction) combinations of the
     * previous `ar & p.v(); ar & p.omega();`:
     *   (MOVE,  SAVE): write v (then omega) to the archive.
     *   (UPDATE,SAVE): identical -- `ar & p.v()` never consulted the policy.
     *   (MOVE,  LOAD): read v (then omega) FROM the archive, overwriting.
     *   (UPDATE,LOAD): identical -- overwrite, NO accumulation (no `+=`).
     * So on LOAD we always write the received value INTO the particle
     * regardless of policy; on SAVE we always read the value FROM the particle.
     * Never bind field-storage memory to the archive directly (it may be a
     * proxy, which cannot be serialized). Wire layout is byte-packed.
     *
     * When a MomentumRowContext is present AND the particle is attached
     * (store_row() >= 0), read/write the velocity (and angular-velocity)
     * columns directly through the raw column pointers rather than through the
     * per-field accessor proxies -- this is the mixed-parts hot path. An
     * unattached particle (resort-step fresh ghost) or an absent context falls
     * back to the accessor. Local copies are made for/from the archive, so the
     * assignment semantics are preserved. */
    auto const row = p.store_row();
    bool const use_ctx = mom_ctx != nullptr and row >= 0;
#ifndef NDEBUG
    // When use_ctx is true the particle must belong to the active store: a
    // stale mom_ctx from a different store generation or a different store
    // object would silently corrupt the wrong memory region.
    if (use_ctx) {
      assert(p.store() == active_particle_store() and
             "serialize_and_reduce: use_ctx row used but p.store() != active "
             "store");
    }
#endif
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d velocity{};
      ar & velocity;
      if (use_ctx) {
        auto *v = mom_ctx->vel +
                  static_cast<std::size_t>(row) * mom_ctx->vel_row_stride;
        v[0 * mom_ctx->vel_comp_stride] = velocity[0];
        v[1 * mom_ctx->vel_comp_stride] = velocity[1];
        v[2 * mom_ctx->vel_comp_stride] = velocity[2];
      } else {
        p.v() = velocity;
      }
    } else {
      Utils::Vector3d velocity{};
      if (use_ctx) {
        auto const *v = mom_ctx->vel +
                        static_cast<std::size_t>(row) * mom_ctx->vel_row_stride;
        velocity = {v[0 * mom_ctx->vel_comp_stride],
                    v[1 * mom_ctx->vel_comp_stride],
                    v[2 * mom_ctx->vel_comp_stride]};
      } else {
        velocity = p.v();
      }
      ar & velocity;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d angular_velocity{};
      ar & angular_velocity;
      if (use_ctx) {
        auto *w = mom_ctx->omega +
                  static_cast<std::size_t>(row) * mom_ctx->omega_row_stride;
        w[0 * mom_ctx->omega_comp_stride] = angular_velocity[0];
        w[1 * mom_ctx->omega_comp_stride] = angular_velocity[1];
        w[2 * mom_ctx->omega_comp_stride] = angular_velocity[2];
      } else {
        p.omega() = angular_velocity;
      }
    } else {
      Utils::Vector3d angular_velocity{};
      if (use_ctx) {
        auto const *w = mom_ctx->omega + static_cast<std::size_t>(row) *
                                             mom_ctx->omega_row_stride;
        angular_velocity = {w[0 * mom_ctx->omega_comp_stride],
                            w[1 * mom_ctx->omega_comp_stride],
                            w[2 * mom_ctx->omega_comp_stride]};
      } else {
        angular_velocity = p.omega();
      }
      ar & angular_velocity;
    }
#endif
  }
  if (data_parts & GHOSTTRANS_FORCE) {
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d force{};
      ar & force;
      if (policy == ReductionPolicy::UPDATE) {
        p.force() += force;
      } else {
        p.force() = force;
      }
    } else {
      Utils::Vector3d force = p.force();
      ar & force;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d torque{};
      ar & torque;
      if (policy == ReductionPolicy::UPDATE) {
        p.torque() += torque;
      } else {
        p.torque() = torque;
      }
    } else {
      Utils::Vector3d torque = p.torque();
      ar & torque;
    }
#endif
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (data_parts & GHOSTTRANS_RATTLE) {
    // rattle_correction() is a ParticleStore column accessor (returns a
    // VectorReference / value proxy, not a raw Vector3d&), so a local Vector3d
    // is used for/from the archive -- never bind the proxy to the archive
    // (mirrors the FORCE branch above). Wire layout is one Utils::Vector3d
    // (24 B).
    if (policy == ReductionPolicy::UPDATE and
        direction == SerializationDirection::LOAD) {
      Utils::Vector3d correction{};
      ar & correction;
      p.rattle_correction() += correction;
    } else if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d correction{};
      ar & correction;
      p.rattle_correction() = correction;
    } else {
      // SAVE. rattle_correction() asserts an attached particle (it is a
      // column-only observable -- like force()). The only
      // caller that reaches here with a DETACHED particle is the
      // SerializationSizeCalculator (a throwaway Particle p{}), which cares
      // only about the byte count; a zero placeholder there is correct because
      // the sizer never inspects the value. Real ghost packs always run on
      // attached particles, so the live column value is written.
      Utils::Vector3d correction = (p.store() != nullptr)
                                       ? Utils::Vector3d(p.rattle_correction())
                                       : Utils::Vector3d{0., 0., 0.};
      ar & correction;
    }
  }
#endif
}

/* ------------------------------------------------------------------------- *
 *  Columnar fast paths for the two per-step ghost communications.
 *
 *  The two hot, per-integration-step ghost transfers carry exactly one field
 *  group each: POSITION (ghost position distribution) and FORCE (ghost force
 *  reduction). For these the ParticleStore state/observable columns can be
 *  streamed to/from the CommBuf in bulk over contiguous store-row ranges,
 *  bypassing the per-Particle proxy construction and per-field archive
 *  dispatch that the generic serialize_and_reduce path incurs.
 *
 *  Wire layout is BYTE-IDENTICAL to serialize_and_reduce (the MemcpyArchive
 *  packs each field as its raw component bytes, tightly, no padding):
 *    POSITION: pos (3 double) | image (3 int)
 *              [| quat (4 double)] [| pos_last_time_step (3 double)]
 *    FORCE:    force (3 double) [| torque (3 double)]
 *  so buffer sizes (calc_transmit_size, untouched) and the GHOST_RDCE bitwise
 *  reduction over raw doubles stay identical.
 *
 *  Applicability (see ghost_communicator): the fast path is taken only when
 *  data_parts is EXACTLY GHOSTTRANS_POSITION or GHOSTTRANS_FORCE, AND the
 *  ParticleStore is not dirty. FORCE reduction always runs on a clean store
 *  (calculate_forces syncs the store, then reduces with no intervening
 *  topology change). POSITION distribution (ghosts_count / ghosts_update) runs
 *  on a clean store for non-resort steps. On resort steps the store is rebuilt
 *  between ghosts_count and ghosts_update: resort-step ghost updates therefore
 *  run with an attached (clean) store and the fast path engages there. The
 *  dirty window (mark_particle_store_dirty before decomposition resort,
 *  re-synced at calculate_forces) is where the mixed-ctx path in
 *  serialize_and_reduce is exercised for freshly attached ghosts. We simply
 *  fall back to the per-particle path when the store is dirty.
 * ------------------------------------------------------------------------- */

/** @brief The ParticleStore backing the current cell structure, or nullptr if
 *  no system/cell structure is set up (e.g. unit-test harnesses). */
static ParticleStore *active_particle_store() {
  auto const &system = System::get_system();
  if (system.cell_structure == nullptr) {
    return nullptr;
  }
  return &system.cell_structure->particle_store();
}

/** @brief Contiguous store-row range [first_row, first_row + size) for a
 *  part_list, or std::nullopt if the list is empty or its first particle is
 *  detached (store()==nullptr / store_row() < 0).
 *
 *  Contiguity is a STRUCTURAL invariant, not something we scan for in release:
 *  ensure_particle_store_synchronized assigns rows in cell-traversal order
 *  (every local cell's particles, then every ghost cell's particles), and every
 *  ghost-comm part_list is exactly one cell's ParticleList, so each cell's
 *  particles occupy one consecutive row block. The fast path is only reachable
 *  on a clean store (columnar_eligible), where this holds by construction. We
 *  therefore derive the range from the FIRST particle alone
 *  (first_row = first->store_row(); size = list size) with no per-particle
 * scan. A detached first particle (empty list or dirty/fresh-ghost store)
 * yields nullopt, detectable BEFORE any caller mutates the store, and callers
 * then fall back to the generic per-particle communication for the whole comm.
 *
 *  In debug builds ONLY, the whole list is validated up front (before any
 * caller reads/writes a column) and a contiguity/attachment violation asserts
 * with a clear message. */
static std::optional<std::pair<int, std::size_t>>
contiguous_store_rows(Cell const &cell) {
  // The cell's committed rows ARE the contiguous store range
  // [offset, offset+count) by construction (the permutation rebuild lays cells
  // back-to-back), so contiguity is structural rather than checked. An
  // empty range yields nullopt. This runs on a clean store (no pending-removed
  // rows), so the raw range is the live range.
  auto const size = cell.count();
  if (size == 0u) {
    return std::nullopt;
  }
  auto const first_row = static_cast<int>(cell.offset());
  return std::make_pair(first_row, size);
}

/** @brief Resolve the contiguous store-row range for every part_list in a comm,
 *  deciding fast-path applicability BEFORE any caller mutates the store.
 *
 *  Returns std::nullopt (fall back to the generic path for the WHOLE comm) if
 *  any non-empty part_list has a detached first particle; empty lists yield an
 *  entry with size 0. Reads only each list's first particle (no per-particle
 *  scan): contiguity within a list is a structural invariant on a clean store
 *  (see contiguous_store_rows). Because the fall-back decision is taken here,
 *  before any column is read or written, the mutating callers (FORCE `+=`,
 *  POSITION assign) can never partially apply and then fall back -- which for
 *  FORCE would double-count. */
static std::vector<std::pair<int, std::size_t>> const *
columnar_resolve_ranges(GhostCommunication const &ghost_comm,
                        ParticleStore const &store) {
  /* Rows only change at a rebuild (store generation bump); the resolved
   * ranges per communication are cached until then, so steady-state steps
   * pay neither the per-list first-particle dereference (a cache-miss
   * pointer chase) nor a heap allocation. */
  if (ghost_comm.cached_ranges_valid and
      ghost_comm.cached_store_generation == store.generation() and
      ghost_comm.cached_store == &store) {
#ifndef NDEBUG
    // Contiguity canary: recompute ranges from scratch and assert equality.
    // This verifies that the cache has not gone stale while the store identity
    // and generation both remained constant (e.g. a mid-step realloc that
    // preserves the generation counter). The debug overhead is acceptable; the
    // release path skips all of this.
    std::vector<std::pair<int, std::size_t>> check_ranges;
    check_ranges.reserve(ghost_comm.part_lists.size());
    for (auto part_list : ghost_comm.part_lists) {
      auto const range = contiguous_store_rows(*part_list);
      if (not range.has_value()) {
        if (part_list->rows().empty()) {
          check_ranges.emplace_back(0, 0u);
          continue;
        }
        assert(false and
               "cached range became detached: cache invariant violated");
        return nullptr;
      }
      check_ranges.push_back(*range);
    }
    assert(check_ranges == ghost_comm.cached_ranges and
           "cached columnar ranges diverged from live contiguous_store_rows");
#endif
    return &ghost_comm.cached_ranges;
  }
  ghost_comm.cached_ranges_valid = false;
  auto &ranges = ghost_comm.cached_ranges;
  ranges.clear();
  ranges.reserve(ghost_comm.part_lists.size());
  for (auto part_list : ghost_comm.part_lists) {
    auto const range = contiguous_store_rows(*part_list);
    if (not range.has_value()) {
      if (part_list->rows().empty()) {
        ranges.emplace_back(0, 0u);
        continue;
      }
      return nullptr; // empty/unsynced cell: fall back before mutating
    }
    ranges.push_back(*range);
  }
  ghost_comm.cached_store_generation = store.generation();
  ghost_comm.cached_store = &store;
  ghost_comm.cached_ranges_valid = true;
  return &ghost_comm.cached_ranges;
}

/** @brief Whether data_parts selects a columnar-eligible per-step comm on a
 *  clean store. When false, callers use the generic per-particle path. */
static bool columnar_eligible(unsigned int data_parts) {
  if (data_parts != GHOSTTRANS_POSITION and data_parts != GHOSTTRANS_FORCE) {
    return false;
  }
  auto const *store = active_particle_store();
  return store != nullptr and not store->is_dirty();
}

/* ------------------------------------------------------------------------- *
 *  Per-communication row-access context for the MIXED-parts POSITION branch.
 *
 *  The per-step ghost update slot NEVER carries a pure GHOSTTRANS_POSITION
 *  transfer: get_global_ghost_flags() always yields POSITION|PROPERTIES, so
 *  the whole-communication POSITION-only fast path above cannot engage for it.
 *  It therefore goes through the generic per-particle serialize_and_reduce,
 *  whose POSITION branch would otherwise route every pos/image/quat/pos_last
 *  field through a per-field VectorReference/QuaternionReference proxy (SAVE:
 *  proxy construct + materialize per field; LOAD: proxy construct + write per
 *  field), ~2k ghosts x both sides x per step.
 *
 *  The context captures, ONCE per ghost_communicator/cell_cell_transfer
 *  invocation (only when the store is present and clean, WITHOUT the
 *  data_parts==POSITION restriction of columnar_eligible), the raw host base
 *  pointers + row/component strides for the four POSITION columns. Any attached
 *  particle's row-r vector then lives at base + r*row_stride, its components at
 *  component_stride apart (component_stride==1 and the row contiguous under the
 *  store's LayoutRight, matching the accessor proxies' invariants). The generic
 *  POSITION branch reads/writes these directly by row instead of through the
 *  accessor, producing a byte-identical wire layout.
 *
 *  Passed to serialize_and_reduce as an optional pointer (nullptr = accessor
 *  fallback), and consulted only when p.store_row() >= 0 (an attached
 *  particle). Detached ghosts during a resort-step update (store()==nullptr /
 *  store_row() < 0) and any dirty-store situation fall back to the accessor,
 *  which falls back to the accessor.
 * ------------------------------------------------------------------------- */

/** @brief Build a PositionRowContext for the current clean store, or
 *  std::nullopt if there is no store or it is dirty (accessor fallback). This
 *  intentionally omits the data_parts==POSITION gate of columnar_eligible: it
 *  is meant for the mixed-parts (POSITION|PROPERTIES) per-step update. */
static std::optional<PositionRowContext> make_position_row_context() {
  auto *store = active_particle_store();
  if (store == nullptr or store->is_dirty() or
      store->number_of_particles() == 0u) {
    return std::nullopt;
  }
  auto pos = store->position_view();
  auto img = store->image_box_view();
  PositionRowContext ctx{};
  ctx.pos = &pos(0, 0);
  ctx.img = &img(0, 0);
  ctx.pos_row_stride = pos.stride(0);
  ctx.pos_comp_stride = pos.stride(1);
  ctx.img_row_stride = img.stride(0);
  ctx.img_comp_stride = img.stride(1);
#ifdef ESPRESSO_ROTATION
  auto quat = store->quaternion_view();
  ctx.quat = &quat(0, 0);
  ctx.quat_row_stride = quat.stride(0);
  ctx.quat_comp_stride = quat.stride(1);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto plast = store->position_last_time_step_view();
  ctx.plast = &plast(0, 0);
  ctx.plast_row_stride = plast.stride(0);
  ctx.plast_comp_stride = plast.stride(1);
#endif
  return ctx;
}

/** @brief Build a ForceRowContext for the current store (caller has already
 *  established store validity/cleanliness). Returns a zero-initialized context
 *  when the store is empty; callers' loops skip n==0 ranges so the null
 *  pointers are never dereferenced, but taking &view(0,0) on an empty view is
 *  UB under bounds-checking. */
static ForceRowContext make_force_row_context(ParticleStore &store) {
  if (store.number_of_particles() == 0u) {
    return ForceRowContext{};
  }
  ForceRowContext ctx{};
  auto force = store.force_view();
  ctx.force = &force(0, 0);
  ctx.force_row_stride = force.stride(0);
  ctx.force_comp_stride = force.stride(1);
#ifdef ESPRESSO_ROTATION
  auto torque = store.torque_view();
  ctx.torque = &torque(0, 0);
  ctx.torque_row_stride = torque.stride(0);
  ctx.torque_comp_stride = torque.stride(1);
#endif
  return ctx;
}

/** @brief Build a PositionRowContext unconditionally from a valid store
 *  (helper for the columnar bulk paths, which have already checked
 *  cleanliness via columnar_eligible/columnar_resolve_ranges). Returns a
 *  zero-initialized context when the store is empty; callers' loops skip n==0
 *  ranges so the null pointers are never dereferenced, but taking &view(0,0)
 *  on an empty view is UB under bounds-checking. */
static PositionRowContext
make_position_row_context_unchecked(ParticleStore &store) {
  if (store.number_of_particles() == 0u) {
    return PositionRowContext{};
  }
  PositionRowContext ctx{};
  auto pos = store.position_view();
  auto img = store.image_box_view();
  ctx.pos = &pos(0, 0);
  ctx.img = &img(0, 0);
  ctx.pos_row_stride = pos.stride(0);
  ctx.pos_comp_stride = pos.stride(1);
  ctx.img_row_stride = img.stride(0);
  ctx.img_comp_stride = img.stride(1);
#ifdef ESPRESSO_ROTATION
  auto quat = store.quaternion_view();
  ctx.quat = &quat(0, 0);
  ctx.quat_row_stride = quat.stride(0);
  ctx.quat_comp_stride = quat.stride(1);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto plast = store.position_last_time_step_view();
  ctx.plast = &plast(0, 0);
  ctx.plast_row_stride = plast.stride(0);
  ctx.plast_comp_stride = plast.stride(1);
#endif
  return ctx;
}

/** @brief Build a MomentumRowContext for the current clean store, or
 *  std::nullopt if there is no store or it is dirty (accessor fallback).
 *  Mirror of make_position_row_context for the velocity/angular-velocity
 *  columns; meant for the mixed-parts (POSITION|PROPERTIES|MOMENTUM) per-step
 *  update.
 *
 */
static std::optional<MomentumRowContext> make_momentum_row_context() {
  // The ParticleStore velocity/angular-velocity columns are authoritative
  // (Particle::v()/omega() read/write them), so this context reads/writes the
  // same memory as the accessors. The MOMENTUM ghost value path routes through
  // the columnar row-context when the store is clean and attached, exactly
  // like the POSITION path.
  static constexpr bool velocity_columns_authoritative = true;
  auto *store = active_particle_store();
  if (not velocity_columns_authoritative or store == nullptr or
      store->is_dirty() or store->number_of_particles() == 0u) {
    return std::nullopt;
  }
  MomentumRowContext ctx{};
  auto vel = store->velocity_view();
  ctx.vel = &vel(0, 0);
  ctx.vel_row_stride = vel.stride(0);
  ctx.vel_comp_stride = vel.stride(1);
#ifdef ESPRESSO_ROTATION
  auto omega = store->angular_velocity_view();
  ctx.omega = &omega(0, 0);
  ctx.omega_row_stride = omega.stride(0);
  ctx.omega_comp_stride = omega.stride(1);
#endif
  return ctx;
}

/** @brief Build a ParameterRowContext for the current clean store, or
 *  std::nullopt if there is no store or it is dirty (accessor fallback).
 *  Mirror of make_momentum_row_context for the PROPERTIES group: the hot
 *  parameter columns (raw base pointers + strides) and a ParticleStore* for the
 *  cold POD sidecars. Meant for the mixed-parts (POSITION|PROPERTIES|MOMENTUM)
 *  per-step update and the resort-driven / reaction-driven PROPERTIES comms.
 *
 */
static std::optional<ParameterRowContext> make_parameter_row_context() {
  // The ParticleStore parameter columns and host sidecars are authoritative
  // (the Particle parameter accessors read/write them), so this context
  // reads/writes the same memory as the accessors. The PROPERTIES ghost value
  // path routes through the columnar row-context when the store is clean and
  // attached, exactly like the POSITION/MOMENTUM paths.
  auto *store = active_particle_store();
  if (store == nullptr or store->is_dirty() or
      store->number_of_particles() == 0u) {
    return std::nullopt;
  }
  ParameterRowContext ctx{};
  ctx.store = store;
  ctx.id = &store->id_view()(0);
  ctx.mol_id = &store->mol_id_view()(0);
  ctx.type = &store->type_view()(0);
  ctx.propagation = &store->propagation_view()(0);
#ifdef ESPRESSO_ROTATION
  ctx.rotation = &store->rotation_view()(0);
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  auto rinertia = store->rinertia_view();
  ctx.rinertia = &rinertia(0, 0);
  ctx.rinertia_row_stride = rinertia.stride(0);
  ctx.rinertia_comp_stride = rinertia.stride(1);
#endif
#endif
#ifdef ESPRESSO_MASS
  ctx.mass = &store->mass_view()(0);
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  ctx.q = &store->q_view()(0);
#endif
#ifdef ESPRESSO_DIPOLES
  ctx.dipm = &store->dipm_view()(0);
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  auto mu_E = store->mu_E_view();
  ctx.mu_E = &mu_E(0, 0);
  ctx.mu_E_row_stride = mu_E.stride(0);
  ctx.mu_E_comp_stride = mu_E.stride(1);
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto gamma = store->gamma_view();
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  ctx.gamma = &gamma(0, 0);
  ctx.gamma_row_stride = gamma.stride(0);
  ctx.gamma_comp_stride = gamma.stride(1);
#else
  ctx.gamma = &gamma(0);
#endif
#ifdef ESPRESSO_ROTATION
  auto gamma_rot = store->gamma_rot_view();
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  ctx.gamma_rot = &gamma_rot(0, 0);
  ctx.gamma_rot_row_stride = gamma_rot.stride(0);
  ctx.gamma_rot_comp_stride = gamma_rot.stride(1);
#else
  ctx.gamma_rot = &gamma_rot(0);
#endif
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  ctx.ext_flag = &store->ext_flag_view()(0);
  auto ext_force = store->ext_force_view();
  ctx.ext_force = &ext_force(0, 0);
  ctx.ext_force_row_stride = ext_force.stride(0);
  ctx.ext_force_comp_stride = ext_force.stride(1);
#ifdef ESPRESSO_ROTATION
  auto ext_torque = store->ext_torque_view();
  ctx.ext_torque = &ext_torque(0, 0);
  ctx.ext_torque_row_stride = ext_torque.stride(0);
  ctx.ext_torque_comp_stride = ext_torque.stride(1);
#endif
#endif
  return ctx;
}

/* No make_momentum_row_context_unchecked / pack_momentum_range /
 * unpack_momentum_range helpers: the only PURE GHOSTTRANS_MOMENTUM comm is
 * RATTLE's correct_velocity_shake ghosts_update(DATA_PART_MOMENTUM), a
 * low-frequency correction loop, not a per-step hot path like the POSITION
 * distribution or the FORCE reduction that justify their columnar bulk paths.
 * The mixed-parts MomentumRowContext
 * in serialize_and_reduce covers every MOMENTUM comm (pure and mixed) without a
 * separate whole-communication columnar path; if profiling later flags the
 * velocity shake, the _unchecked builder + pack/unpack_momentum_range + a
 * columnar_eligible(GHOSTTRANS_MOMENTUM) branch can be added as for FORCE. */

namespace {
/** @brief Byte cursor over the CommBuf that packs/unpacks scalars via memcpy,
 *  exactly like the MemcpyArchive. Using memcpy (not a reinterpret_cast +
 *  dereference) is REQUIRED: the tightly packed per-particle wire layout is not
 *  8-byte aligned (e.g. 36 B/particle for POSITION without ROTATION), so a
 *  direct double* dereference would be misaligned UB (and trip UBSan). */
template <class CharPtr> struct ByteCursor {
  CharPtr ptr;
  template <class T> void put(T const &value) {
    std::memcpy(ptr, &value, sizeof(T));
    ptr += sizeof(T);
  }
  template <class T> T get() {
    T value;
    std::memcpy(&value, ptr, sizeof(T));
    ptr += sizeof(T);
    return value;
  }
};
using WriteCursor = ByteCursor<char *>;
using ReadCursor = ByteCursor<char const *>;
} // namespace

/** @brief Bulk-pack POSITION for a contiguous row range into @p out (SAVE).
 *  Reproduces serialize_and_reduce's SAVE branch exactly: with a ghost shift,
 *  apply +shift then fold_position (writing the folded image); otherwise copy
 *  the stored pos/image. Then quat (ROTATION) and pos_last_time_step
 *  (BOND_CONSTRAINT) verbatim. */
static char *pack_position_range(char *out, PositionRowContext const &ctx,
                                 int first_row, std::size_t n,
                                 BoxGeometry const &box_geo,
                                 Utils::Vector3d const *ghost_shift) {
  WriteCursor cursor{out};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = static_cast<std::size_t>(first_row) + i;
    auto const *p = ctx.pos + row * ctx.pos_row_stride;
    auto const *im = ctx.img + row * ctx.img_row_stride;
    Utils::Vector3d position{p[0], p[ctx.pos_comp_stride],
                             p[2u * ctx.pos_comp_stride]};
    Utils::Vector3i image_box{im[0], im[ctx.img_comp_stride],
                              im[2u * ctx.img_comp_stride]};
    if (ghost_shift != nullptr) {
      position += *ghost_shift;
      box_geo.fold_position(position, image_box);
    }
    cursor.put(position[0]);
    cursor.put(position[1]);
    cursor.put(position[2]);
    cursor.put(image_box[0]);
    cursor.put(image_box[1]);
    cursor.put(image_box[2]);
#ifdef ESPRESSO_ROTATION
    auto const *q = ctx.quat + row * ctx.quat_row_stride;
    cursor.put(q[0]);
    cursor.put(q[ctx.quat_comp_stride]);
    cursor.put(q[2u * ctx.quat_comp_stride]);
    cursor.put(q[3u * ctx.quat_comp_stride]);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    auto const *pl = ctx.plast + row * ctx.plast_row_stride;
    cursor.put(pl[0]);
    cursor.put(pl[ctx.plast_comp_stride]);
    cursor.put(pl[2u * ctx.plast_comp_stride]);
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-unpack POSITION for a contiguous row range from @p in (LOAD).
 *  Assigns pos/image (and quat, pos_last_time_step) directly into columns. */
static char const *unpack_position_range(char const *in,
                                         PositionRowContext const &ctx,
                                         int first_row, std::size_t n) {
  ReadCursor cursor{in};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = static_cast<std::size_t>(first_row) + i;
    auto *p = ctx.pos + row * ctx.pos_row_stride;
    auto *im = ctx.img + row * ctx.img_row_stride;
    p[0] = cursor.get<double>();
    p[ctx.pos_comp_stride] = cursor.get<double>();
    p[2u * ctx.pos_comp_stride] = cursor.get<double>();
    im[0] = cursor.get<int>();
    im[ctx.img_comp_stride] = cursor.get<int>();
    im[2u * ctx.img_comp_stride] = cursor.get<int>();
#ifdef ESPRESSO_ROTATION
    auto *q = ctx.quat + row * ctx.quat_row_stride;
    q[0] = cursor.get<double>();
    q[ctx.quat_comp_stride] = cursor.get<double>();
    q[2u * ctx.quat_comp_stride] = cursor.get<double>();
    q[3u * ctx.quat_comp_stride] = cursor.get<double>();
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    auto *pl = ctx.plast + row * ctx.plast_row_stride;
    pl[0] = cursor.get<double>();
    pl[ctx.plast_comp_stride] = cursor.get<double>();
    pl[2u * ctx.plast_comp_stride] = cursor.get<double>();
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-pack FORCE for a contiguous row range into @p out (SAVE). */
static char *pack_force_range(char *out, ForceRowContext const &ctx,
                              int first_row, std::size_t n) {
  WriteCursor cursor{out};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = static_cast<std::size_t>(first_row) + i;
    auto const *f = ctx.force + row * ctx.force_row_stride;
    cursor.put(f[0]);
    cursor.put(f[ctx.force_comp_stride]);
    cursor.put(f[2u * ctx.force_comp_stride]);
#ifdef ESPRESSO_ROTATION
    auto const *t = ctx.torque + row * ctx.torque_row_stride;
    cursor.put(t[0]);
    cursor.put(t[ctx.torque_comp_stride]);
    cursor.put(t[2u * ctx.torque_comp_stride]);
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-unpack FORCE for a contiguous row range from @p in.
 *  With @p accumulate (UPDATE policy) the received value is ADDED; otherwise
 *  (MOVE policy) it is assigned. Matches serialize_and_reduce's FORCE LOAD. */
static char const *unpack_force_range(char const *in,
                                      ForceRowContext const &ctx, int first_row,
                                      std::size_t n, bool accumulate) {
  ReadCursor cursor{in};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = static_cast<std::size_t>(first_row) + i;
    auto *f = ctx.force + row * ctx.force_row_stride;
    if (accumulate) {
      f[0] += cursor.get<double>();
      f[ctx.force_comp_stride] += cursor.get<double>();
      f[2u * ctx.force_comp_stride] += cursor.get<double>();
    } else {
      f[0] = cursor.get<double>();
      f[ctx.force_comp_stride] = cursor.get<double>();
      f[2u * ctx.force_comp_stride] = cursor.get<double>();
    }
#ifdef ESPRESSO_ROTATION
    auto *t = ctx.torque + row * ctx.torque_row_stride;
    if (accumulate) {
      t[0] += cursor.get<double>();
      t[ctx.torque_comp_stride] += cursor.get<double>();
      t[2u * ctx.torque_comp_stride] += cursor.get<double>();
    } else {
      t[0] = cursor.get<double>();
      t[ctx.torque_comp_stride] = cursor.get<double>();
      t[2u * ctx.torque_comp_stride] = cursor.get<double>();
    }
#endif
  }
  return cursor.ptr;
}

static auto calc_transmit_size(BoxGeometry const &box_geo,
                               unsigned data_parts) {
  std::size_t properties_size = 0ul;
  if (data_parts & GHOSTTRANS_PROPRTS) {
    /* Compositional, compile-time-constant size of the PROPERTIES wire layout,
     * matching the field order/ifdef structure of serialize_and_reduce's
     * PROPERTIES branch exactly. Each `ar & <field>` packs sizeof(field)
     * tightly (the MemcpyArchive uses memcpy for trivially/bitwise-serializable
     * types, no alignment padding between fields), and
     * SerializationSizeCalculator likewise accumulates sizeof(T), so this
     * constant matches the sizer output byte-for-byte. The three cold PODs are
     * bitwise-serializable, so each is packed whole as sizeof(T) INCLUDING
     * internal struct padding (e.g. VirtualSitesRelativeParameters is 80 B: int
     * + 4 B pad + double + 2 quaternions, NOT 76). dip_fld and magnetodynamics
     * are NOT part of the PROPERTIES ghost group (see the branch), so they are
     * excluded here too. */
    properties_size = 4ul * sizeof(int); // id, mol_id, type, propagation
#ifdef ESPRESSO_ROTATION
    properties_size += sizeof(std::uint8_t); // rotation
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    properties_size += sizeof(Utils::Vector3d); // rinertia
#endif
#endif
#ifdef ESPRESSO_MASS
    properties_size += sizeof(double); // mass
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    properties_size += sizeof(double); // q
#endif
#ifdef ESPRESSO_DIPOLES
    properties_size += sizeof(double); // dipm
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    properties_size += sizeof(Utils::Vector3d); // mu_E
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    properties_size += sizeof(VirtualSitesRelativeParameters); // vs_relative
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    properties_size += sizeof(ParticleStore::GammaValue); // gamma
#ifdef ESPRESSO_ROTATION
    properties_size += sizeof(ParticleStore::GammaValue); // gamma_rot
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    properties_size += sizeof(std::uint8_t);    // ext_flag / fixed
    properties_size += sizeof(Utils::Vector3d); // ext_force
#ifdef ESPRESSO_ROTATION
    properties_size += sizeof(Utils::Vector3d); // ext_torque
#endif
#endif
#ifdef ESPRESSO_ENGINE
    properties_size += sizeof(ParticleParametersSwimming); // swimming
#endif
    data_parts &= ~static_cast<unsigned>(GHOSTTRANS_PROPRTS);
    if (data_parts == 0u) {
      return properties_size;
    }
  }
  std::size_t force_size = 0ul;
  if (data_parts & GHOSTTRANS_FORCE) {
#ifdef ESPRESSO_ROTATION
    force_size = 6ul * sizeof(double);
#else
    force_size = 3ul * sizeof(double);
#endif
    data_parts &= ~static_cast<unsigned>(GHOSTTRANS_FORCE);
    if (data_parts == 0u) {
      return properties_size + force_size;
    }
  }
  std::size_t position_size = 0ul;
  if (data_parts & GHOSTTRANS_POSITION) {
    /* Compositional, compile-time-constant size of the POSITION wire layout:
     * pos (Vector3d, 24 B) + image_box (Vector3i, 12 B)
     * [+ quat (Quaternion<double>, 32 B) under ROTATION]
     * [+ pos_last_time_step (Vector3d, 24 B) under BOND_CONSTRAINT].
     * The MemcpyArchive packs every field tightly (exactly sizeof(T), no
     * alignment padding), and SerializationSizeCalculator likewise accumulates
     * sizeof(T), so this constant matches the sizer output exactly. */
    position_size = 3ul * sizeof(double) + 3ul * sizeof(int);
#ifdef ESPRESSO_ROTATION
    position_size += sizeof(Utils::Quaternion<double>);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    position_size += 3ul * sizeof(double);
#endif
    data_parts &= ~static_cast<unsigned>(GHOSTTRANS_POSITION);
    if (data_parts == 0u) {
      return properties_size + force_size + position_size;
    }
  }
  std::size_t momentum_size = 0ul;
  if (data_parts & GHOSTTRANS_MOMENTUM) {
    /* Compositional, compile-time-constant size of the MOMENTUM wire layout:
     * v (Vector3d, 24 B) [+ omega (Vector3d, 24 B) under ROTATION].
     * serialize_and_reduce's MOMENTUM branch archives one Vector3d (velocity)
     * and, under ROTATION, a second (angular velocity); each is 3 tightly
     * packed doubles, so this constant matches the sizer's accumulation of
     * sizeof(T) exactly. */
    momentum_size = 3ul * sizeof(double);
#ifdef ESPRESSO_ROTATION
    momentum_size += 3ul * sizeof(double);
#endif
    data_parts &= ~static_cast<unsigned>(GHOSTTRANS_MOMENTUM);
    if (data_parts == 0u) {
      return properties_size + force_size + position_size + momentum_size;
    }
  }
  /* Only GHOSTTRANS_RATTLE (under ESPRESSO_BOND_CONSTRAINT) can still reach
   * here: PROPRTS/FORCE/POSITION/MOMENTUM are masked above, PARTNUM is handled
   * in the ghost_comm overload, and BONDS travels in the separate bond buffer.
   * The RATTLE branch archives exactly one Utils::Vector3d (rattle_correction)
   * on MOVE/SAVE, which the MemcpyArchive packs tightly as 3 doubles = 24 B --
   * the same compositional constant the POSITION/MOMENTUM legs above use, so
   * the generic SerializationSizeCalculator fallback (and a live Particle) is
   * not needed here. */
  std::size_t rattle_size = 0ul;
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (data_parts & GHOSTTRANS_RATTLE) {
    rattle_size = 3ul * sizeof(double);
  }
#endif
  /* Accounting check: every part must be handled above or contribute zero
   * bytes to the main transfer buffer. GHOSTTRANS_BONDS is tolerated because
   * the bond payload travels in the dedicated ragged bond buffer
   * (CommBuf::bondbuf), sized at pack time; a communicator's data_parts
   * legitimately carries the BONDS bit (e.g. the exchange communicator when
   * bonded interactions and a thermostat are active) -- caught by the
   * Debug-with-asserts checkpoint tests on 4 MPI ranks (upstream fedora CI).
   * RATTLE is tolerated only when its ifdef branch above accounts for it. */
  [[maybe_unused]] constexpr auto tolerated_parts =
      static_cast<unsigned>(GHOSTTRANS_BONDS)
#ifdef ESPRESSO_BOND_CONSTRAINT
      | static_cast<unsigned>(GHOSTTRANS_RATTLE)
#endif
      ;
  assert((data_parts & ~tolerated_parts) == 0u);
  return properties_size + force_size + position_size + momentum_size +
         rattle_size;
}

static auto calc_transmit_size(GhostCommunication const &ghost_comm,
                               BoxGeometry const &box_geo,
                               unsigned int data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM)
    return sizeof(unsigned int) * ghost_comm.part_lists.size();

  auto const n_part = boost::accumulate(
      ghost_comm.part_lists, std::size_t{0},
      [](std::size_t sum, auto part_list) { return sum + part_list->count(); });
  return n_part * calc_transmit_size(box_geo, data_parts);
}

/** @brief Columnar bulk fill of @p send_buffer for POSITION/FORCE on a clean
 *  store. Returns false (leaving the buffer untouched) if any part_list breaks
 *  the contiguity/attachment precondition, so the caller can fall back. */
static bool columnar_prepare_send_buffer(CommBuf &send_buffer,
                                         GhostCommunication const &ghost_comm,
                                         BoxGeometry const &box_geo,
                                         unsigned int data_parts) {
  auto &store = *active_particle_store();
  auto const *ranges = columnar_resolve_ranges(ghost_comm, store);
  if (ranges == nullptr) {
    return false; // decided before touching the buffer/store
  }
  // data_parts is a single fixed value for the whole columnar-eligible call
  // (POSITION xor FORCE), so build the one plain-value row context the branch
  // needs rather than a per-branch std::optional. Using the value type
  // directly (like the make_*_row_context_unchecked callers in
  // columnar_put_recv_buffer / columnar_add_forces_from_recv_buffer) keeps
  // gcc-14 from mis-tracking the std::optional payload as maybe-uninitialized
  // (the optional-payload false-positive family).
  auto *cursor = send_buffer.data();
  if (data_parts == GHOSTTRANS_POSITION) {
    auto const position_ctx = make_position_row_context_unchecked(store);
    for (auto const &[first_row, n] : *ranges) {
      if (n == 0u) {
        continue;
      }
      cursor = pack_position_range(cursor, position_ctx, first_row, n, box_geo,
                                   &ghost_comm.shift);
    }
  } else {
    auto const force_ctx = make_force_row_context(store);
    for (auto const &[first_row, n] : *ranges) {
      if (n == 0u) {
        continue;
      }
      cursor = pack_force_range(cursor, force_ctx, first_row, n);
    }
  }
  assert(static_cast<std::size_t>(cursor - send_buffer.data()) ==
         send_buffer.size());
  return true;
}

static void prepare_send_buffer(CommBuf &send_buffer,
                                GhostCommunication const &ghost_comm,
                                BoxGeometry const &box_geo,
                                unsigned int data_parts) {

  /* reallocate send buffer */
  send_buffer.resize(calc_transmit_size(ghost_comm, box_geo, data_parts));
  send_buffer.bonds().clear();

  /* Columnar fast path for the two per-step comms on a clean store. */
  if (columnar_eligible(data_parts) and
      columnar_prepare_send_buffer(send_buffer, ghost_comm, box_geo,
                                   data_parts)) {
    return;
  }

  auto archiver = Utils::MemcpyOArchive{send_buffer.make_span()};

  /* Mixed-parts POSITION fast path: read pos/image/quat/pos_last straight from
   * the (clean) store columns by row instead of via per-field accessor proxies.
   * nullptr when POSITION is absent or the store is dirty/absent (accessor
   * fallback, per-particle). */
  auto const pos_ctx = (data_parts & GHOSTTRANS_POSITION)
                           ? make_position_row_context()
                           : std::nullopt;
  auto const *pos_ctx_ptr = pos_ctx.has_value() ? &(*pos_ctx) : nullptr;
  /* Mixed-parts MOMENTUM fast path: read velocity/omega straight from the
   * (clean) store columns by row. nullptr when MOMENTUM is absent, the store is
   * dirty/absent -- accessor fallback, per-particle. */
  auto const mom_ctx = (data_parts & GHOSTTRANS_MOMENTUM)
                           ? make_momentum_row_context()
                           : std::nullopt;
  auto const *mom_ctx_ptr = mom_ctx.has_value() ? &(*mom_ctx) : nullptr;
  /* Mixed-parts PROPERTIES fast path: read the parameter columns/POD sidecars
   * straight from the (clean) store by row. nullptr when PROPERTIES is absent,
   * the store is dirty/absent -- accessor fallback, per-particle. */
  auto const param_ctx = (data_parts & GHOSTTRANS_PROPRTS)
                             ? make_parameter_row_context()
                             : std::nullopt;
  auto const *param_ctx_ptr = param_ctx.has_value() ? &(*param_ctx) : nullptr;

  /* Construct archive that pushes back to the bond buffer */
  namespace io = boost::iostreams;
  io::stream<io::back_insert_device<std::vector<char>>> os{
      io::back_inserter(send_buffer.bonds())};
  boost::archive::binary_oarchive bond_archiver{os};

  /* put in data */
  for (auto part_list : ghost_comm.part_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      // Total count (committed rows + staged): a source cell resized earlier in
      // this pass as a ghost destination has its ghosts staged, not yet
      // committed. The receiver resizes to this pending count.
      assert(part_list->size() <= std::numeric_limits<unsigned int>::max());
      auto np = static_cast<unsigned int>(part_list->size());
      archiver << np;
    } else {
      for (auto &p : part_list->particles()) {
        serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                             SerializationDirection::SAVE, box_geo,
                             &ghost_comm.shift, pos_ctx_ptr, mom_ctx_ptr,
                             param_ctx_ptr);
        if (data_parts & GHOSTTRANS_BONDS) {
          bond_archiver << p.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_written() == send_buffer.size());
}

static void prepare_recv_buffer(CommBuf &recv_buffer,
                                GhostCommunication const &ghost_comm,
                                BoxGeometry const &box_geo,
                                unsigned int data_parts) {
  /* reallocate recv buffer */
  recv_buffer.resize(calc_transmit_size(ghost_comm, box_geo, data_parts));
  /* clear bond buffer */
  recv_buffer.bonds().clear();
}

/** @brief Columnar bulk unpack of @p recv_buffer into columns for POSITION on
 *  a clean store. Returns false (untouched) on a broken precondition. */
static bool columnar_put_recv_buffer(CommBuf &recv_buffer,
                                     GhostCommunication const &ghost_comm,
                                     unsigned int data_parts) {
  assert(data_parts == GHOSTTRANS_POSITION);
  static_cast<void>(data_parts);
  auto &store = *active_particle_store();
  auto const *ranges = columnar_resolve_ranges(ghost_comm, store);
  if (ranges == nullptr) {
    return false; // decided before writing any column
  }
  auto const position_ctx = make_position_row_context_unchecked(store);
  auto const *cursor = recv_buffer.data();
  for (auto const &[first_row, n] : *ranges) {
    if (n == 0u) {
      continue;
    }
    cursor = unpack_position_range(cursor, position_ctx, first_row, n);
  }
  assert(static_cast<std::size_t>(cursor - recv_buffer.data()) ==
         recv_buffer.size());
  return true;
}

static void put_recv_buffer(CommBuf &recv_buffer,
                            GhostCommunication const &ghost_comm,
                            BoxGeometry const &box_geo,
                            unsigned int data_parts) {
  /* Columnar fast path: pure POSITION on a clean store. */
  if (data_parts == GHOSTTRANS_POSITION and columnar_eligible(data_parts) and
      columnar_put_recv_buffer(recv_buffer, ghost_comm, data_parts)) {
    recv_buffer.bonds().clear();
    return;
  }

  /* put back data */
  auto archiver = Utils::MemcpyIArchive{recv_buffer.make_span()};

  if (data_parts & GHOSTTRANS_PARTNUM) {
    for (auto part_list : ghost_comm.part_lists) {
      unsigned int np;
      archiver >> np;
      CellParticleStorage::resize_ghost_storage(*part_list, np);
    }
  } else {
    /* Mixed-parts POSITION fast path (see prepare_send_buffer): write
     * pos/image/quat/pos_last straight into the (clean) store columns by row.
     * nullptr => accessor fallback (POSITION absent, or dirty/absent store).
     * On resort steps the store is REBUILT before ghosts_update, so this path
     * engages there with the attached, clean store. */
    auto const pos_ctx = (data_parts & GHOSTTRANS_POSITION)
                             ? make_position_row_context()
                             : std::nullopt;
    auto const *pos_ctx_ptr = pos_ctx.has_value() ? &(*pos_ctx) : nullptr;
    /* Mixed-parts MOMENTUM fast path (see prepare_send_buffer): write
     * velocity/omega straight into the (clean) store columns by row.
     * nullptr => accessor fallback. */
    auto const mom_ctx = (data_parts & GHOSTTRANS_MOMENTUM)
                             ? make_momentum_row_context()
                             : std::nullopt;
    auto const *mom_ctx_ptr = mom_ctx.has_value() ? &(*mom_ctx) : nullptr;
    /* Mixed-parts PROPERTIES fast path (see prepare_send_buffer): write the
     * parameter columns/POD sidecars straight into the (clean) store by row.
     * nullptr => accessor fallback. */
    auto const param_ctx = (data_parts & GHOSTTRANS_PROPRTS)
                               ? make_parameter_row_context()
                               : std::nullopt;
    auto const *param_ctx_ptr = param_ctx.has_value() ? &(*param_ctx) : nullptr;
    for (auto part_list : ghost_comm.part_lists) {
      for (auto &p : part_list->particles()) {
        serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                             SerializationDirection::LOAD, box_geo, nullptr,
                             pos_ctx_ptr, mom_ctx_ptr, param_ctx_ptr);
      }
    }
    if (data_parts & GHOSTTRANS_BONDS) {
      namespace io = boost::iostreams;
      io::stream<io::array_source> bond_stream(io::array_source{
          recv_buffer.bonds().data(), recv_buffer.bonds().size()});
      boost::archive::binary_iarchive bond_archiver(bond_stream);

      for (auto part_list : ghost_comm.part_lists) {
        for (auto &p : part_list->particles()) {
          bond_archiver >> p.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_read() == recv_buffer.size());

  recv_buffer.bonds().clear();
}

#ifdef ESPRESSO_BOND_CONSTRAINT
static void
add_rattle_correction_from_recv_buffer(CommBuf &recv_buffer,
                                       const GhostCommunication &ghost_comm) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{recv_buffer.make_span()};
  for (auto &part_list : ghost_comm.part_lists) {
    for (Particle &part : part_list->particles()) {
      // The RATTLE correction is a ParticleStore column; the wire carries
      // one Utils::Vector3d. Reduce it into the local column through the
      // accessor.
      Utils::Vector3d correction{};
      archiver >> correction;
      part.rattle_correction() += correction;
    }
  }
}
#endif

/** @brief Columnar bulk force accumulation from @p recv_buffer on a clean
 *  store. Returns false (untouched) on a broken precondition. */
static bool
columnar_add_forces_from_recv_buffer(CommBuf &recv_buffer,
                                     GhostCommunication const &ghost_comm) {
  /* The fall-back decision is made BEFORE any `+=`: a mid-list fall-back after
   * a partial accumulate would double-count on the generic retry. */
  auto &store = *active_particle_store();
  auto const *ranges = columnar_resolve_ranges(ghost_comm, store);
  if (ranges == nullptr) {
    return false;
  }
  auto const force_ctx = make_force_row_context(store);
  auto const *cursor = recv_buffer.data();
  for (auto const &[first_row, n] : *ranges) {
    if (n == 0u) {
      continue;
    }
    cursor = unpack_force_range(cursor, force_ctx, first_row, n,
                                /*accumulate=*/true);
  }
  assert(static_cast<std::size_t>(cursor - recv_buffer.data()) ==
         recv_buffer.size());
  return true;
}

static void add_forces_from_recv_buffer(CommBuf &recv_buffer,
                                        const GhostCommunication &ghost_comm) {
  /* Columnar fast path: FORCE reduction always runs on a clean store. */
  if (columnar_eligible(GHOSTTRANS_FORCE) and
      columnar_add_forces_from_recv_buffer(recv_buffer, ghost_comm)) {
    return;
  }

  /* put back data */
  auto archiver = Utils::MemcpyIArchive{recv_buffer.make_span()};
  for (auto &part_list : ghost_comm.part_lists) {
    for (Particle &part : part_list->particles()) {
      Utils::Vector3d force{};
      archiver >> force;
      part.force() += force;
#ifdef ESPRESSO_ROTATION
      Utils::Vector3d torque{};
      archiver >> torque;
      part.torque() += torque;
#endif
    }
  }
}

/** @brief Columnar POSITION cell-to-cell transfer: shift+fold the source rows
 *  and assign into the destination rows (quat, pos_last_time_step verbatim).
 *  Matches the generic LOCL SAVE(shift)+LOAD(assign) semantics: the shift is
 *  ALWAYS applied (the generic LOCL SAVE always passes a non-null shift), and
 * it is the zero vector for a local cell-to-cell transfer (no periodic wrap),
 * so the `+= shift` is a no-op before the fold in that common case rather than
 * a branch. */
static void locl_transfer_position(PositionRowContext const &ctx, int src_first,
                                   int dst_first, std::size_t n,
                                   BoxGeometry const &box_geo,
                                   Utils::Vector3d const &shift) {
  for (std::size_t i = 0u; i < n; ++i) {
    auto const s = static_cast<std::size_t>(src_first) + i;
    auto const d = static_cast<std::size_t>(dst_first) + i;
    auto const *ps = ctx.pos + s * ctx.pos_row_stride;
    auto const *is = ctx.img + s * ctx.img_row_stride;
    auto *pd = ctx.pos + d * ctx.pos_row_stride;
    auto *id = ctx.img + d * ctx.img_row_stride;
    Utils::Vector3d position{ps[0], ps[ctx.pos_comp_stride],
                             ps[2u * ctx.pos_comp_stride]};
    Utils::Vector3i image_box{is[0], is[ctx.img_comp_stride],
                              is[2u * ctx.img_comp_stride]};
    position += shift;
    box_geo.fold_position(position, image_box);
    pd[0] = position[0];
    pd[ctx.pos_comp_stride] = position[1];
    pd[2u * ctx.pos_comp_stride] = position[2];
    id[0] = image_box[0];
    id[ctx.img_comp_stride] = image_box[1];
    id[2u * ctx.img_comp_stride] = image_box[2];
#ifdef ESPRESSO_ROTATION
    auto const *qs = ctx.quat + s * ctx.quat_row_stride;
    auto *qd = ctx.quat + d * ctx.quat_row_stride;
    qd[0] = qs[0];
    qd[ctx.quat_comp_stride] = qs[ctx.quat_comp_stride];
    qd[2u * ctx.quat_comp_stride] = qs[2u * ctx.quat_comp_stride];
    qd[3u * ctx.quat_comp_stride] = qs[3u * ctx.quat_comp_stride];
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    auto const *ls = ctx.plast + s * ctx.plast_row_stride;
    auto *ld = ctx.plast + d * ctx.plast_row_stride;
    ld[0] = ls[0];
    ld[ctx.plast_comp_stride] = ls[ctx.plast_comp_stride];
    ld[2u * ctx.plast_comp_stride] = ls[2u * ctx.plast_comp_stride];
#endif
  }
}

/** @brief Columnar FORCE cell-to-cell transfer: accumulate source rows into
 *  destination rows (UPDATE policy: += force and torque). */
static void locl_transfer_force(ForceRowContext const &ctx, int src_first,
                                int dst_first, std::size_t n) {
  for (std::size_t i = 0u; i < n; ++i) {
    auto const s = static_cast<std::size_t>(src_first) + i;
    auto const d = static_cast<std::size_t>(dst_first) + i;
    auto const *fs = ctx.force + s * ctx.force_row_stride;
    auto *fd = ctx.force + d * ctx.force_row_stride;
    fd[0] += fs[0];
    fd[ctx.force_comp_stride] += fs[ctx.force_comp_stride];
    fd[2u * ctx.force_comp_stride] += fs[2u * ctx.force_comp_stride];
#ifdef ESPRESSO_ROTATION
    auto const *ts = ctx.torque + s * ctx.torque_row_stride;
    auto *td = ctx.torque + d * ctx.torque_row_stride;
    td[0] += ts[0];
    td[ctx.torque_comp_stride] += ts[ctx.torque_comp_stride];
    td[2u * ctx.torque_comp_stride] += ts[2u * ctx.torque_comp_stride];
#endif
  }
}

/** @brief Columnar LOCL fast path. Returns false on a broken precondition
 *  (so the caller falls back to the per-particle path).
 *
 *  All src/dst row ranges are resolved (reading only each list's first
 *  particle, no per-particle scan) BEFORE any transfer, so the FORCE `+=`
 *  transfer can never partially apply and then fall back (which would
 *  double-count on the generic retry). */
static bool columnar_cell_cell_transfer(GhostCommunication const &ghost_comm,
                                        BoxGeometry const &box_geo,
                                        unsigned int data_parts) {
  auto const offset = ghost_comm.part_lists.size() / 2;
  // Resolve every src/dst range pair; bail before mutating on break.
  std::vector<std::pair<int, int>> pairs; // (src_first, dst_first); n implicit
  std::vector<std::size_t> sizes;
  pairs.reserve(offset);
  sizes.reserve(offset);
  for (std::size_t pl = 0; pl < offset; pl++) {
    auto *src_list = ghost_comm.part_lists[pl];
    auto *dst_list = ghost_comm.part_lists[pl + offset];
    assert(src_list->rows().size() == dst_list->rows().size());
    if (src_list->rows().empty()) {
      pairs.emplace_back(0, 0);
      sizes.push_back(0u);
      continue;
    }
    auto const src_range = contiguous_store_rows(*src_list);
    auto const dst_range = contiguous_store_rows(*dst_list);
    if (not src_range.has_value() or not dst_range.has_value()) {
      return false;
    }
    pairs.emplace_back(src_range->first, dst_range->first);
    sizes.push_back(src_range->second);
  }
  // Apply the transfers. data_parts is a single fixed value for the
  // whole columnar-eligible call (POSITION xor FORCE), so build the one
  // plain-value row context the branch needs rather than a per-branch
  // std::optional. Using the value type directly (like the
  // make_*_row_context_unchecked callers in columnar_put_recv_buffer /
  // columnar_add_forces_from_recv_buffer) keeps gcc-14 from mis-tracking the
  // std::optional payload as maybe-uninitialized (the optional-payload
  // false-positive family).
  auto &store = *active_particle_store();
  if (data_parts == GHOSTTRANS_POSITION) {
    auto const position_ctx = make_position_row_context_unchecked(store);
    for (std::size_t pl = 0; pl < offset; pl++) {
      auto const n = sizes[pl];
      if (n == 0u) {
        continue;
      }
      auto const [src_first, dst_first] = pairs[pl];
      locl_transfer_position(position_ctx, src_first, dst_first, n, box_geo,
                             ghost_comm.shift);
    }
  } else {
    auto const force_ctx = make_force_row_context(store);
    for (std::size_t pl = 0; pl < offset; pl++) {
      auto const n = sizes[pl];
      if (n == 0u) {
        continue;
      }
      auto const [src_first, dst_first] = pairs[pl];
      locl_transfer_force(force_ctx, src_first, dst_first, n);
    }
  }
  return true;
}

static void cell_cell_transfer(GhostCommunication const &ghost_comm,
                               BoxGeometry const &box_geo,
                               unsigned int data_parts) {
  /* Columnar fast path for the two per-step comms on a clean store. */
  if (columnar_eligible(data_parts) and
      columnar_cell_cell_transfer(ghost_comm, box_geo, data_parts)) {
    return;
  }

  CommBuf buffer;
  if (!(data_parts & GHOSTTRANS_PARTNUM)) {
    buffer.resize(calc_transmit_size(box_geo, data_parts));
  }
  /* Mixed-parts POSITION fast path (see prepare_send_buffer): both source and
   * destination particles read/write pos/image/quat/pos_last through the
   * (clean) store columns by their own row. nullptr => accessor fallback. */
  auto const pos_ctx = (data_parts & GHOSTTRANS_POSITION)
                           ? make_position_row_context()
                           : std::nullopt;
  auto const *pos_ctx_ptr = pos_ctx.has_value() ? &(*pos_ctx) : nullptr;
  /* Mixed-parts MOMENTUM fast path (see prepare_send_buffer): both source and
   * destination particles read/write velocity/omega through the (clean) store
   * columns by their own row. nullptr => accessor fallback. */
  auto const mom_ctx = (data_parts & GHOSTTRANS_MOMENTUM)
                           ? make_momentum_row_context()
                           : std::nullopt;
  auto const *mom_ctx_ptr = mom_ctx.has_value() ? &(*mom_ctx) : nullptr;
  /* Mixed-parts PROPERTIES fast path (see prepare_send_buffer): both source and
   * destination particles read/write the parameter columns/POD sidecars through
   * the (clean) store by their own row. nullptr => accessor fallback. */
  auto const param_ctx = (data_parts & GHOSTTRANS_PROPRTS)
                             ? make_parameter_row_context()
                             : std::nullopt;
  auto const *param_ctx_ptr = param_ctx.has_value() ? &(*param_ctx) : nullptr;
  /* transfer data */
  auto const offset = ghost_comm.part_lists.size() / 2;
  for (std::size_t pl = 0; pl < offset; pl++) {
    auto *src_list = ghost_comm.part_lists[pl];
    auto *dst_list = ghost_comm.part_lists[pl + offset];

    if (data_parts & GHOSTTRANS_PARTNUM) {
      // Use the total count (committed rows + staged) of the source: within a
      // single PARTNUM pass a cell resized earlier as a ghost destination has
      // its new ghosts STAGED (not yet committed to rows), and a downstream
      // ghost layer reading it as a source must see that pending size.
      CellParticleStorage::resize_ghost_storage(*dst_list, src_list->size());
    } else {
      // Contiguous store-row ranges; clean store, so index directly.
      auto const src_offset = src_list->offset();
      auto const dst_offset = dst_list->offset();
      auto const n_rows = src_list->count();
      assert(n_rows == dst_list->count());
      auto &store = *active_particle_store();

      // Reuse two cached views across the row loop, REBOUND per row via
      // attach_to_store (two handle-field writes) instead of constructing a
      // fresh Particle per row per side. This loop runs on every ghost update.
      Particle p1, p2;
      for (std::size_t i = 0; i < n_rows; i++) {
        auto ar_out = Utils::MemcpyOArchive{buffer.make_span()};
        auto ar_in = Utils::MemcpyIArchive{buffer.make_span()};
        // Views over the source and destination store rows.
        p1.attach_to_store(store, static_cast<int>(src_offset + i));
        p2.attach_to_store(store, static_cast<int>(dst_offset + i));
        serialize_and_reduce(ar_out, p1, data_parts, ReductionPolicy::UPDATE,
                             SerializationDirection::SAVE, box_geo,
                             &ghost_comm.shift, pos_ctx_ptr, mom_ctx_ptr,
                             param_ctx_ptr);
        serialize_and_reduce(ar_in, p2, data_parts, ReductionPolicy::UPDATE,
                             SerializationDirection::LOAD, box_geo, nullptr,
                             pos_ctx_ptr, mom_ctx_ptr, param_ctx_ptr);
        if (data_parts & GHOSTTRANS_BONDS) {
          p2.bonds() = p1.bonds();
        }
      }
    }
  }
}

static bool is_send_op(int comm_type, int node, int this_node) {
  return ((comm_type == GHOST_SEND) || (comm_type == GHOST_RDCE) ||
          (comm_type == GHOST_BCST && node == this_node));
}

static bool is_recv_op(int comm_type, int node, int this_node) {
  return ((comm_type == GHOST_RECV) ||
          (comm_type == GHOST_BCST && node != this_node) ||
          (comm_type == GHOST_RDCE && node == this_node));
}

static bool is_prefetchable(GhostCommunication const &ghost_comm,
                            int this_node) {
  int const comm_type = ghost_comm.type & GHOST_JOBMASK;
  int const prefetch = ghost_comm.type & GHOST_PREFETCH;
  int const node = ghost_comm.node;
  return is_send_op(comm_type, node, this_node) && prefetch;
}

static bool is_poststorable(GhostCommunication const &ghost_comm,
                            int this_node) {
  int const comm_type = ghost_comm.type & GHOST_JOBMASK;
  int const poststore = ghost_comm.type & GHOST_PSTSTORE;
  int const node = ghost_comm.node;
  return is_recv_op(comm_type, node, this_node) && poststore;
}

void ghost_communicator(GhostCommunicator const &gcr,
                        BoxGeometry const &box_geo, unsigned int data_parts) {
  if (GHOSTTRANS_NONE == data_parts)
    return;

  static CommBuf send_buffer, recv_buffer;

  auto const &comm = gcr.mpi_comm;

  for (auto cit = gcr.communications.cbegin(); cit != gcr.communications.cend();
       ++cit) {
    auto const &ghost_comm = *cit;
    int const comm_type = ghost_comm.type & GHOST_JOBMASK;

    if (comm_type == GHOST_LOCL) {
      cell_cell_transfer(ghost_comm, box_geo, data_parts);
      continue;
    }

    int const prefetch = ghost_comm.type & GHOST_PREFETCH;
    int const poststore = ghost_comm.type & GHOST_PSTSTORE;
    int const node = ghost_comm.node;

    /* prepare send buffer if necessary */
    if (is_send_op(comm_type, node, comm.rank())) {
      /* ok, we send this step, prepare send buffer if not yet done */
      if (!prefetch) {
        prepare_send_buffer(send_buffer, ghost_comm, box_geo, data_parts);
      }
      // Check prefetched send buffers (must also hold for buffers allocated
      // in the previous lines.)
      assert(send_buffer.size() ==
             calc_transmit_size(ghost_comm, box_geo, data_parts));
    } else if (prefetch) {
      /* we do not send this time, let's look for a prefetch */
      auto prefetch_ghost_comm =
          std::find_if(std::next(cit), gcr.communications.cend(),
                       [this_node = comm.rank()](auto const &other_ghost_comm) {
                         return is_prefetchable(other_ghost_comm, this_node);
                       });

      if (prefetch_ghost_comm != gcr.communications.end())
        prepare_send_buffer(send_buffer, *prefetch_ghost_comm, box_geo,
                            data_parts);
    }

    /* recv buffer for recv and multinode operations to this node */
    if (is_recv_op(comm_type, node, comm.rank()))
      prepare_recv_buffer(recv_buffer, ghost_comm, box_geo, data_parts);

    /* transfer data */
    // Use two send/recvs in order to avoid having to serialize CommBuf
    // (which consists of already serialized data).
    switch (comm_type) {
    case GHOST_RECV:
      comm.recv(node, REQ_GHOST_SEND, recv_buffer.data(),
                static_cast<int>(recv_buffer.size()));
      comm.recv(node, REQ_GHOST_SEND, recv_buffer.bonds());
      break;
    case GHOST_SEND:
      comm.send(node, REQ_GHOST_SEND, send_buffer.data(),
                static_cast<int>(send_buffer.size()));
      comm.send(node, REQ_GHOST_SEND, send_buffer.bonds());
      break;
    case GHOST_BCST:
      if (node == comm.rank()) {
        boost::mpi::broadcast(comm, send_buffer.data(),
                              static_cast<int>(send_buffer.size()), node);
        boost::mpi::broadcast(comm, send_buffer.bonds(), node);
      } else {
        boost::mpi::broadcast(comm, recv_buffer.data(),
                              static_cast<int>(recv_buffer.size()), node);
        boost::mpi::broadcast(comm, recv_buffer.bonds(), node);
      }
      break;
    case GHOST_RDCE:
      if (node == comm.rank())
        boost::mpi::reduce(
            comm, reinterpret_cast<double *>(send_buffer.data()),
            static_cast<int>(send_buffer.size() / sizeof(double)),
            reinterpret_cast<double *>(recv_buffer.data()), std::plus<double>{},
            node);
      else
        boost::mpi::reduce(
            comm, reinterpret_cast<double *>(send_buffer.data()),
            static_cast<int>(send_buffer.size() / sizeof(double)),
            std::plus<double>{}, node);
      break;
    }

    // recv op; write back data directly, if no PSTSTORE delay is requested.
    if (is_recv_op(comm_type, node, comm.rank())) {
      if (!poststore) {
        /* forces have to be added, the rest overwritten. Exception is RDCE,
         * where the addition is integrated into the communication. */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(recv_buffer, ghost_comm);
#ifdef ESPRESSO_BOND_CONSTRAINT
        else if (data_parts == GHOSTTRANS_RATTLE && comm_type != GHOST_RDCE)
          add_rattle_correction_from_recv_buffer(recv_buffer, ghost_comm);
#endif
        else
          put_recv_buffer(recv_buffer, ghost_comm, box_geo, data_parts);
      }
    } else if (poststore) {
      /* send op; write back delayed data from last recv, when this was a
       * prefetch send. */
      /* find previous action where we recv and which has PSTSTORE set */
      auto poststore_ghost_comm = std::find_if(
          std::make_reverse_iterator(cit), gcr.communications.crend(),
          [this_node = comm.rank()](auto const &other_ghost_comm) {
            return is_poststorable(other_ghost_comm, this_node);
          });

      if (poststore_ghost_comm != gcr.communications.rend()) {
        assert(recv_buffer.size() ==
               calc_transmit_size(*poststore_ghost_comm, box_geo, data_parts));
        /* as above */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(recv_buffer, *poststore_ghost_comm);
#ifdef ESPRESSO_BOND_CONSTRAINT
        else if (data_parts == GHOSTTRANS_RATTLE && comm_type != GHOST_RDCE)
          add_rattle_correction_from_recv_buffer(recv_buffer,
                                                 *poststore_ghost_comm);
#endif
        else
          put_recv_buffer(recv_buffer, *poststore_ghost_comm, box_geo,
                          data_parts);
      }
    }
  }
}
