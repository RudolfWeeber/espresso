/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifdef SHARED_MEMORY_PARALLELISM

#include <Cabana_VerletList.hpp>

namespace Cabana {
// ONLY FOR 2D LAYOUT, OTHERWISE NEIGHBOR LIST INTERFACE IMPLEMENTATION WILL
// CAUSE PROBLEMS (NOT IMPLEMENTED)
template <class MemorySpace, class AlgorithmTag, class LayoutTag,
          class BuildTag = TeamVectorOpTag>
class CustomVerletList
    : public VerletList<MemorySpace, AlgorithmTag, LayoutTag, BuildTag> {
public:
  using Base = VerletList<MemorySpace, AlgorithmTag, LayoutTag, BuildTag>;

  // Default constructor
  CustomVerletList() : Base() {}

  // Custom constructor
  template <class PositionSlice>
  CustomVerletList(PositionSlice x, const std::size_t begin,
                   const std::size_t end, const std::size_t max_neigh,
                   const std::size_t thread_number) {
    initializeData(x.size(), max_neigh, thread_number);
  }
  virtual ~CustomVerletList() {};

private:
  Kokkos::View<int *, MemorySpace> max_thread;
  Kokkos::View<int **, MemorySpace> counts_thread;
  Kokkos::View<int ***, MemorySpace> neighbors_thread;

public:
  Kokkos::View<int *, MemorySpace> counts;
  Kokkos::View<int **, MemorySpace> neighbors;

  // Method to initialize _data without filling neighbors
  KOKKOS_INLINE_FUNCTION
  void initializeData(const std::size_t num_particles,
                      const std::size_t max_neigh,
                      const std::size_t thread_number) {
    counts = Kokkos::View<int *, MemorySpace>("num_neighbors", num_particles);
    neighbors = Kokkos::View<int **, MemorySpace>(
        Kokkos::ViewAllocateWithoutInitializing("neighbors"), num_particles,
        max_neigh);
    max_thread = Kokkos::View<int *, MemorySpace>("max_thread", thread_number);
    for (int tid = 0; tid < thread_number; ++tid) {
      max_thread(tid) = 1;
    }
  }

  // Method to dynamically expand the size of max_neighbors
  // This function may be vaiolated Kokkos's rule.
  // Kokkos::View should not be created in Kokkos::parallel.
  // However, this function is called from addNeighbor used
  // in the Kokkos::parallel.
  KOKKOS_INLINE_FUNCTION
  void expandMaxNeighbors(const std::size_t new_max_neigh) {
    // Create a new view with the larger size
    Kokkos::View<int **, MemorySpace> new_neighbors(
        Kokkos::ViewAllocateWithoutInitializing("neighbors"),
        neighbors.extent(0), new_max_neigh);

    // Copy existing data to the new view
    Kokkos::parallel_for("copy_neighbors", neighbors.extent(0),
                         [=, this](const int i) {
                           for (std::size_t j = 0; j < counts(i); ++j) {
                             new_neighbors(i, j) = neighbors(i, j);
                           }
                         });

    // Replace the old view with the new view
    neighbors = new_neighbors;
  }

  // Method to add a neighbor
  KOKKOS_INLINE_FUNCTION
  void addNeighbor(const int tid, int pid, int nid) {
    if (counts(pid) + 1 > max_thread(tid))
      std::swap(pid, nid);
    std::size_t count = Kokkos::atomic_fetch_add(&counts(pid), 1);
    if (count >= neighbors.extent(1)) {
      // expandMaxNeighbors(neighbors.extent(1) * 2);
      throw std::runtime_error(
          "Number of count is larger than VerletList size.");
    }
    neighbors(pid, count) = nid;
    if (counts(pid) > max_thread(tid))
      max_thread(tid) = counts(pid);
  }

  // Thread safe but non atomic method to add a neighbor
  KOKKOS_INLINE_FUNCTION
  void addNeighborNonAtomic(const int tid, int pid, int nid) {
    if (counts(pid) + 1 > max_thread(tid))
      std::swap(pid, nid);
    std::size_t count = counts(pid);
    counts(pid) += 1;
    if (count >= neighbors.extent(1)) {
      throw std::runtime_error(
          "Number of count in one thread is larger than VerletList size.");
    }
    neighbors(pid, count) = nid;
    if (counts(pid) > max_thread(tid))
      max_thread(tid) = counts(pid);
  }

  // Reduction of counts and neighbor in all threads
  KOKKOS_INLINE_FUNCTION
  void reduction() {
    // Kokkos::RangePolicy<execution_space> policy(0, counts.extent(0));
    int thread_number = counts_thread.extent(0);
    Kokkos::parallel_for(
        "reduction_neighbor", counts.extent(0), [&](const int pid) {
          counts(pid) = 0;
          for (int tid = 0; tid < thread_number; ++tid) {
            std::size_t offset = counts(pid);
            counts(pid) += counts_thread(tid, pid);
            if (counts(pid) >= neighbors.extent(1)) {
              throw std::runtime_error(
                  "Number of count is larger than VerletList size.");
            }
            for (int cid = 0; cid < counts_thread(tid, pid); ++cid) {
              neighbors(pid, offset + cid) = neighbors_thread(tid, pid, cid);
            }
          }
        });
    Kokkos::fence();
  }

  // Find max counts
  KOKKOS_INLINE_FUNCTION
  std::size_t get_max_counts() {
    std::size_t max_counts = 0;
    std::size_t ave_counts = 0;
    for (int pid = 0; pid < counts.extent(0); ++pid) {
      if (max_counts < counts(pid))
        max_counts = counts(pid);
      ave_counts += counts(pid);
    }
    if (counts.extent(0) != 0) {
      std::cout << "max:" << max_counts
                << " ave:" << ave_counts / counts.extent(0) << std::endl;
    }
    return max_counts;
  }
};

template <class MemorySpace, class AlgorithmTag, class BuildTag>
class NeighborList<
    CustomVerletList<MemorySpace, AlgorithmTag, VerletLayout2D, BuildTag>> {
public:
  //! Kokkos memory space.
  using memory_space = MemorySpace;
  //! Neighbor list type.
  using list_type =
      CustomVerletList<MemorySpace, AlgorithmTag, VerletLayout2D, BuildTag>;

  //! Get the total number of neighbors across all particles.
  KOKKOS_INLINE_FUNCTION
  static std::size_t totalNeighbor(const list_type &list) {
    std::size_t total_n = 0;
    std::size_t num_p = list.counts.size();
    for (std::size_t i = 0; i < num_p; ++i)
      total_n += list.counts(i);
    return total_n;
  }

  //! Get the maximum number of neighbors per particle.
  KOKKOS_INLINE_FUNCTION
  static std::size_t maxNeighbor(const list_type &list) {
    // Stored during neighbor search.
    return list.max_n;
  }

  //! Get the number of neighbors for a given particle index.
  KOKKOS_INLINE_FUNCTION
  static std::size_t numNeighbor(const list_type &list,
                                 const std::size_t particle_index) {
    return list.counts(particle_index);
  }

  //! Get the id for a neighbor for a given particle index and the index of
  //! the neighbor relative to the particle.
  KOKKOS_INLINE_FUNCTION
  static std::size_t getNeighbor(const list_type &list,
                                 const std::size_t particle_index,
                                 const std::size_t count) {
    return list.neighbors(particle_index, count);
  }
};

} // namespace Cabana

#endif
