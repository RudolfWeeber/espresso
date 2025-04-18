
#include <array>
#include <atomic>
#include <cstddef>
#include <iostream>

#include "Vector.hpp"

namespace Utils {

template <typename T, std::size_t N> class AtomicVector {
  static_assert(std::atomic<T>::is_always_lock_free,
                "T must be lock-free atomic");

public:
  AtomicVector(const AtomicVector &other) { set(other.get()); }
  AtomicVector &operator=(const AtomicVector &other) {
    set(other.get());
    return *this;
  }

  AtomicVector() { clear(); }

  void clear() {
    for (std::size_t i = 0; i < N; ++i) {
      data_[i].store(0, std::memory_order_relaxed);
    }
  }

  void set(const Vector<T, N> &other) {
    for (std::size_t i = 0; i < N; ++i) {
      data_[i].store(other[i], std::memory_order_relaxed);
    }
  }

  void add(const Vector<T, N> &other) {
    for (std::size_t i = 0; i < N; ++i) {
      data_[i].fetch_add(other[i], std::memory_order_relaxed);
    }
  }

  Vector<T, N> get() const {
    Vector<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
      result[i] = data_[i].load(std::memory_order_relaxed);
    }
    return result;
  }

private:
  std::atomic<T> data_[N];
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int /*version*/) {
    if constexpr (Archive::is_saving::value) {
      Vector<T, N> snapshot = get();
      ar & snapshot;
    } else {
      Utils::Vector<T, N> temp;
      ar & temp;
      set(temp); // std::array<T,N>::data() returns T*
    }
  }
};

} // namespace Utils
