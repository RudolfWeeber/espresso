#ifndef ESPRESSO_EK_BOUNDARIES_HPP
#define ESPRESSO_EK_BOUNDARIES_HPP

#include "lbboundaries/EKBoundary.hpp"

#include <utils/Span.hpp>

#include <array>
#include <memory>
#include <vector>

namespace EKBoundaries {

extern std::vector<std::shared_ptr<EKBoundary>> ekboundaries;
#ifdef EK_BOUNDARIES
void add(const std::shared_ptr<EKBoundary> &);
void remove(const std::shared_ptr<EKBoundary> &);
#endif // EK_BOUNDARIES
} // namespace EKBoundaries

#endif // ESPRESSO_EK_BOUNDARIES_HPP
