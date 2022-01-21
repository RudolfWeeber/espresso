#ifndef LEES_EDWARDS_PROTOCOL_HPP
#define LEES_EDWARDS_PROTOCOL_HPP

#include "integrate.hpp"

#include <utils/Vector.hpp>

#include <boost/variant.hpp>

#include <cmath>

/** \file lees_edwards.hpp
 *
 */

namespace LeesEdwards {

// Protocols determining shear rate and positional offset as a function of time
/** Lees Edwards protocol for un-altered periodic boundary conditions */
struct Off {
  double shear_velocity(double time) const { return 0.; }
  double pos_offset(double time) const { return 0.; }
};

/** Lees-Edwards protocol for linear shearing */
struct LinearShear {
  LinearShear() : m_initial_pos_offset{0}, m_shear_velocity{0}, m_time_0{0} {}
  LinearShear(double initial_offset, double shear_velocity, double time_0)
      : m_initial_pos_offset{initial_offset},
        m_shear_velocity{shear_velocity}, m_time_0{time_0} {}
  double shear_velocity(double time) const { return m_shear_velocity; }
  double pos_offset(double time) const {
    return m_initial_pos_offset + (time - m_time_0) * m_shear_velocity;
  }
  double m_initial_pos_offset;
  double m_shear_velocity;
  double m_time_0;
};

/** Lees-Edwards protocol for osciallatory shearing */
struct OscillatoryShear {
  OscillatoryShear() : m_amplitude{0}, m_omega{0}, m_time_0{0} {}
  OscillatoryShear(double amplitude, double omega, double time_0)
      : m_amplitude{amplitude}, m_omega{omega}, m_time_0{time_0} {}
  double pos_offset(double time) const {
    return m_amplitude * std::sin(m_omega * (time - m_time_0));
  }
  double shear_velocity(double time) const {
    return m_omega * m_amplitude * std::cos(m_omega * (time - m_time_0));
  }
  double m_amplitude;
  double m_omega;
  double m_time_0;
};

/** Type which holds the currently active protocol */
using ActiveProtocol = boost::variant<Off, LinearShear, OscillatoryShear>;

class PosOffsetGetter : public boost::static_visitor<double> {
public:
  PosOffsetGetter(double time) : m_time{time} {}
  template <typename T> double operator()(const T &protocol) const {
    return protocol.pos_offset(m_time);
  }

private:
  double m_time;
};

inline double get_pos_offset(double time, const ActiveProtocol &protocol) {
  return boost::apply_visitor(PosOffsetGetter(time), protocol);
}

/** Visitor to get shear velocity from the Lees-Edwards protocol */
class ShearVelocityGetter : public boost::static_visitor<double> {
public:
  ShearVelocityGetter(double time) : m_time{time} {}
  template <typename T> double operator()(const T &protocol) const {
    return protocol.shear_velocity(m_time);
  }

private:
  double m_time;
};

/** Calculation of current velocity*/
inline double get_shear_velocity(double time, const ActiveProtocol &protocol) {
  return boost::apply_visitor(ShearVelocityGetter(time), protocol);
}

} // namespace LeesEdwards

#endif
