/*
 * Copyright (C) 2014-2019 The ESPResSo project
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

#ifndef ERROR_HANDLING_RUNTIMEERRORCOLLECTOR_HPP
#define ERROR_HANDLING_RUNTIMEERRORCOLLECTOR_HPP

#include "RuntimeError.hpp"

#include <boost/mpi/communicator.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace ErrorHandling {

class RuntimeErrorCollector {
public:
  explicit RuntimeErrorCollector(boost::mpi::communicator comm);
  ~RuntimeErrorCollector();

  void message(RuntimeError message);
  void message(const RuntimeError &message);
  void message(RuntimeError::ErrorLevel level, const std::string &msg,
               const char *function, const char *file, int line);

  void warning(const std::string &msg, const char *function, const char *file,
               int line);
  void warning(const char *msg, const char *function, const char *file,
               int line);
  void warning(const std::ostringstream &mstr, const char *function,
               const char *file, int line);

  void error(const std::string &msg, const char *function, const char *file,
             int line);
  void error(const char *msg, const char *function, const char *file, int line);
  void error(const std::ostringstream &mstr, const char *function,
             const char *file, int line);

  /**
   * \brief Get the number of all flying messages on all nodes.
   *
   * @return Total number of messages.
   */
  int count() const;

  /**
   * \brief Get the number of messages that have at least severity
   * @p level on this node.
   *
   * @param level Severity filter.
   * @return Number of messages that match the filter.
   */
  int count(RuntimeError::ErrorLevel level);

  /**
   * @brief Reset error messages.
   */
  void clear();

  /**
   * @brief Flush error messages to standard error.
   */
  void flush();

  std::vector<RuntimeError> gather();
  void gather_local();

  const boost::mpi::communicator &comm() const { return m_comm; }

private:
  std::vector<RuntimeError> m_errors;
  boost::mpi::communicator m_comm;
};

} // namespace ErrorHandling

#endif
