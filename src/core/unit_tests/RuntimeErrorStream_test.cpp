/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

/* Unit tests for the ErrorHandling::RuntimeErrorStream class. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE RuntimeErrorStream test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "error_handling/RuntimeError.hpp"
#include "error_handling/RuntimeErrorCollector.hpp"
#include "error_handling/RuntimeErrorStream.hpp"

#include <boost/mpi.hpp>

#include <cstddef>
#include <cstring>
#include <new>
#include <string>

using ErrorHandling::RuntimeError;
using ErrorHandling::RuntimeErrorCollector;
using ErrorHandling::RuntimeErrorStream;
using ErrorLevel = ErrorHandling::RuntimeError::ErrorLevel;

/*
 * Regression test for the RuntimeErrorStream copy constructor (bug-sweep #58).
 *
 * The user-defined copy constructor must copy every member, including the
 * severity m_level. It used to omit m_level from its member initializer list,
 * so the copy's level was left whatever happened to be in the storage; the
 * destructor reads m_level and forwards it to RuntimeErrorCollector::message().
 *
 * Test mechanics:
 *  - The copy is placement-new constructed into a raw stack buffer. Before
 *    construction we deterministically write a *valid but non-ERROR* byte
 *    pattern (WARNING == 0) over the whole buffer. The copy constructor does
 *    not touch m_level on the unfixed branch, so on unfixed code the copy's
 *    m_level reads back as WARNING, while on fixed code it is overwritten with
 *    the source's ERROR. The written-then-read byte is fully defined (no
 *    reliance on indeterminate-value semantics).
 *  - The copy constructor copies the stream content via
 *    `m_buff << rhs.m_buff.rdbuf()`. Because the source is an *output*
 *    stringstream, its get area is empty, so this copies nothing and the
 *    copy's message text is always empty. We therefore identify the copy's
 *    record not by its (empty) text but by a unique file name that the copy
 *    constructor *does* propagate (m_file is copied). The source's record
 *    carries the same unique file plus a non-empty marker text, so the copy
 *    is the record whose file matches and whose text is empty.
 *
 * On unfixed code the copy's record level is WARNING (mismatch -> failure);
 * on fixed code it is ERROR (match -> pass).
 */
BOOST_AUTO_TEST_CASE(copy_ctor_preserves_level) {
  boost::mpi::communicator world;

  RuntimeErrorCollector rec(world);

  std::string const unique_file = "bug58_unique_file_marker";
  std::string const source_text = "bug58_source_text";

  {
    // Raw aligned storage, deterministically pre-filled with bytes that
    // encode the valid ErrorLevel WARNING (== 0). On the unfixed copy
    // constructor m_level is never written, so it reads back as WARNING.
    alignas(
        RuntimeErrorStream) unsigned char storage[sizeof(RuntimeErrorStream)];
    static_assert(static_cast<int>(ErrorLevel::WARNING) == 0,
                  "test assumes WARNING maps to all-zero bytes");
    std::memset(storage, 0, sizeof(storage));

    // Source stream: ERROR level, unique file, non-empty marker text.
    RuntimeErrorStream source(rec, ErrorLevel::ERROR, unique_file, 42,
                              "Test_function");
    source << source_text;

    // Copy-construct into the pre-filled storage. This invokes the copy
    // constructor on a genuine lvalue (never elided), exercising exactly
    // the code path that dropped m_level on the unfixed branch.
    auto *copy = new (storage) RuntimeErrorStream(source);

    // Manually invoke the destructor: it reads m_level and records a
    // message (empty text, unique file). On the unfixed branch m_level is
    // still WARNING; on the fixed branch it is the source's ERROR.
    copy->~RuntimeErrorStream();
  }

  world.barrier();

  if (world.rank() == 0) {
    auto const results = rec.gather();

    // The copy's record is the one with the unique file and empty text
    // (the copy constructor copies m_file but, due to the output-stream
    // rdbuf, never carries over the message text).
    int copy_count = 0;
    ErrorLevel copy_level{};
    for (auto const &err : results) {
      if (err.file() == unique_file && err.what().empty()) {
        ++copy_count;
        copy_level = err.level();
      }
    }

    BOOST_REQUIRE_EQUAL(copy_count, 1);
    // The copy must preserve the source's ERROR level. On the unfixed copy
    // constructor m_level keeps the pre-filled WARNING value, so this fails.
    BOOST_CHECK_EQUAL(static_cast<int>(copy_level),
                      static_cast<int>(ErrorLevel::ERROR));
  } else {
    rec.gather_local();
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
