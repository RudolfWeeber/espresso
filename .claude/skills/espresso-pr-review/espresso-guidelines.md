# ESPResSo review guidelines

This file is the only source of project rules for the `espresso-pr-review` skill.
Reviewers cite its headings when they report a finding.

How to extend it:
add one bullet per check, phrased as something a reviewer can verify in a diff.
Link the source (wiki page, issue, discussion) where one exists.
Leave out anything `clang-format`, `autopep8`, `pylint` or `clang-tidy` already catch.

Background on the code base and its infrastructure: [`AGENT.md`](../../../AGENT.md).
Human review process and style guide:
[wiki "Code review"](https://github.com/espressomd/espresso/wiki/Code-review).

## Goals

### Understandability

A developer who did not write the code can follow it without asking the author.

* Names say what a thing is or does; one term per concept, consistent with the surrounding code.
* Functions do one thing. Long functions with several levels of nesting are split.
* Non-obvious physics, numerics and parallel logic carry a short comment or a reference to the paper or method.
* Comments describe the current state of the code, not its history.
* The PR description explains why the change is needed and what the user-visible effect is.

### Maintainability

The change is cheap to keep working while the rest of the code evolves.

* No duplicated logic. A second copy of an algorithm, a parameter check or a communication pattern is a finding.
* Dependencies point in the existing direction: Python → script interface → core → utils.
  The core does not include script-interface or Python headers.
* No new global state. New code receives existing globals by reference instead of reaching for them.
* Each commit is self-contained, so it can be reverted or cherry-picked.
* New code paths are covered by tests
  (the wiki "Code review" page sets a codecov patch target of 99 %).
* Optional dependencies are guarded by the existing feature and CMake options, not by ad-hoc macros.

### Use of existing infrastructure

New code builds on what ESPResSo already provides instead of reinventing it.
`AGENT.md` lists the main building blocks.

* Script-interface classes derive from `AutoParameters` and use the existing container types
  (`ObjectList`, `ObjectMap`, `ObjectContainer`) rather than hand-written parameter handling.
* Vector math uses `Utils::Vector` and the helpers in `src/utils`.
* Parallel communication uses `boost::mpi` collectives and the existing ghost and callback mechanisms
  where boost::mpi or those mechanisms cover the operation.
  Raw MPI calls exist on purpose in performance-critical code (FFT, MPI-IO, mesh communication).
* Errors in the core use `runtimeErrorMsg()` or exceptions that the script interface propagates to Python.
* Compile-time features use `src/config/features.def`; build options use `ESPRESSO_BUILD_WITH_*`.
* Tests reuse the helpers in `testsuite/python` (e.g. `unittest_decorators`, `tests_common`)
  and are registered with `python_test(...)` in `testsuite/python/CMakeLists.txt`.
* A reviewer who suspects a reinvention names the existing component and its path.

## Checks

### C++ core

* No `malloc`/`free`, no bare owning pointers; use STL containers and smart pointers.
* No `const_cast` except at C-library boundaries.
* Everything that can be `const` is `const`.
* Prefer standard or Boost.Range algorithms over raw loops where it improves clarity.
* No new macros where a function, template or `constexpr` works.
* No mixing of signed and unsigned integer types.
* Code that runs on all ranks performs the same collective calls in the same order on every rank.
* Rank-local data (local particles, local cells) is not treated as global, and vice versa.
* Ghost particles are updated before they are read, and forces on ghosts are folded back.
* Changes to the particle data, cell system or integrator keep results identical with 1 and several MPI ranks.
* Units and conventions match the surrounding code (ESPResSo is unit-agnostic; no hard-coded physical constants).

### Script interface

* New parameters are registered via `add_parameters()` with correct read-only flags.
* Parameter validation happens once, in one place, with an error message that names the parameter.
* New objects and parameters survive checkpointing
  (`testsuite/python/save_checkpoint.py`, `test_checkpoint.py`).
* Python and C++ class names and parameter names match.

### Python API

* Modules are imported with their namespace (`import espressomd.lb`, not `from espressomd import lb`); no `eval()`.
* User-visible API changes are consistent with existing naming in `espressomd`.
* Breaking changes are intentional, documented, and listed in `CHANGELOG.md`.
* Default values do not change silently.

### Behavior changes

* A behavior change that the PR does not announce in its description is a defect, not a documentation issue.
* This includes boundary conventions (inclusive vs. exclusive limits), defaults, units, signs,
  and which quantities an observable or analysis result includes.

### Tests

* New features and bug fixes come with tests that fail without the change.
* New checks extend the existing test file for the feature instead of adding a new file,
  unless the new file covers a clearly separate concern.
* Tests check physical or numerical results against a reference, not only that code runs.
* Parallel code is tested with more than one MPI rank (`MAX_NUM_PROC` in `python_test`).
* Tests that need optional features use the `unittest_decorators` skip helpers.
* New test files are registered in the matching `CMakeLists.txt`.
* C++ logic that can be tested in isolation has a unit test in `src/core/unit_tests` or the matching `tests` folder.

### Docs

* User-visible changes are documented in `doc/sphinx`.
* New public C++ functions and classes have Doxygen comments.
* Samples and tutorials that use a changed API are updated.
* `CHANGELOG.md` has an entry for user-visible changes.

### Build and CI

* New dependencies are optional unless the team agreed otherwise, and are guarded by an `ESPRESSO_BUILD_WITH_*` option.
* New source files are added to the right CMake target.
* Changes to CI configuration keep all existing jobs meaningful (no silently skipped tests).
