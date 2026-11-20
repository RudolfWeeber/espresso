# Risk levels

The risk rating tells a human reviewer how much attention a PR needs.
It has two parts:

1. A **floor** from the path table: the highest level any changed file reaches.
2. An optional **raise** after verification, when verified findings show wider effects.

The rating never goes below the floor.
The report states the level, one sentence on what could break and who would notice,
and the files that set the floor.

## Path table

Patterns are globs relative to the repository root.
For each changed file, the first matching row wins.
The team edits this table by PR when the code layout changes.

| Pattern | Level |
|---|---|
| `src/core/integrate*` | High |
| `src/core/integrators/**` | High |
| `src/core/cell_system/**` | High |
| `src/core/cells.*` | High |
| `src/core/Particle*.hpp` | High |
| `src/core/ghosts*` | High |
| `src/core/ghosts/**` | High |
| `src/core/communication.*` | High |
| `src/core/MpiCallbacks.hpp` | High |
| `src/core/forces*` | High |
| `src/core/system/**` | High |
| `src/script_interface/*.hpp` | High |
| `src/script_interface/auto_parameters/**` | High |
| `src/utils/**` | High |
| `CMakeLists.txt` | High |
| `cmake/**` | High |
| `.gitlab-ci.yml` | High |
| `.github/**` | High |
| `maintainer/CI/**` | High |
| `src/**` | Medium |
| `doc/**` | Low |
| `testsuite/**` | Low |
| `samples/**` | Low |
| `maintainer/**` | Low |
| `*.md` | Low |
| `**` | Medium |

`*` matches within one path component; `**` matches any number of components.
So `CMakeLists.txt` and `*.md` match only top-level files;
files in subdirectories fall through to the row for their directory.

## Level definitions

### Low

Docs, tests, samples, or code confined to one opt-in feature.
Few users are affected, and mistakes are easy to spot.

### Medium

One physics method, interaction, or script-interface class;
or user-visible changes (API, defaults, numerical results) within one method.

### High

Shared infrastructure that every simulation passes through
(integrator, cell system, particle data, ghost and MPI communication,
script-interface base classes, build system, CI);
or anything that silently changes results, parallel correctness, or performance.

## Raising the rating

After verification, the coordinator raises the rating to High when a confirmed finding
shows that the PR itself causes one of these
(defects the PR leaves in place do not count):

* results change silently, without an error or a documented change;
* parallel correctness is affected (rank-dependent results, deadlocks, missing ghost updates);
* performance of code paths shared by many simulations changes.

A Low-floor PR whose change reaches beyond one method rises to Medium.
The report names the item that caused the raise.
The coordinator never lowers the floor.
