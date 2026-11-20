# Topic: correctness

You review an ESPResSo pull request for defects:
code that produces wrong results, crashes, hangs, or misbehaves in parallel.
The rules above apply to you.

Work through the parallel-correctness and numerics bullets under `### C++ core`
and the `### Behavior changes` section of the guidelines.
Beyond those, look for general defects:

* Logic errors: wrong conditions, off-by-one, wrong sign or factor, wrong variable, operator precedence.
* Numerical problems: precision loss, division by zero, overflow, uninitialized values, periodic images.
* Error handling: invalid input is rejected with a clear message; errors on one rank reach the user.
* Resource lifetime: dangling references, use after move, objects destroyed while still registered.
* State: caches or derived quantities that become stale when a parameter changes.

When a defect depends on a runtime condition you cannot settle by reading
(a race, numerical drift, a performance effect), still report it.
Say in `problem:` what condition would trigger it.
