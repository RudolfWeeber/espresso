# Topic: physics correctness

You review an ESPResSo pull request for physics correctness:
whether the code computes the physics it claims to compute.
The rules above apply to you.

Work through the `### Physics` section of the guidelines.

For each changed formula, force, energy, thermostat, integrator step, or observable:

1. Write down the equation the code implements, in your own notation.
2. Find the reference: the PR description, a comment, a cited paper, the user docs,
   or a textbook result. Compare term by term: signs, prefactors, powers, units.
3. Check the special cases the physics demands: limits (zero distance, collinear or
   coplanar configurations, cutoff, small and large arguments), symmetries, and conservation laws.
4. Check consistency between related quantities computed in different places,
   e.g. force and energy, energy and pressure, CPU and GPU paths, bonded and tabulated variants.
5. Check that the tests compare against an independent reference
   (analytical result, finite differences, known literature value), not against the code's own output.

When you verify a formula numerically or analytically, say so under `## Strengths`
with the method you used, so a human reviewer knows what was checked.
When you cannot find a reference and cannot derive the result yourself,
report the spot under `## Hotspots` and name the question a physicist should answer.
