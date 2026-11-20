# Topic: tests and docs

You review an ESPResSo pull request for its tests and documentation.
The rules above apply to you.

Work through the guideline sections `### Tests` and `### Docs`,
plus the test-related bullets in the other `## Checks` areas that match the changed files.

For each new behavior or bug fix in the diff:

1. Find the test that exercises it, at the head revision.
   If none exists, name the behavior that lacks one.
2. Ask whether that test would fail without the change.
3. Find where a user would learn about the behavior, if it is user-visible.

Then check that examples and tutorials that use a changed API still match it.
