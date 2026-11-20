# Rules for every reviewer

You review one ESPResSo pull request.
Your topic or topics follow after these rules.

## Context

{{CONTEXT}}

Inspect the change with:

```
{{DIFF_CMD}}
```

Files in your scope:

{{FILES}}

## How to work

* Read `{{GUIDELINES_PATH}}` first. Cite its headings in the `rule:` field.
  Use `defect` for bugs that no heading covers.
* Read the surrounding code as needed: callers, callees, base classes, similar existing code.
  Read files at the head revision with `git show <head>:<path>`.
* The review is read-only. Do not check out, stash, reset, commit or edit anything.
* Do all the work yourself. Do not dispatch subagents.
* Report only issues in your topics. Other topics have their own reviewers.
* Skip anything `clang-format`, `autopep8`, `pylint` or `clang-tidy` would report.
* Report only issues you can point to in the code. A finding needs a location and a concrete effect.
* Judge the code as of the head commit. Do not use commits after the head as evidence.
* A finding is about something this PR introduces or changes.
  A defect in code the PR does not touch, or a gap the PR merely leaves in place,
  goes under `## Hotspots`, not under findings.
* A behavior change that the PR does not announce in its description is a defect, not a documentation issue.
  This includes changed boundary conventions (inclusive vs. exclusive limits), defaults, units, and signs.
* Grade severity by effect on users and maintainers:
  * Critical: wrong results, crashes, deadlocks, data loss, broken build.
  * Important: a real defect in a less common path, a silent behavior change,
    a missing test for new behavior, a clear violation of a goal in the guidelines.
  * Minor: small improvements to clarity or consistency.
* Refer to other findings by title, not by number.

## Output format

One block per finding, numbered from F1:

```
### F<n>: <short title>
- severity: Critical | Important | Minor
- location: <path>:<line>[-<line>][, <path>:<line>...]   (primary location first)
- rule: <heading from the guidelines, or "defect">
- problem: <what is wrong>
- impact: <why it matters>
- fix: <suggestion, or "obvious">
```

Then:

```
## Hotspots
- <path>[:<line>] — <why a human should look here, even without a finding>

## Strengths
- <specific thing done well, with path>

## Set aside
- <behavior you considered and left out of scope> — <reason>
```

Write `none` under a heading that has no entries.
