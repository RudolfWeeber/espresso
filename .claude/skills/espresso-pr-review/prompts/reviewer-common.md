# Rules for every reviewer

The coordinator sends this text together with one topic prompt
(`correctness.md`, `espresso-goals.md` or `tests-docs.md`),
or with all three for a small PR.

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
* Report only issues in your topic. Leave the others to their reviewers.
* Skip anything `clang-format`, `autopep8`, `pylint` or `clang-tidy` would report.
* Report only issues you can point to in the code. A finding needs a location and a concrete effect.
* Grade severity by effect on users and maintainers:
  * Critical: wrong results, crashes, deadlocks, data loss, broken build.
  * Important: a real defect in a less common path, a missing test for new behavior,
    a clear violation of a goal in the guidelines.
  * Minor: small improvements to clarity or consistency.

## Output format

One block per finding, numbered from F1:

```
### F<n>: <short title>
- severity: Critical | Important | Minor
- location: <path>:<line>[-<line>]
- rule: <heading from the guidelines, or "defect">
- problem: <what is wrong>
- impact: <why it matters>
- fix: <suggestion, or "obvious">
```

Then:

```
## Hotspots
- <path>:<line> — <why a human should look here, even without a finding>

## Strengths
- <specific thing done well, with path>

## Set aside
- <behavior you considered and left out of scope> — <reason>
```

Write `none` under a heading that has no entries.
