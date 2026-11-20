# ESPResSo PR review skill: design

## Intent

The ESPResSo core team wants a Claude Code skill that reviews pull requests.
It addresses three problems at once:

* PRs reach human reviewers with avoidable problems.
* Reviewers spend too long finding their way through large diffs.
* Too few people are available to review.

The skill therefore produces three things:
a risk rating, a reading guide for the human reviewer, and verified findings.
It judges PRs against ESPResSo-specific goals that the team maintains in one editable file.

## Users

* PR authors, who run it before requesting review.
* Maintainers, who run it on other people's PRs.

Both get the same report (decided: one report layout, no modes).

## Goals

* Report starts with a risk rating, then a reading guide, then findings.
* Findings are checked against the ESPResSo goals:
  understandability by humans, maintainability, use of existing infrastructure.
* The team extends goals and checks by editing one Markdown file.
* Output is a local Markdown file.
* Findings are verified against the code before they enter the report.

## Non-goals

* Style issues that `clang-format`, `autopep8`, `pylint` or `clang-tidy` already catch.
* Building or running tests. The review is static.
* Replacing the human review or the CI.

## Anti-goals

* Posting to GitHub without an explicit request from the user.
  Explicit means the user asks in words, e.g. "post this to the PR".
* Noisy reports. False alarms make the team stop trusting the tool.
* Judging by taste instead of by the written guidelines or by concrete defects.
* Modifying the user's checkout (no checkout of the PR branch, no stash, no reset).

## Constraints

* Lives in the repo at `.claude/skills/espresso-pr-review/`, shared with the team.
* Works for any maintainer's Claude Code setup.
  It must not depend on the superpowers plugin being installed; it adapts its ideas.
* Input: a GitHub PR number (fetched with `gh`) or a local branch / commit range.
  Base is the merge-base with the upstream `python` branch unless the user gives one.
* PR refs are fetched without checking them out,
  e.g. `git fetch <upstream> pull/<N>/head:refs/review/pr-<N>`.
* Report path: `.claude/reviews/<pr-N or branch-name>.md`.
  `.claude/reviews/` is added to `.gitignore`.
* `espresso-guidelines.md` is the only source of review rules at runtime.
  It is seeded once by distilling `AGENT.md`, `CONTRIBUTING.md` and the wiki page "Code review";
  it links to them instead of copying long passages.
  Reviewers may read the local `AGENT.md` for infrastructure details.
  The skill never fetches the wiki at runtime (no network dependency). *(Our call.)*

## Approaches considered

* **A. One reviewer subagent.** Cheapest; mirrors superpowers `requesting-code-review`.
  On large PRs one reader skims, and the ESPResSo goals compete with bug hunting.
* **B. Coordinator plus focused reviewers (chosen).**
  Best coverage of avoidable problems.
  The ESPResSo goals get a dedicated reviewer.
  A verification step keeps false alarms out. Costs more tokens.
* **C. Inline, no subagents.** Simplest; runs out of context on large PRs.

The team picked B.

## Design

### Files

```
.claude/skills/espresso-pr-review/
  SKILL.md                  # coordinator procedure, invocation, posting rule
  espresso-guidelines.md    # team-editable goals and checks
  risk-levels.md            # definition of Low / Medium / High
  report-template.md        # layout of the output file
  prompts/
    correctness.md          # reviewer: bugs, MPI/parallel correctness
    espresso-goals.md       # reviewer: understandability, maintainability, infrastructure reuse
    tests-docs.md           # reviewer: tests, docs, samples, changelog
    verify.md               # verifier: confirm or reject one batch of findings
```

Prompts never contain project rules.
They tell reviewers to read `espresso-guidelines.md`.
Editing the guidelines therefore cannot break the procedure.

### Coordinator procedure (SKILL.md)

1. **Collect context.**
   Diff, commit list, PR title and description, linked issues, changed-file list with line counts.
2. **Rate risk** using `risk-levels.md`.
   A path-to-level map gives a floor: the highest level any changed file reaches.
   Record the files that set it.
   The coordinator raises the rating after step 5 if verified findings show
   silent result changes, parallel-correctness or performance effects,
   and states why. It never lowers the floor. *(Our call.)*
3. **Draft the reading guide** from the diff structure:
   the order to read files or commits in, one line per step on what it does,
   and a list of spots needing close attention.
4. **Dispatch reviewers in parallel.**
   Small PRs: one reviewer that covers all three prompts.
   Larger PRs: the three focused reviewers.
   Each gets the context from step 1, the guidelines, and its prompt. Never the coordinator's history.
   Reviewers get the diff and the file list and read surrounding code on demand.
   Very large PRs (above a second threshold) are additionally split by subsystem,
   each subsystem group getting its own reviewer set. *(Our call.)*
5. **Verify findings.**
   Each finding is checked against the code; rejected ones are dropped.
   Findings a reviewer set aside as out of scope are listed, not silently dropped.
   Findings that static review can neither confirm nor reject
   (MPI races, numerical drift, performance claims) move to the reading guide's
   attention list, marked "needs human check". They do not count as findings
   and do not affect the verdict. *(Our call.)*
6. **Merge** reviewer notes on hotspots into the reading guide's attention list.
7. **Write the report** to `.claude/reviews/…` and print its path plus a three-line summary.
8. **Post only on explicit request**, after showing the user what will be posted.

### Risk levels (risk-levels.md)

The file has two parts: a team-editable path-to-level table (the floor),
and the level definitions below, which the coordinator uses to raise the rating.

* **Low:** docs, tests, samples, or code confined to one opt-in feature.
  Few users affected; mistakes are easy to spot.
* **Medium:** one physics method, interaction, or script-interface class;
  or user-visible changes (API, defaults, numerical results) within one method.
* **High:** shared infrastructure every simulation passes through
  (integrator, cell system, particle data, ghost/MPI communication,
  script-interface base classes, build system, CI);
  or anything that silently changes results, parallel correctness, or performance.

Each rating carries one sentence: what could break and who would notice.

### Guidelines file (espresso-guidelines.md)

Two sections:

* **Goals:** the three goals, each with a few lines on what it means in practice.
* **Checks:** short bullets grouped by area:
  C++ core, script interface, Python API, tests, docs, build.

Seed content comes from `AGENT.md`, `CONTRIBUTING.md` and the wiki page "Code review".
A header explains how to add a check (one bullet, one PR).

### Report layout (report-template.md)

1. Header: PR number or branch, title, base and head SHA, date.
2. **Risk:** level, reason sentence, files that set it.
3. **What the PR does:** three to five sentences.
4. **Reading guide:** ordered steps; attention list.
5. **Findings:** Critical / Important / Minor.
   Each with `file:line`, what is wrong, why it matters, the goal or check it violates, a fix if not obvious.
6. **Strengths:** short, specific.
7. **Set aside:** behaviors considered and left out of scope, one line each.
8. **Verdict:** ready / ready with fixes / not ready, with one or two sentences.

## Validation

A smoke test, not a full evaluation suite. *(Our call.)*
The builder picks two or three merged ESPResSo PRs
whose review threads found real problems, plus one clean PR.
The skill should report the known problems and stay quiet on the clean one.
Findings that miss or add noise feed back into the prompts.
A curated evaluation set is a possible follow-up, out of scope here.

## Left to the builder

* The line-count threshold between "small" and "larger" PRs.
* The exact `git`/`gh` commands and how the upstream remote is detected.
* Verification batching (one verifier per reviewer, or per finding).
* Wording of prompts, and the seed checks drawn from the sources above.
* Format of the GitHub post when the user asks for one.
