# ESPResSo PR Review Skill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A team-shared Claude Code skill that reviews ESPResSo PRs and writes a local Markdown report with risk rating, reading guide, and verified findings.

**Architecture:** A coordinator procedure (`SKILL.md`) collects PR context, rates risk from a path table, drafts a reading guide, dispatches focused reviewer subagents, verifies their findings, and fills a report template. All project rules live in one team-editable file; prompts reference it and contain no rules.

**Tech Stack:** Claude Code skills (Markdown + YAML frontmatter), `git`, `gh` CLI.

**Spec:** `docs/superpowers/specs/2026-11-20-espresso-pr-review-skill-design.md`

## Global Constraints

* Skill directory: `.claude/skills/espresso-pr-review/`.
* Report path: `.claude/reviews/<pr-N or branch-name>.md`; `.claude/reviews/` is git-ignored.
* Never post to GitHub unless the user asks in words (e.g. "post this to the PR"); show the text before posting.
* Never modify the user's checkout: no checkout, stash, reset, or branch switch. PR heads are fetched to `refs/review/pr-<N>`.
* Base defaults to `git merge-base <upstream>/python <head>`. Upstream = the remote whose URL contains `espressomd/espresso` (here `origin`).
* `espresso-guidelines.md` is the only rules source at runtime. No wiki fetch at runtime.
* No dependency on the superpowers plugin.
* Prompts contain no project rules; they tell reviewers to read `espresso-guidelines.md`.
* Do not report issues that `clang-format`, `autopep8`, `pylint`, or `clang-tidy` catch.
* Report section order: header, Risk, What the PR does, Reading guide, Findings, Strengths, Set aside, Verdict.
* Risk levels: Low, Medium, High. Path table sets a floor; coordinator may raise after verification, never lower.
* Unconfirmable findings go to the reading guide attention list as "needs human check"; they are not findings and do not affect the verdict.
* Comments and docs describe the current state, no change history.

## Review Focus

1. **PR from a fork whose head ref is gone or force-pushed:** the skill fetches `pull/<N>/head` and reports the SHA it reviewed in the header. Pinned in Task 4, Step 2.
2. **User runs the skill in a checkout with uncommitted changes:** review still uses refs only and leaves the working tree untouched. Pinned in Task 4, Step 3 (`git status` before/after identical).
3. **Branch name containing `/` (e.g. `me/comm`):** report filename replaces `/` with `_`. Pinned in Task 4, Step 3.
4. **Diff touching only docs/tests:** rating is Low and the single-reviewer path runs. Pinned in Task 5 (clean PR case).
5. **User says "post it":** skill shows the text and asks once before `gh` posts; it never posts on "looks good" alone. Pinned in Task 4, Step 4.

---

## File Structure

```
.gitignore                                      # add .claude/reviews/
.claude/skills/espresso-pr-review/
  SKILL.md                  # coordinator procedure (Task 4)
  espresso-guidelines.md    # goals + checks, team-editable (Task 1)
  risk-levels.md            # path table + level definitions (Task 2)
  report-template.md        # output layout (Task 3)
  prompts/
    correctness.md          # reviewer (Task 3)
    espresso-goals.md       # reviewer (Task 3)
    tests-docs.md           # reviewer (Task 3)
    verify.md               # verifier (Task 3)
```

## Shared Formats (used by Tasks 3 and 4)

Reviewer output, one block per finding:

```
### F<n>: <short title>
- severity: Critical | Important | Minor
- location: <path>:<line>[-<line>]
- rule: <goal or check heading from espresso-guidelines.md, or "defect">
- problem: <what is wrong>
- impact: <why it matters>
- fix: <suggestion, or "obvious">
```

Followed by `## Hotspots` (bullets: `<path>:<line> — <why a human should look>`) and `## Set aside` (bullets: behavior considered, reason).

Verifier output per finding: `F<n>: confirmed | rejected | needs-human-check — <one-line reason with path:line evidence>`.

---

### Task 1: Guidelines file and gitignore

**Files:**
- Create: `.claude/skills/espresso-pr-review/espresso-guidelines.md`
- Modify: `.gitignore` (append `.claude/reviews/`)

**Interfaces:**
- Produces: headings `## Goals` with subsections `### Understandability`, `### Maintainability`, `### Use of existing infrastructure`; heading `## Checks` with subsections `### C++ core`, `### Script interface`, `### Python API`, `### Tests`, `### Docs`, `### Build and CI`. Reviewers cite these headings in the `rule:` field.

- [ ] **Step 1: Gather sources.** Read `AGENT.md`, `CONTRIBUTING.md`, `.github/PULL_REQUEST_TEMPLATE.md`. Fetch https://github.com/espressomd/espresso/wiki/Code-review once (WebFetch or web-fetch agent).
- [ ] **Step 2: Write the file.** Top: a 3-line "How to extend" note (one bullet per check, one PR; cite a source link where one exists). Goals: 3–6 lines each on what the goal means in practice for ESPResSo. Checks: short bullets, each phrased as something a reviewer can verify in a diff. Link to `AGENT.md` for infrastructure details instead of copying. Exclude anything formatters/linters catch.
- [ ] **Step 3: Verify.**
  Run: `grep -c '^### ' .claude/skills/espresso-pr-review/espresso-guidelines.md` → Expected: `9`.
  Run: `git check-ignore .claude/reviews/x.md` → Expected: prints the path.
- [ ] **Step 4: Commit** (`git add` both files; message "Add ESPResSo review guidelines for PR review skill").

### Task 2: Risk levels

**Files:**
- Create: `.claude/skills/espresso-pr-review/risk-levels.md`

**Interfaces:**
- Produces: `## Path table` (Markdown table, columns `Pattern | Level`, glob patterns relative to repo root, first match wins, ordered High → Medium → Low, final row `* | Medium`); `## Level definitions` with `### Low`, `### Medium`, `### High` using the spec's wording verbatim; `## Raising the rating` listing the triggers (silent result change, parallel correctness, performance) and the rule "never lower the floor".

- [ ] **Step 1: Write the file.** High rows at minimum: `src/core/integrate*`, `src/core/integrators/**`, `src/core/cell_system/**`, `src/core/cells.*`, `src/core/Particle*.hpp`, `src/core/ghosts*`, `src/core/ghosts/**`, `src/core/communication.*`, `src/core/MpiCallbacks.hpp`, `src/core/forces*`, `src/core/system/**`, `src/script_interface/*.hpp`, `src/script_interface/auto_parameters/**`, `src/utils/**`, `CMakeLists.txt`, `cmake/**`, `.gitlab-ci.yml`, `.github/**`, `maintainer/CI/**`. Medium: `src/**`. Low: `doc/**`, `testsuite/**`, `samples/**`, `maintainer/**`, `*.md`. Check each High pattern matches at least one file.
- [ ] **Step 2: Verify.**
  Run: for each High pattern, `git ls-files '<pattern>' | head -1` → Expected: non-empty for every row.
- [ ] **Step 3: Commit.**

### Task 3: Report template and prompts

**Files:**
- Create: `.claude/skills/espresso-pr-review/report-template.md`
- Create: `.claude/skills/espresso-pr-review/prompts/{correctness,espresso-goals,tests-docs,verify}.md`

**Interfaces:**
- Consumes: guideline headings (Task 1), Shared Formats above.
- Produces: placeholders the coordinator fills, identical in all prompts: `{{CONTEXT}}` (PR title, description, linked issues, base/head SHA, commit list, changed-file stats), `{{DIFF_CMD}}` (the `git diff` command to run), `{{FILES}}` (file subset for this reviewer; whole list unless split), `{{GUIDELINES_PATH}}`. Verifier also gets `{{FINDINGS}}`. Template placeholders: `{{TITLE}}`, `{{REF}}`, `{{BASE}}`, `{{HEAD}}`, `{{DATE}}`, plus one per section.

- [ ] **Step 1: Write `report-template.md`.** Sections in the Global Constraints order. Reading guide has `### Read in this order` (numbered: path or commit, one line on what it does) and `### Pay attention to` (bullets; "needs human check" items marked so). Findings grouped by severity, each in the Shared Format minus the `F<n>` id plus the `rule`.
- [ ] **Step 2: Write the three reviewer prompts.** Each: role, scope (its topic only), "read `{{GUIDELINES_PATH}}` first", read-only rule, "read surrounding code on demand", no subagents, skip linter-catchable issues, output exactly the Shared Format. Topics:
  * `correctness.md`: logic bugs, MPI/parallel correctness (rank-local vs global state, collective call symmetry, ghost updates), numerical correctness, error handling, resource lifetime.
  * `espresso-goals.md`: the three Goals sections plus the area Checks.
  * `tests-docs.md`: tests exercise new behavior incl. parallel runs where relevant, user docs in `doc/sphinx` for user-visible changes, samples/tutorials still valid, CHANGELOG-worthy changes flagged.
- [ ] **Step 3: Write `verify.md`.** Instructions: for each finding, open the code at `location`, try to confirm or refute; output the Verifier format. "needs-human-check" only when static reading cannot settle it (races, drift, performance); otherwise decide.
- [ ] **Step 4: Verify.**
  Run: `grep -L 'GUIDELINES_PATH' .claude/skills/espresso-pr-review/prompts/*.md` → Expected: only `verify.md` (or empty if it references guidelines too).
  Run: `grep -niE 'clang-format|autopep8' .claude/skills/espresso-pr-review/prompts/correctness.md` → Expected: the skip-linter line is present.
- [ ] **Step 5: Commit.**

### Task 4: Coordinator SKILL.md

**Files:**
- Create: `.claude/skills/espresso-pr-review/SKILL.md`

**Interfaces:**
- Consumes: everything from Tasks 1–3.
- Produces: frontmatter `name: espresso-pr-review`; `description:` starting "Use when reviewing an ESPResSo pull request or branch…". Invocation `/espresso-pr-review <PR number | branch | base..head>`.

- [ ] **Step 1: Write SKILL.md** with the spec's 8 coordinator steps. Decisions to encode:
  * Input resolution: integer → PR: `gh pr view <N> --repo espressomd/espresso --json title,body,headRefOid,baseRefName,files,commits,closingIssuesReferences`, then `git fetch <upstream> pull/<N>/head:refs/review/pr-<N>`. Branch → local ref. `a..b` → used as given.
  * Thresholds: single reviewer if ≤ 400 changed lines; three reviewers above; additionally split by top-level subsystem (`src/core/<dir>`, `src/script_interface`, `src/python`, `src/walberla_bridge`, `testsuite`, `doc`, other) if > 3000 lines or > 40 files.
  * Verification: one verifier subagent per reviewer output.
  * Report filename: `pr-<N>.md`, or branch name with `/` → `_`.
  * Final message: report path, risk level, finding counts by severity, verdict.
  * Posting: only on an explicit request; show the text; then `gh pr comment <N> --repo espressomd/espresso --body-file <report>`.
- [ ] **Step 2: Dry run on a PR.** In a fresh subagent: "Follow `.claude/skills/espresso-pr-review/SKILL.md` for PR 5470." Expected: `.claude/reviews/pr-5470.md` exists; header head SHA equals `gh pr view 5470 --json headRefOid`; all 8 sections present in order.
- [ ] **Step 3: Dry run on a branch, with a dirty tree.** `echo x >> AUTHORS`; capture `git status --porcelain` and `git rev-parse HEAD`; run the skill on `me/comm` (or any slash-named local branch). Expected: `.claude/reviews/me_comm.md` exists; status and HEAD unchanged; then `git checkout AUTHORS`.
- [ ] **Step 4: Check posting rule.** Grep SKILL.md for the explicit-request rule and the show-before-post step. Expected: both present; no `gh pr comment`/`gh pr review` appears outside the posting section.
- [ ] **Step 5: Commit.**

### Task 5: Smoke test on past PRs

**Files:**
- Modify (if needed): prompts and guidelines from Tasks 1–3.

- [ ] **Step 1: Pick PRs.** Use `gh pr list --repo espressomd/espresso --state merged --limit 100` and `gh pr view <N> --comments` to find 2–3 merged PRs where review comments led to fix commits, plus 1 docs/tests-only PR without problems. Record numbers and the known problems in the commit message of Step 4.
- [ ] **Step 2: Run the skill on each PR at the head SHA reviewed before the fix** (use `gh api repos/espressomd/espresso/pulls/<N>/commits` to find it; pass `base..sha`).
- [ ] **Step 3: Compare.** Expected: each known problem appears as a finding or a "needs human check" item; the clean PR has Risk Low and no Critical/Important findings. If not, adjust prompts or guidelines and rerun that PR (max two rounds; report remaining misses to the user).
- [ ] **Step 4: Commit** prompt changes with a message listing PRs used and the outcome.
