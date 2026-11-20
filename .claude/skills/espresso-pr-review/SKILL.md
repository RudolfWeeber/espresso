---
name: espresso-pr-review
description: Use when reviewing an ESPResSo pull request or branch, e.g. "/espresso-pr-review 5470" or "review branch my_feature". Writes a local Markdown report with a risk rating, a reading guide for human reviewers, and verified findings against the ESPResSo guidelines. Never posts to GitHub unless explicitly asked.
---

# ESPResSo PR review

You coordinate a review of one ESPResSo pull request or branch.
The result is a local report for the PR author and for human reviewers:
a risk rating, a reading guide, and findings that were checked against the code.

**Usage:** `/espresso-pr-review <PR number | branch | base..head> [--output <path>]`

`--output` writes the report to the given path instead of `.claude/reviews/`
(used by `maintainer/ai_review.sh`).

## Ground rules

* **Never post to GitHub** unless the user asks for it in words, e.g. "post this to the PR".
  Approval of the report ("looks good") is not a request to post.
* **Never modify the checkout.** No `checkout`, `switch`, `stash`, `reset`, `commit`, or edits to tracked files.
  Fetch into `refs/review/*` and read with `git show`, `git diff`, `git log`.
* **Project rules live in `espresso-guidelines.md`.** Do not add rules of your own.
* **Do not write findings yourself.** Read the diff only to draft the reading guide.
  Findings come from reviewers and pass through verification.
* **Run in the main session**, not as a subagent: the skill dispatches its own subagents and writes the report file.

All paths below are relative to this skill's directory, `<repo>/.claude/skills/espresso-pr-review/`,
where `<repo>` is `git rev-parse --show-toplevel`.

The commands below use shell variables (`BASE`, `HEAD`) for readability.
If your environment refuses command substitution, run each command alone and use the printed values.
Keep temporary files (filled prompts, context) in a folder unique to this run,
e.g. `${TMPDIR:-/tmp}/espresso-pr-review-<name>-<head7>/`.

A branch or range review covers committed code only; uncommitted changes in the working tree are not reviewed.
Fetched refs stay under `refs/review/`; remove them with
`git for-each-ref refs/review --format='delete %(refname)' | git update-ref --stdin`.

## Step 1: Collect context

**Upstream remote:** the first remote in `git remote -v` whose URL contains `espressomd/espresso`.
If there is none, use the URL `https://github.com/espressomd/espresso` in place of a remote name.

**Resolve the input:**

* **PR number `N`:**
  ```
  gh pr view N --repo espressomd/espresso \
    --json number,title,body,headRefOid,baseRefName,commits,files,closingIssuesReferences
  git fetch <upstream> pull/N/head:refs/review/pr-N <baseRefName>:refs/review/pr-N-base
  ```
  `HEAD_REF=refs/review/pr-N`. Check that `git rev-parse $HEAD_REF` equals `headRefOid`;
  if not, report the mismatch to the user and stop.
  `BASE=$(git merge-base refs/review/pr-N-base $HEAD_REF)`.
  For each linked issue, read its title and body with `gh issue view <n> --repo espressomd/espresso`.
* **Branch name:** `HEAD_REF=<branch>`.
  `git fetch <upstream> python:refs/review/python-base`, then
  `BASE=$(git merge-base refs/review/python-base <branch>)`.
  Title and description come from the commit messages.
* **Range `a..b`:** `BASE=a`, `HEAD_REF=b`, used as given.
  For a branch or range, read issues that the commit messages reference ("Fixes #n") with `gh issue view`.

Then collect:

* `HEAD=$(git rev-parse $HEAD_REF)`
* `git log --format='%h %s' $BASE..$HEAD`
* `git diff --stat $BASE $HEAD` and `git diff --shortstat $BASE $HEAD` (changed lines = insertions + deletions)

Write `{{CONTEXT}}` in this layout:

```
Title: <title>
Reference: PR #<N> | branch <name> | range <a..b>
Repository root: <repo>   (run all git commands here)
Base: <BASE sha>
Head: <HEAD sha>

Description:
<PR description, or the commit messages for a branch>

Linked issues:
- #<n> <title>: <first paragraph>

Commits:
<git log output>

Diff stat:
<git diff --stat output>
```
`{{DIFF_CMD}}` is `git diff <BASE> <HEAD>`.
`{{GUIDELINES_PATH}}` is the absolute path of `espresso-guidelines.md`.

## Step 2: Rate the risk floor

Apply the path table in `risk-levels.md` to every changed file
(`git diff --name-only $BASE $HEAD`).
`*` matches any characters within one path component, `**` matches any number of whole components,
and the first matching row wins.
The floor is the highest level reached. Remember the files that reached it.

## Step 3: Draft the reading guide

From the diff stat, the commit list, and a skim of the diff, write:

* **Read in this order:** a numbered list of code locations, each with line numbers and/or the function or class:
  `path:start-end` (`function_name()`). One line per item on what changed there and why it matters.
  Default order: new or changed interfaces (headers, base classes) → core logic →
  script interface → Python API → tests → docs.
  Group trivial or mechanical changes (renames, formatting, boilerplate) into one final item.
  When commits are self-contained steps, order by commit and give locations within each.
* **Pay attention to:** a first draft of the places where human judgement is most valuable.
  This is not a list of changed files. Include only what is unclear, controversial, risky or brittle:
  design decisions a maintainer should agree with, behavior changes, physics or numerics
  that need an expert eye, assumptions that may break later, parallel or performance effects.
  Phrase each item as the question or decision the human reviewer should settle, with its location.

## Step 4: Dispatch reviewers

Size the review by changed lines and changed files. The first matching row wins:

| Size | Reviewers |
|---|---|
| > 3000 lines or > 40 files | 4 reviewers per subsystem group |
| risk floor is Low | 1 reviewer with all four topics |
| > 400 lines | 4 reviewers, one per topic |
| otherwise | 1 reviewer with all four topics |

When the split would need more than 16 reviewers, tell the user the count and ask before dispatching.

**Subsystem groups** for the largest PRs: each `src/core/<dir>` (top-level files in `src/core` form one group),
`src/script_interface`, `src/python`, `src/walberla_bridge`, `testsuite`, `doc`, and everything else.
Merge groups with fewer than 100 changed lines into a neighbor so no reviewer gets a trivial share.

**Each reviewer prompt** is `prompts/reviewer-common.md` with its placeholders filled,
followed by the topic prompt(s): `prompts/correctness.md`, `prompts/physics.md`,
`prompts/espresso-goals.md`, `prompts/tests-docs.md`.
A subsystem group without physics code (e.g. `doc`, build files) skips the physics reviewer.
`{{FILES}}` is the reviewer's file list (all changed files unless split).
Send all reviewers in one message so they run in parallel, as `general-purpose` subagents.
Never pass your own conversation history.

## Step 5: Verify findings

For each reviewer output that contains findings, dispatch one verifier with `prompts/verify.md`,
`{{FINDINGS}}` set to that reviewer's finding blocks, verbatim. Run verifiers in parallel.

Apply the verdicts:

* `confirmed`: keep the reviewer's wording, with the verifier's corrections applied.
  If the verifier corrected the severity, file the finding under the new severity
  and append "(severity adjusted in verification)" to its title.
* `rejected`: drop.
* `needs-human-check`: move to the reading guide's attention list,
  prefixed "**Needs human check:**". These are not findings and do not affect the recommendation.

Merge duplicates across reviewers: same location and same problem become one finding with the higher severity.
Findings that overlap only partly stay separate.

## Step 6: Finalize risk and reading guide

* Raise the risk per `risk-levels.md` ("Raising the rating") if a confirmed finding calls for it.
  Never lower the floor.
  When the floor is already High, there is nothing to raise; name the finding in the risk sentence instead.
* Curate the attention list. Merge your draft with the reviewers' hotspots and the needs-human-check items.
  Keep only items where a human's judgement adds something the review could not settle;
  drop items that merely restate a change, a finding, or a file name.
  Aim for at most seven items, most important first.
  Prefix an item with "**Policy decision:**" only when the PR makes a deliberate choice with real trade-offs
  that the maintainers have not visibly agreed on (in the PR description or a linked issue), such as:
  * a change that breaks existing user scripts or data without a migration path;
  * a change of physics model or convention (e.g. how a potential, unit or sign is defined);
  * a changed default that alters simulation results;
  * a new required dependency, or dropped support for a feature, platform or toolchain.

  Not policy decisions: bug fixes that restore documented or evidently intended behavior,
  cosmetic changes to output or messages, internal refactoring, and test or doc changes.
  These go in the findings or the attention list without the prefix.
* Write one sentence on what could break and who would notice.
* Write the summary (what the PR does) from the context and the reading guide.
* Merge the reviewers' `## Strengths` lists. Drop duplicates and generic praise.
* Merge the reviewers' `## Set aside` lists. Keep only behaviors that were left unjudged;
  drop entries that say something was checked and found fine.

## Step 7: Write the report

Fill `report-template.md`, follow its filling rules, and remove the rules comment.
Write it to the `--output` path if one was given.
Otherwise write it to `<repo>/.claude/reviews/<name>.md`, creating the folder if needed:

* PR: `pr-<N>.md`
* branch: the branch name with `/` replaced by `_`, plus `.md`
* range: `<base7>_<head7>.md`

An existing report with the same name is overwritten; the header records the reviewed SHA.

**AI review recommendation:**
`Not ready` if any Critical finding;
`Ready with fixes` if any Important finding;
otherwise `Ready`.
Append "— policy decision required" when the attention list has at least one "**Policy decision:**" item,
and name the decisions in the recommendation's reasoning.

Tell the user, in this order:
the report path, the risk level, the number of findings per severity, and the AI review recommendation.

## Step 8: Post only on request

Only when the user asks for posting in words:

1. Write the text to post (the whole report, or the part the user named)
   to `post-body.md` in this run's temp folder.
2. Show the user exactly that text and ask once for confirmation.
3. Post that file: `gh pr comment <N> --repo espressomd/espresso --body-file <temp folder>/post-body.md`.

For a branch or range review there is no PR to post to; say so.
