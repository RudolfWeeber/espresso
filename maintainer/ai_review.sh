#!/usr/bin/env bash
#
# Copyright (C) 2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Run the AI-based PR review (Claude Code skill `espresso-pr-review`)
# non-interactively and print the risk level, the AI review recommendation
# and the path of the report. The review is never posted to GitHub.

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") pr <number>
       $(basename "$0") branch <name>

Runs the AI-based review of a GitHub pull request or a local branch with
Claude Code and writes the report to .claude/reviews/ in the repository.
Requires the 'claude' CLI; PR reviews also require an authenticated 'gh'.
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
  usage
  exit 0
fi
if [ "$#" -ne 2 ]; then
  usage >&2
  exit 2
fi

mode="$1"
target="$2"
repo_root="$(git rev-parse --show-toplevel)"

case "${mode}" in
  pr)
    if ! [[ "${target}" =~ ^[0-9]+$ ]]; then
      echo "error: PR number must be an integer, got '${target}'" >&2
      exit 2
    fi
    report="${repo_root}/.claude/reviews/pr-${target}.md"
    ;;
  branch)
    if ! git -C "${repo_root}" rev-parse --verify --quiet "${target}^{commit}" >/dev/null; then
      echo "error: branch '${target}' not found" >&2
      exit 2
    fi
    report="${repo_root}/.claude/reviews/${target//\//_}.md"
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac

if ! command -v claude >/dev/null; then
  echo "error: the 'claude' CLI is not installed" >&2
  exit 1
fi

work_dir="$(mktemp -d "${TMPDIR:-/tmp}/ai_review.XXXXXX")"
log="${work_dir}/claude.log"
output="${work_dir}/report.md"

echo "Running AI-based review of ${mode} ${target} (this takes several minutes)..." >&2

# Read-only git and gh access; posting, pushing and checkout changes are denied.
if ! (cd "${repo_root}" && claude -p "/espresso-pr-review ${target} --output ${output}" \
    --permission-mode acceptEdits --add-dir "${work_dir}" \
    --allowedTools "Read" "Write" "Edit" "Glob" "Grep" "Agent" "Skill" \
      "Bash(git:*)" "Bash(gh pr view:*)" "Bash(gh issue view:*)" \
      "Bash(mkdir:*)" "Bash(cp:*)" "Bash(python3:*)" "Bash(ls:*)" \
    --disallowedTools "Bash(gh pr comment:*)" "Bash(gh pr review:*)" \
      "Bash(gh api:*)" "Bash(git push:*)" "Bash(git checkout:*)" \
      "Bash(git switch:*)" "Bash(git reset:*)" "Bash(git stash:*)" \
      "Bash(git commit:*)" \
    >"${log}" 2>&1); then
  echo "error: claude failed; last lines of ${log}:" >&2
  tail -n 5 "${log}" >&2
  exit 1
fi

if [ ! -s "${output}" ]; then
  echo "error: no report was written; last lines of ${log}:" >&2
  tail -n 5 "${log}" >&2
  exit 1
fi
mkdir -p "$(dirname "${report}")"
cp "${output}" "${report}"

risk="$(sed -n 's/^## Risk: *//p' "${report}" | head -n 1)"
recommendation="$(sed -n 's/^\* \*\*AI review recommendation:\*\* *//p' "${report}" | head -n 1)"

echo "Risk:           ${risk:-unknown}"
echo "Recommendation: ${recommendation:-unknown}"
echo "Report:         ${report}"
rm -rf "${work_dir}"
