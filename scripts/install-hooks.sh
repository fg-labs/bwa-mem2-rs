#!/usr/bin/env bash
# Install git pre-commit hook for bwa-mem3-rs.
# Usage: ./scripts/install-hooks.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
HOOKS_SOURCE="${SCRIPT_DIR}/hooks"
HOOKS_DEST="${REPO_ROOT}/.git/hooks"

if [[ ! -f "${REPO_ROOT}/Cargo.toml" ]]; then
    echo "error: must be run from the bwa-mem3-rs repository" >&2
    exit 1
fi
if [[ ! -d "${REPO_ROOT}/.git" ]]; then
    echo "error: .git directory not found" >&2
    exit 1
fi

mkdir -p "${HOOKS_DEST}"

for hook in "${HOOKS_SOURCE}"/*; do
    [[ -f "${hook}" ]] || continue
    hook_name="$(basename "${hook}")"
    dest="${HOOKS_DEST}/${hook_name}"
    if [[ -e "${dest}" ]]; then
        if [[ -L "${dest}" ]]; then
            rm "${dest}"
        else
            mv "${dest}" "${dest}.backup"
            echo "backed up existing ${hook_name} to ${hook_name}.backup"
        fi
    fi
    ln -s "${hook}" "${dest}"
    chmod +x "${hook}"
    echo "installed: ${hook_name}"
done

echo ""
echo "hooks installed. bypass with 'git commit --no-verify' when needed."
