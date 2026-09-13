#!/usr/bin/env bash
#
# Install the FLEKS agent skills from this repository into the CodeBuddy user
# skills directory.
#
# By default the skills are installed as SYMLINKS back into the repository, so
# the installed copy can never drift from the versioned one. Use --copy if the
# runtime does not follow symlinks.
#
# Usage:
#   tools/install_skill.sh [options]
#
# Options:
#   --dest DIR   Install into DIR (default: $CODEBUDDY_SKILLS_DIR or
#                $HOME/.codebuddy/skills)
#   --copy       Copy files instead of symlinking
#   --force      Replace an already installed skill
#   -h, --help   Show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC_DIR="${REPO_ROOT}/.agent/skills"
DEST_DIR="${CODEBUDDY_SKILLS_DIR:-${HOME}/.codebuddy/skills}"

MODE="symlink"
FORCE=0

usage() { sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
  case "$1" in
    --dest)  DEST_DIR="$2"; shift 2 ;;
    --copy)  MODE="copy"; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [ ! -d "${SRC_DIR}" ]; then
  echo "No skills directory found at ${SRC_DIR}" >&2
  exit 1
fi

mkdir -p "${DEST_DIR}"

installed=0
skipped=0

for src in "${SRC_DIR}"/*; do
  [ -d "${src}" ] || continue
  name="$(basename "${src}")"
  if [ ! -f "${src}/SKILL.md" ]; then
    echo "skip   ${name} (no SKILL.md)"
    skipped=$((skipped + 1))
    continue
  fi

  dest="${DEST_DIR}/${name}"
  if [ -e "${dest}" ] || [ -L "${dest}" ]; then
    if [ "${FORCE}" -eq 0 ]; then
      echo "skip   ${name} (already installed; use --force to replace)"
      skipped=$((skipped + 1))
      continue
    fi
    rm -rf "${dest}"
  fi

  if [ "${MODE}" = "copy" ]; then
    cp -R "${src}" "${dest}"
    echo "copied ${name} -> ${dest}"
  else
    ln -s "${src}" "${dest}"
    echo "linked ${name} -> ${dest}"
  fi
  installed=$((installed + 1))
done

echo ""
echo "Installed ${installed} skill(s) into ${DEST_DIR} (${MODE} mode), skipped ${skipped}."
