#!/usr/bin/env python3
"""Check the consistency of the FLEKS agent/documentation tree.

Run with no arguments from anywhere inside the repository:

    python3 tools/check_docs.py

The checks encode the rules agreed for the documentation layout:

1. There is exactly one ``AGENT.md``, at the repository root, and it stays a
   short router (not a manual).
2. Every agent skill has a ``SKILL.md`` with ``name`` and ``description``
   frontmatter; every workflow has a ``description``.
3. Every reference file under ``skills/fleks-expert/references/`` is reachable
   from its ``SKILL.md``.
4. Backticked repository paths in the agent docs resolve to something real
   (catches the ``Doc/`` vs ``doc/`` class of drift).
5. Documentation paths use the lowercase ``doc/`` directory.
6. The canonical documents exist.

Exits with status 1 and a list of problems when a check fails.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

MAX_AGENT_MD_LINES = 150
DOC_EXTENSIONS = {".md", ".tex", ".py", ".sh", ".xml", ".yml", ".yaml", ".in"}
PATH_EXTENSIONS = DOC_EXTENSIONS | {
    ".cpp",
    ".h",
    ".hpp",
    ".f90",
    ".F90",
    ".pdf",
    ".txt",
    ".json",
}
SCAN_DIRS = ["doc", ".agent"]
SCAN_FILES = ["AGENT.md", "README.md", "CONTRIBUTING.md", "PARAM.XML"]
SCAN_GLOBS = ["tests/*/README.md"]

CANONICAL_FILES = [
    "doc/DEVELOPING.md",
    "doc/Coding_standards.md",
    "doc/Algorithm.tex",
    "tests/README.md",
    "PARAM.XML",
]

# Backticked strings that are not repository paths.
PATH_SKIP_CHARS = set("#$*()=<>|& \t{}")
PATH_RE = re.compile(r"`([^`\n]+)`")

# Paths that only exist after a build, a run, or in a parent SWMF tree.
GENERATED_PATHS = {
    "doc/USERMANUAL.pdf",
    "doc/Algorithm.pdf",
    "bin/FLEKS.exe",
    "bin/converter.exe",
    "./FLEKS.exe",
    "src/libFLEKS.a",
    "compile_commands.json",
    "include/show_git_info.h",
    "html/index.html",
    "SWMF/PC/FLEKS",
    "../../util/AMREX/InstallDir/",
    "-I../include",
    "FLEKS1/",
    "output/test16/",
    "run_test/runlog",
    "run_test/RESULTS/3d/PC/",
    "tests/summary.md",
    "tests/performance_summary.md",
    "GM/BATSRUS/srcInterface/ModGridDescriptor.f90",
}

# Files that live in the parent SWMF/BATSRUS tree, and filename patterns used
# in recipes (they are real names but not paths inside this repository).
EXTERNAL_PATHS = {
    "Source.h",
    "GM_couple_pc.f90",
    "GM_wrapper.f90",
    "ModUserMars.f90",
    "PC_wrapper.f90",
}

# Placeholders used in "how to add a file" recipes.
PLACEHOLDERS = ("NewFeature", "NewSource", "MyIC", "New")
UPPER_DOC_RE = re.compile(r"(?<![\w-])Doc/")
FRONTMATTER_RE = re.compile(r"\A---\n(.*?)\n---\n", re.DOTALL)


def markdown_files() -> list[Path]:
    files: list[Path] = []
    for name in SCAN_FILES:
        path = REPO_ROOT / name
        if path.is_file():
            files.append(path)
    for directory in SCAN_DIRS:
        base = REPO_ROOT / directory
        if base.is_dir():
            files.extend(sorted(base.rglob("*.md")))
    for pattern in SCAN_GLOBS:
        files.extend(sorted(REPO_ROOT.glob(pattern)))
    return sorted(set(files))


def frontmatter(text: str) -> dict[str, str]:
    match = FRONTMATTER_RE.match(text)
    if not match:
        return {}
    fields: dict[str, str] = {}
    for line in match.group(1).splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip()
    return fields


def resolve(candidate: str, source: Path) -> bool:
    """True when a backticked candidate points at something that exists.

    A path may be written relative to the file that mentions it
    (``references/standards.md``), relative to the repository root
    (``.agent/skills/...``) or as a bare name resolved by search
    (``validate.py``).
    """
    path = candidate.rstrip("/")
    allowed = GENERATED_PATHS | EXTERNAL_PATHS
    if not path or path in allowed or candidate in allowed:
        return True
    for base in (source.parent, REPO_ROOT):
        if (base / path).exists():
            return True
    name = Path(path).name
    if any(REPO_ROOT.glob(f"**/{name}")):
        return True
    return False


def looks_like_path(candidate: str) -> bool:
    if any(char in PATH_SKIP_CHARS for char in candidate):
        return False
    if "://" in candidate or candidate.startswith(("-", "http")):
        return False
    if any(holder in candidate for holder in PLACEHOLDERS):
        return False
    if candidate.endswith("/"):
        return True
    # Anything else must name a file with a recognized extension, so that
    # inline math such as `5/3` or `dAy/dt` is not mistaken for a path.
    return Path(candidate).suffix in PATH_EXTENSIONS


def check_agent_md(errors: list[str]) -> None:
    root = REPO_ROOT / "AGENT.md"
    if not root.is_file():
        errors.append("AGENT.md is missing from the repository root")
        return
    lines = root.read_text(encoding="utf-8").splitlines()
    if len(lines) > MAX_AGENT_MD_LINES:
        errors.append(
            f"AGENT.md is {len(lines)} lines; it must stay a router "
            f"(<= {MAX_AGENT_MD_LINES}). Move detail into doc/DEVELOPING.md."
        )
    for stray in sorted(REPO_ROOT.rglob("AGENT.md")):
        if stray != root:
            errors.append(
                f"unexpected AGENT.md at {stray.relative_to(REPO_ROOT)}; "
                "only the root router is allowed"
            )


def check_skills(errors: list[str]) -> None:
    skills_dir = REPO_ROOT / ".agent" / "skills"
    if not skills_dir.is_dir():
        errors.append(".agent/skills/ is missing")
        return
    for skill in sorted(skills_dir.iterdir()):
        if not skill.is_dir():
            continue
        skill_md = skill / "SKILL.md"
        if not skill_md.is_file():
            errors.append(f"{skill.name}: missing SKILL.md")
            continue
        fields = frontmatter(skill_md.read_text(encoding="utf-8"))
        for key in ("name", "description"):
            if not fields.get(key):
                errors.append(f"{skill.name}/SKILL.md: frontmatter lacks '{key}'")

    workflows_dir = REPO_ROOT / ".agent" / "workflows"
    for workflow in sorted(workflows_dir.glob("*.md")) if workflows_dir.is_dir() else []:
        if not frontmatter(workflow.read_text(encoding="utf-8")).get("description"):
            errors.append(f"workflows/{workflow.name}: frontmatter lacks 'description'")


def check_references(errors: list[str]) -> None:
    skill_md = REPO_ROOT / ".agent" / "skills" / "fleks-expert" / "SKILL.md"
    refs_dir = skill_md.parent / "references"
    if not skill_md.is_file() or not refs_dir.is_dir():
        return
    skill_text = skill_md.read_text(encoding="utf-8")
    for reference in sorted(refs_dir.glob("*.md")):
        if f"references/{reference.name}" not in skill_text:
            errors.append(
                f"references/{reference.name} is not linked from "
                "skills/fleks-expert/SKILL.md"
            )


def check_paths(errors: list[str]) -> None:
    for path in markdown_files():
        text = path.read_text(encoding="utf-8")
        rel = path.relative_to(REPO_ROOT)
        for match in UPPER_DOC_RE.finditer(text):
            line = text[: match.start()].count("\n") + 1
            errors.append(f"{rel}:{line}: use the lowercase 'doc/' directory")
        for candidate in PATH_RE.findall(text):
            if not looks_like_path(candidate):
                continue
            if not resolve(candidate, path):
                line = text[: text.find(f"`{candidate}`")].count("\n") + 1
                errors.append(f"{rel}:{line}: broken path reference `{candidate}`")


def check_canonical(errors: list[str]) -> None:
    for name in CANONICAL_FILES:
        if not (REPO_ROOT / name).exists():
            errors.append(f"canonical document missing: {name}")


def main() -> int:
    errors: list[str] = []
    check_agent_md(errors)
    check_skills(errors)
    check_references(errors)
    check_paths(errors)
    check_canonical(errors)

    if not errors:
        print("docs check: OK")
        return 0

    print(f"docs check: {len(errors)} problem(s)")
    for error in errors:
        print(f"  - {error}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
