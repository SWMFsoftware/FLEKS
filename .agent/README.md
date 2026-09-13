# .agent/ — Agent Skills, Workflows and Knowledge Base

This directory holds the machine-facing documentation for FLEKS. It is
versioned with the code: **the repository is the single source of truth**, and
the local CodeBuddy installation is derived from it with
`tools/install_skill.sh`.

## Layout

| Path | Contents |
|---|---|
| `agents/fleks.md` | The dedicated implementation agent: context-heavy exploration and long build/test loops, reporting concise summaries. |
| `skills/fleks-expert/` | Knowledge base: `SKILL.md` router + `references/` loaded on demand (architecture, file layout, standards, testing, coupling). |
| `skills/build-fleks/`, `add-new-source/`, `code-cleanup/`, `debug-session/`, `generate-docs/` | Task recipes — one skill per recurring activity. |
| `workflows/` | Multi-step guides meant to be invoked explicitly (`add-param`, `add-coupling-var`, `run-test`). |

## Install into CodeBuddy

```bash
tools/install_skill.sh                 # symlink into ~/.codebuddy/skills
tools/install_skill.sh --copy          # copy instead (if symlinks are not followed)
tools/install_skill.sh --force         # replace already installed skills
tools/install_skill.sh --dest /path/to/skills
```

Symlinks are the default so the installed copy cannot drift from this one.

## Authoring rules

1. **Task skills stay task-shaped.** A skill describes *how to do one thing*;
   shared facts (architecture, conventions, test catalogue) belong in
   `skills/fleks-expert/references/` and are referenced, not repeated.
2. **Every `SKILL.md` needs frontmatter** with `name` and a `description` that
   says *when* to use the skill.
3. **Keep `SKILL.md` short** (≤ ~120 lines). Anything not needed on every
   invocation goes into `references/`.
4. **Human-facing docs live in `doc/` and `README.md`**, not here. If something
   is useful to both audiences, write it once in `doc/` and link to it.
5. **No per-directory `AGENT.md` files.** The root `AGENT.md` is a router;
   detailed content belongs in `doc/` or in a reference here.

Run `python3 tools/check_docs.py` (also in CI) after doc changes: it enforces
the single-root `AGENT.md` rule, skill/workflow frontmatter, that every
reference is linked from `SKILL.md`, and that backticked repository paths
resolve.
