# .agents/ — Agent Skills, Workflows and Knowledge Base

This directory holds the machine-facing documentation for FLEKS. It is
versioned with the code: **the repository is the single source of truth**, and
every tool either reads it directly or links to it (see *Using these skills*).

## Layout

| Path | Contents |
|---|---|
| `agents/fleks.md` | The dedicated implementation agent: context-heavy exploration and long build/test loops, reporting concise summaries. |
| `skills/fleks-expert/` | Knowledge base: `SKILL.md` router + `references/` loaded on demand (architecture, file layout, standards, testing, coupling). |
| `skills/build-fleks/`, `add-new-source/`, `code-cleanup/`, `debug-session/`, `generate-docs/` | Task recipes — one skill per recurring activity. |
| `workflows/` | Multi-step guides meant to be invoked explicitly (`add-param`, `add-coupling-var`, `run-test`). |

## Using these skills

`SKILL.md` (YAML frontmatter + progressive disclosure) is an open format, but
each tool looks for skills in its own place:

| Tool | Discovery | What to do |
|---|---|---|
| Antigravity | `.agents/skills/` in the workspace | nothing — this is the native layout |
| Claude Code | `.claude/skills/` (project) or `~/.claude/skills/` | link the skills you want |
| CodeBuddy | `~/.codebuddy/skills/` | link the skills you want |
| Codex | project instructions file (AGENTS.md by convention) | add a pointer to `AGENT.md` and `.agents/` |

Linking is a single command per skill; symlinks keep the installed copy from
drifting out of the repository:

```bash
# Claude Code (project-level)
mkdir -p .claude/skills
for s in .agents/skills/*/; do ln -sfn "$PWD/$s" ".claude/skills/$(basename "$s")"; done

# CodeBuddy (user-level)
for s in .agents/skills/*/; do ln -sfn "$PWD/$s" "$HOME/.codebuddy/skills/$(basename "$s")"; done
```

Use copies instead if a tool does not follow symlinks.

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
