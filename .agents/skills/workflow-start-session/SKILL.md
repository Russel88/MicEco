---
name: workflow-start-session
description: Start an engineering session with context and task loading Use when the user asks to start session in this repository or refers to the old `/StartSession` command.
---

# Start Session

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# StartSession

Start a new engineering session by loading project context and active tasks.

## Arguments

- `focus`: Optional session focus (`feature` | `bug` | `refactor` | `maintenance`)

## Execution

### Step 1: Gather context (single parallel batch)

Make ALL of these calls simultaneously in one tool-call round:

1. **Read** `/docs/planning/session-state.json` — previous session state
2. **Read** `/docs/planning/CURRENT-STATE.md` — project status summary
3. **Run** `git log --oneline -10` — recent commits
4. **Run** `git status` — working tree state
5. **Query task management** for issues assigned to me in the configured project (use project/team filters from AGENTS.md)
6. **Run** version check (include in the same parallel batch):
   ```bash
   LOCAL=$(cat .codex-workflow/VERSION 2>/dev/null || echo "unknown"); REMOTE=$(curl -sf --max-time 2 https://raw.githubusercontent.com/Clinical-Microbiomics/codex-workflow/main/base/VERSION || echo ""); echo "$LOCAL|$REMOTE"
   ```

> If `session-state.json` does not exist, skip it — you will create one at the end.

### Step 2: Display session summary

If the version check returned a remote version newer than the local version, include this at the top of the summary:
```
> Workflow update available: {local} → {remote}. Run `codex-workflow update` to update.
```
If the curl failed or versions match, show nothing — do not mention versions.

Combine the results into a single summary. Format:

```
## Session Started

**Date:** YYYY-MM-DD | **Branch:** {branch} | **Focus:** {focus or "general"}

### Previous Session
- Last: {date} — {focus}
- Completed: {completedTasks from session-state}
- Carried over: {pendingTasks from session-state}

### Active Issues (from task management)

| Priority | ID | Title | State |
|----------|----|-------|-------|
| ... | ... | ... | ... |

### Recommended First Task
{Highest priority in-progress or unblocked issue}

### Working Tree
{Clean or summary of uncommitted changes}
```

### Step 3: Write session plan file

Write the displayed summary to `/docs/planning/session-plan-YYYY-MM-DD.md` so it is available as a human-readable reference during the session. This file is cleaned up by `workflow-end-session`.

### Step 4: Update session state

Update `/docs/planning/session-state.json` with:
- New `sessionId` (format: `YYYY-MM-DD-focus`)
- `startedAt` timestamp
- `status: "active"`
- `focus` if provided
- Preserve `pendingTasks`, `backlogTasks`, `integrations`, `project`, and `testAccounts` from previous state

## What NOT to Do

These are the specific anti-patterns that cause slow execution:

- **Do NOT spawn subagents or invoke task definitions** — use direct tool calls only
- **Do NOT verify MCP connections separately** — the task management query IS the verification; if it works, the connection is healthy
- **Do NOT read AGENTS.md** — it is already loaded in your system context
- **Do NOT scan .codex-workflow/commands/, .codex-workflow/tasks/, or .codex-workflow/agents/** — not needed for session start
- **Do NOT read package.json** or check framework versions
- **Do NOT update CURRENT-STATE.md** — that happens at `workflow-end-session`
- **Do NOT run `gh auth status`** — unnecessary verification

## Output Files

- `/docs/planning/session-plan-YYYY-MM-DD.md` (human-readable session plan, deleted by `workflow-end-session`)
- `/docs/planning/session-state.json` (updated with new session)

## Error Handling

- If task management query fails: display the error, continue with local state only, note "Task management: offline" in the summary
- If session-state.json is missing: create a minimal one from AGENTS.md project metadata
- If CURRENT-STATE.md is missing: skip the "Previous Session" section

## Example

```
`workflow-start-session` focus="feature"
```

## Related

- `workflow-end-session` - Close session with handoff notes and state updates
- `workflow-session-status` - Check state without modification
- `workflow-sync-linear` - Sync with task management manually
- `workflow-setup-project-meta` - Configure project integrations (run once per project)
