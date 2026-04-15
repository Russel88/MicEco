---
description: Close engineering session with handoff notes and cleanup
---

# EndSession

Close the current engineering session with documentation updates and cleanup.

## Arguments

- `message`: Optional session summary message
- `skipSync`: Skip task management synchronization (default: false)
- `codeReview`: Trigger GitHub Actions code review workflow (default: false)
- `securityReview`: Trigger GitHub Actions security review workflow (default: false)
- `cycleUpdate`: Include Cycle status update (default: true if in active Cycle)

**Shorthand:** You can also use natural language like "with code review" or "run security review" in the command arguments.

## Execution

### Step 1: Gather session context (single parallel batch)

Make ALL of these calls simultaneously in one tool-call round:

1. **Read** `/docs/planning/session-state.json` — current session state
2. **Run** `git log --oneline -10` — recent commits
3. **Run** `git status` — working tree state
4. **Run** `git diff --stat HEAD~5..HEAD` — files changed in recent commits

### Step 2: Update task management (if not skipSync)

Update any in-progress issues in task management with comments noting progress. Only update issues that were actively worked on during the session (check commits and session state).

If in active Cycle (cycleUpdate enabled):
- Invoke `cycle-monitoring.updateProgress`
- Records session work against cycle metrics
- Updates cycle progress percentage
- Flags any new blockers against cycle scope

If `codeReview` or `securityReview` is enabled:
- Invoke `ci-integration` task with action: `triggerReview`
- Triggers the appropriate GitHub Actions workflow(s)
- Reports workflow run URL(s)

### Step 3: Display session summary

```
## Session Ended

**Date:** YYYY-MM-DD | **Duration:** {startedAt to now}

### Work Completed
- {List of commits made during session}
- {Issues progressed or completed}

### Carried Over
- {Incomplete tasks moving to next session}

### Handoff Notes
{message if provided, or auto-generated summary of what's next}

### Working Tree
{Clean or summary of uncommitted changes}
```

### Step 4: Write output files

1. **Write** `/docs/planning/session-summary-YYYY-MM-DD.md` — the displayed summary above
2. **Update** `/docs/planning/CURRENT-STATE.md` — refresh with current project status
3. **Update** `/docs/planning/session-state.json`:
   - Set `status: "ended"`, `endedAt` timestamp
   - Move active work to `pendingTasks`
   - Record `completedTasks` from this session
   - Clear `currentTask`

### Step 5: Clean up session plan

**Delete** `/docs/planning/session-plan-*.md` files from the current session. These are temporary working files created by `/StartSession` and are superseded by the session summary.

## What NOT to Do

- **Do NOT spawn subagents or invoke task definitions** — use direct tool calls only
- **Do NOT verify MCP connections** — if sync fails, report the error and continue
- **Do NOT read CLAUDE.md** — it is already loaded in your system context
- **Do NOT run `gh auth status`** — unnecessary verification

## Output Files

- `/docs/planning/session-summary-YYYY-MM-DD.md` (created)
- `/docs/planning/session-state.json` (updated)
- `/docs/planning/CURRENT-STATE.md` (updated)
- `/docs/planning/session-plan-*.md` (deleted — cleanup)

## Error Handling

- If task management sync fails: report the error, continue with local-only cleanup
- If session-state.json is missing: create a minimal ended state from git context
- If no commits were made: note "No commits this session" in the summary

## Example

```bash
/EndSession message="Completed certificate feature, ready for testing"
```

## Related

- `/StartSession` - Initialize a session
- `/SessionStatus` - Check session state
- `/SyncLinear` - Sync with task management manually

## Tasks Invoked

- `session-management.end`
- `mcp-sync.syncLinear`
- `mcp-sync.syncCoda`
- `git-workflow.getState`
- `documentation-sync.updateCurrentState`
- `ci-integration.triggerReview` (if reviews enabled)
- `cycle-monitoring.updateProgress` (if in active cycle)
