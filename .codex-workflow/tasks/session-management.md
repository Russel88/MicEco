---
name: session-management
description: Orchestrate session lifecycle including initialization, state management, and cleanup
---

# Session Management Task

Centralized orchestration for engineering session lifecycle. Manages state persistence, MCP synchronization, and context handoff between sessions.

## Operations

### `start`

Initialize a new engineering session.

**Steps:**
1. Load `/docs/planning/session-state.json`
2. Validate schema version (warn if outdated)
3. Check staleness (warn if `lastUpdated` > 24h ago)
4. Invoke `project-metadata` task with action: `verify`
   - Check if Linear, Coda, GitHub metadata is configured
   - If not configured, prompt to run `/SetupProjectMeta`
   - Verify configured integrations are accessible
5. Invoke `mcp-sync` task with action: `verifyConnections`
6. Invoke `context-loader` task to gather project context
7. Invoke `priority-matrix` task to build task list (using configured project)
8. Initialize new session in state:
   - Generate session ID: `session-YYYY-MM-DD-HHMM`
   - Set `startedAt` timestamp
   - Clear previous session's `inProgress` items
   - Preserve `handoff` data for reference
9. Generate session plan document: `/docs/planning/session-plan-YYYY-MM-DD.md`
10. Save updated state

**Outputs:**
- Session plan markdown file
- Updated `session-state.json`
- Priority matrix summary

### `end`

Close the current engineering session.

**Inputs:**
- `message`: Optional session summary message
- `skipSync`: Skip MCP synchronization (default: false)
- `codeReview`: Trigger GitHub Actions code review (default: false)
- `securityReview`: Trigger GitHub Actions security review (default: false)

**Steps:**
1. Load current `session-state.json`
2. Collect work summaries from session
3. Invoke `mcp-sync` task with action: `syncLinear` to update issue statuses
4. Invoke `mcp-sync` task with action: `syncCoda` if documentation changed
5. Invoke `git-workflow` task to capture repository state
6. Populate `handoff` section:
   - Move incomplete `inProgress` items to `handoff.stashedWork`
   - Capture untracked files
   - Generate `nextSessionPrereqs`
7. Generate session summary: `/docs/planning/session-summary-YYYY-MM-DD-HHMM.md`
8. Update `CURRENT-STATE.md` with latest status
9. If `codeReview` or `securityReview` requested:
   - Invoke `ci-integration` task with action: `triggerReview`
   - Pass `codeReview` and `securityReview` flags
   - Record workflow URLs in session summary
10. Clear `currentSession` object
11. Save updated state

**Outputs:**
- Session summary markdown file
- Updated `session-state.json`
- Updated `CURRENT-STATE.md`
- GitHub Actions workflow URLs (if CI triggered)

### `status`

Report current session state without modification.

**Steps:**
1. Load `session-state.json`
2. Format human-readable status:
   - Session duration
   - Active work items
   - Blockers
   - MCP health
3. Return formatted status

**Outputs:**
- Status summary (no file changes)

### `updateProgress`

Update in-progress work during a session.

**Steps:**
1. Load `session-state.json`
2. Update specified `inProgress` item
3. Add to `recentDecisions` if decision provided
4. Save state

**Inputs:**
- `task`: Task identifier
- `files`: Array of modified files
- `nextStep`: Description of next action
- `decision`: Optional decision record

**Outputs:**
- Updated `session-state.json`

### `addBlocker`

Record a blocker encountered during session.

**Steps:**
1. Load `session-state.json`
2. Add blocker to `activeWork.blockers`
3. Save state

**Inputs:**
- `description`: Blocker description
- `linkedIssue`: Optional issue ID

**Outputs:**
- Updated `session-state.json`

## State Schema

See `/docs/planning/session-state.json` for full schema.

## Dependencies

- **project-metadata**: For verifying project integrations (Linear, Coda, GitHub)
- **mcp-sync**: For task management/documentation interactions
- **context-loader**: For project context
- **priority-matrix**: For task prioritization
- **git-workflow**: For repository state
- **documentation-sync**: For doc reconciliation

## Error Handling

- If `session-state.json` missing: Create with default empty state
- If MCPs unavailable: Continue in offline mode, mark `mcpState.healthy: false`
- If state schema outdated: Migrate to current version

## File Outputs

| Operation | File Created |
|-----------|--------------|
| start | `/docs/planning/session-plan-YYYY-MM-DD.md` |
| end | `/docs/planning/session-summary-YYYY-MM-DD-HHMM.md` |
| end | `/docs/planning/CURRENT-STATE.md` (updated) |
