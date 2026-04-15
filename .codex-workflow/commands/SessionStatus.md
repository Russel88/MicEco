---
description: Check current session state without modification
---

# SessionStatus

View the current session state and project status without making changes.

## Purpose

Provides a quick view of the current session state, active work, blockers, and MCP health without modifying any state.

## Arguments

- `verbose`: Show detailed information (default: false)

## Execution

1. Invoke `session-management` task with action: `status`
   - Loads session state
   - Formats status report
   - Returns without modification

2. Display:
   - Session overview
   - Active work items
   - Blockers
   - MCP health
   - Quick links

## Prerequisites

- None (read-only operation)

## Output Format

```
### Session Status

#### Current Session
- ID: session-2026-01-16-1000
- Started: 2026-01-16 10:00 UTC
- Duration: 2h 30m
- Focus: feature
- Branch: feature/PROJ-211-feature-name

#### Active Work
- [In Progress] Feature implementation
  - Files: lib/feature/implementation.ts
  - Next: Test edge cases

#### Blockers
None currently.

#### MCP Health
- Task Management: OK (synced 5m ago)
- Documentation: OK (synced 5m ago)
- GitHub CLI: OK

#### Recent Decisions
- 2026-01-16: Architecture decision notes

#### Quick Links
- Session State: /docs/planning/session-state.json
- Current State: /docs/planning/CURRENT-STATE.md
- Latest Plan: /docs/planning/session-plan-2026-01-16.md
```

## Example

```
/SessionStatus verbose=true
```

## Verbose Mode

With `verbose=true`, includes:
- Full list of assigned issues
- Detailed MCP connection info
- Git repository state
- All recent decisions

## Related

- `/StartSession` - Start a session
- `/EndSession` - End a session
- `/CheckMCPStatus` - Detailed MCP status
- `/SyncLinear` - Sync with task management

## Tasks Invoked

- `session-management.status`

## Notes

This command is read-only and never modifies session state or any files. Safe to run at any time.
