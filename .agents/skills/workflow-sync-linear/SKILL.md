---
name: workflow-sync-linear
description: Manually sync with task management for issue management Use when the user asks to sync linear in this repository or refers to the old `/SyncLinear` command.
---

# Sync Linear

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# SyncLinear

Manual synchronization with task management system for issue fetching and updates.

## Purpose

Provides direct access to task management MCP operations outside of session management, useful for quick updates or when not in an active session.

## Arguments

- `action`: Sync action (required)
  - `fetch` - Fetch issues
  - `update` - Update an issue
  - `create` - Create new issue
- `team`: Team name or ID (for fetch)
- `issueId`: Issue ID (for update)
- `state`: New state (for update)
- `comment`: Comment to add (for update/create)
- `title`: Issue title (for create)
- `description`: Issue description (for create)

## Execution

1. Invoke `mcp-sync` task with appropriate action
   - For `fetch`: Query for issues
   - For `update`: Update issue state/add comment
   - For `create`: Create new issue

2. Display:
   - Action result
   - Issue details
   - Confirmation

## Prerequisites

- Task management MCP configured and authenticated

## Examples

### Fetch Issues
```
`workflow-sync-linear` action="fetch" team="Core"
```

### Update Issue
```
`workflow-sync-linear` action="update" issueId="PROJ-211" state="In Review" comment="Ready for code review"
```

### Create Issue
```
`workflow-sync-linear` action="create" team="Core" title="Fix timeout bug" description="Handle timeout errors gracefully"
```

## Output Format

### Fetch
```
### Issues (Team: Core)
| ID | Title | State | Priority |
|----|-------|-------|----------|
| PROJ-211 | Feature implementation | In Progress | P1 |
| PROJ-212 | Bug fix | Backlog | P2 |
```

### Update
```
### Issue Updated
- ID: PROJ-211
- New State: In Review
- Comment Added: Yes
```

## Related

- `workflow-start-session` - Full session with sync
- `workflow-end-session` - Session end with updates
- `workflow-check-mcp-status` - Verify connection

## Tasks Invoked

- `mcp-sync.syncLinear.fetchIssues`
- `mcp-sync.syncLinear.updateIssue`
- `mcp-sync.syncLinear.createIssue`
