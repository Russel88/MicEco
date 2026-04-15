---
description: Manually sync with task management for issue management
---

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
/SyncLinear action="fetch" team="Core"
```

### Update Issue
```
/SyncLinear action="update" issueId="PROJ-211" state="In Review" comment="Ready for code review"
```

### Create Issue
```
/SyncLinear action="create" team="Core" title="Fix timeout bug" description="Handle timeout errors gracefully"
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

- `/StartSession` - Full session with sync
- `/EndSession` - Session end with updates
- `/CheckMCPStatus` - Verify connection

## Tasks Invoked

- `mcp-sync.syncLinear.fetchIssues`
- `mcp-sync.syncLinear.updateIssue`
- `mcp-sync.syncLinear.createIssue`
