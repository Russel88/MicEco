---
name: mcp-sync
description: Centralized MCP interaction layer for task management and documentation synchronization
---

# MCP Sync Task

Unified interface for all MCP (Model Context Protocol) interactions. Centralizes task management and documentation operations to ensure consistent error handling, health tracking, and state management.

## Operations

### `verifyConnections`

Check health and authentication status of all MCPs.

**Steps:**
1. Test task management MCP connection:
   - Attempt list operation
   - Verify authentication
   - Check capability schema
2. Test documentation MCP connection:
   - Attempt list operation
   - Verify authentication
3. Test GitHub CLI (not MCP but required):
   - Run `gh auth status`
4. Return status object

**Outputs:**
```json
{
  "taskManagement": {
    "healthy": true,
    "lastChecked": "2026-01-16T10:00:00Z",
    "error": null
  },
  "documentation": {
    "healthy": true,
    "lastChecked": "2026-01-16T10:00:00Z",
    "error": null
  },
  "github": {
    "healthy": true,
    "authenticated": true,
    "error": null
  }
}
```

### `syncLinear`

Synchronize with task management for issue management.

**Sub-operations:**

#### `fetchIssues`
Fetch assigned issues.

**Inputs:**
- `team`: Team name or ID (optional)
- `assignee`: Filter by assignee (default: "me")
- `state`: Filter by state (default: open/active)
- `limit`: Max issues to return (default: 50)

**Steps:**
1. Call task management MCP with filters
2. Parse response into standardized format
3. Return issue list

**Outputs:**
```json
{
  "issues": [
    {
      "id": "PROJ-123",
      "title": "Issue title",
      "state": "In Progress",
      "priority": 2,
      "assignee": "User",
      "labels": ["feature"],
      "dueDate": "2026-01-20"
    }
  ],
  "syncedAt": "2026-01-16T10:00:00Z"
}
```

#### `updateIssue`
Update an issue.

**Inputs:**
- `id`: Issue ID
- `state`: New state (optional)
- `comment`: Comment to add (optional)
- `labels`: Labels to set (optional)

**Steps:**
1. Update issue if state/labels changed
2. Add comment if provided
3. Return confirmation

#### `createIssue`
Create a new issue.

**Inputs:**
- `title`: Issue title
- `description`: Issue description
- `team`: Team name or ID
- `assignee`: Assignee (optional)
- `labels`: Labels (optional)
- `priority`: Priority level (optional)

**Steps:**
1. Create issue via MCP
2. Return created issue details

#### `searchIssues`
Search issues by text query against titles and descriptions.

**Inputs:**
- `query`: Search text to match against issue titles and descriptions
- `team`: Team name or ID (optional, loaded from `project-metadata.getLinearContext` if not provided)
- `project`: Project name or ID (optional, loaded from `project-metadata.getLinearContext` if not provided)
- `limit`: Max issues to return (default: 10)
- `includeCompleted`: Include completed/done issues (default: false)

**Steps:**
1. Load project context from `project-metadata.getLinearContext` if team/project not provided
2. Call task management MCP to list issues within the team/project scope
3. Filter results by matching `query` text against issue titles and descriptions
4. Score results by relevance:
   - Title exact match: highest weight
   - Title partial/word match: high weight
   - Description match: moderate weight
5. Sort results by relevance score (descending)
6. Truncate to `limit` results
7. Return ranked issue list

**Outputs:**
```json
{
  "query": "dark mode",
  "totalMatches": 3,
  "issues": [
    {
      "id": "PROJ-456",
      "title": "Add dark mode support",
      "state": "Backlog",
      "priority": 2,
      "labels": ["feature"],
      "relevanceScore": 0.95,
      "matchedIn": ["title"]
    },
    {
      "id": "PROJ-389",
      "title": "CSS variable refactor",
      "state": "In Progress",
      "priority": 1,
      "labels": ["ui", "feature"],
      "relevanceScore": 0.62,
      "matchedIn": ["description"]
    }
  ],
  "searchedAt": "2026-01-27T10:00:00Z"
}
```

### `syncCoda`

Synchronize with documentation system.

**Sub-operations:**

#### `fetchPRD`
Fetch PRD content from documentation system.

**Inputs:**
- `docId`: Document ID
- `pageId`: Page ID or name

**Steps:**
1. Fetch page content via MCP
2. Parse markdown content
3. Return PRD object

**Outputs:**
```json
{
  "title": "PRD Title",
  "content": "Markdown content...",
  "lastUpdated": "2026-01-15",
  "source": "documentation"
}
```

#### `updatePage`
Update a documentation page.

**Inputs:**
- `docId`: Document ID
- `pageId`: Page ID or name
- `content`: New markdown content

**Steps:**
1. Update page via MCP
2. Return confirmation

#### `appendToPage`
Append content to a page.

**Inputs:**
- `docId`: Document ID
- `pageId`: Page ID or name
- `content`: Content to append

**Steps:**
1. Append content via MCP
2. Return confirmation

### `createInitiative`

Create an Initiative in Linear from PRD data.

**Inputs:**
- `title`: Initiative title
- `description`: Initiative description (from PRD)
- `targetDate`: Target completion date (optional)

**Steps:**
1. Create Initiative via Linear MCP
2. Return Initiative ID and URL

**Outputs:**
```json
{ "initiativeId": "init-123", "url": "https://linear.app/..." }
```

### `createWorkItems`

Batch create Work Items from breakdown.

**Inputs:**
- `items`: Array of work item specs (title, description, team, labels, priority, estimate)
- `projectId`: Linear project ID

**Steps:**
1. Iterate over items
2. Create each issue via Linear MCP
3. Set dependency relationships between items
4. Return created issue IDs

**Outputs:**
```json
{ "created": [{ "id": "PROJ-101", "title": "..." }], "count": 12 }
```

### `createCycle`

Create a Cycle in Linear.

**Inputs:**
- `teamId`: Team ID
- `name`: Cycle name
- `startDate`: Start date
- `endDate`: End date

**Steps:**
1. Create Cycle via Linear MCP
2. Return Cycle ID

**Outputs:**
```json
{ "cycleId": "cycle-123", "url": "https://linear.app/..." }
```

### `assignToCycle`

Assign Work Items to a Cycle.

**Inputs:**
- `cycleId`: Cycle ID
- `issueIds`: Array of issue IDs to assign

**Steps:**
1. Update each issue to assign to cycle
2. Return confirmation

**Outputs:**
```json
{ "assigned": 12, "cycleId": "cycle-123" }
```

### `getCycleProgress`

Fetch Cycle status from Linear.

**Inputs:**
- `cycleId`: Cycle ID (optional, defaults to current cycle)
- `teamId`: Team ID

**Steps:**
1. Fetch cycle details from Linear MCP
2. Fetch all issues in cycle
3. Calculate progress metrics
4. Return progress object

**Outputs:**
```json
{ "cycleId": "cycle-123", "total": 12, "completed": 5, "inProgress": 4, "blocked": 1, "remaining": 2 }
```

### `getBacklog`

Fetch prioritized backlog for planning.

**Inputs:**
- `teamId`: Team ID (optional)
- `projectId`: Project ID (optional)
- `state`: Filter by state (default: "backlog")
- `limit`: Max items (default: 100)

**Steps:**
1. Fetch backlog issues from Linear MCP sorted by priority
2. Return prioritized list

**Outputs:**
```json
{ "items": [{ "id": "PROJ-50", "title": "...", "priority": 1, "estimate": 8 }], "total": 45 }
```

## Configuration

The task works with any configured MCP servers for:
- Task management (Linear, Jira, etc.)
- Documentation (Coda, Notion, etc.)

**Project Scoping:** When `project-metadata` is configured, operations are automatically scoped:
- Linear operations use the configured team and project
- Coda operations use the configured document and page
- Call `project-metadata.getLinearContext` or `project-metadata.getCodaContext` to get these values

## Error Handling

| Error Type | Action |
|------------|--------|
| Authentication failure | Return `healthy: false` with error details |
| Rate limit | Retry with exponential backoff (max 3 attempts) |
| Network error | Mark unhealthy, allow offline mode |
| Invalid response | Log error, return partial data if possible |

## Health Tracking

Health status is stored in `session-state.json`:
```json
{
  "mcpState": {
    "taskManagementLastSync": "2026-01-16T10:05:00Z",
    "documentationLastSync": "2026-01-16T10:05:00Z",
    "taskManagementHealthy": true,
    "documentationHealthy": true
  }
}
```

## Dependencies

- **project-metadata** (optional): For getting configured team/project/document context

## Usage Examples

```
// Verify all connections
invoke mcp-sync.verifyConnections

// Fetch my open issues
invoke mcp-sync.syncLinear.fetchIssues(assignee="me", state="open")

// Update issue status
invoke mcp-sync.syncLinear.updateIssue(id="PROJ-123", state="Done", comment="Completed in session")

// Search for issues by keyword
invoke mcp-sync.syncLinear.searchIssues(query="dark mode", limit=10)

// Search including completed issues
invoke mcp-sync.syncLinear.searchIssues(query="auth", includeCompleted=true)

// Fetch PRD from documentation
invoke mcp-sync.syncCoda.fetchPRD(docId="doc123", pageId="PRD Page")
```
