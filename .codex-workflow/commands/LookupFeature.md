---
name: LookupFeature
description: Search existing backlog items by keyword or phrase
user_invocable: true
---

# /LookupFeature

Search the project backlog for existing issues matching a text query. Returns ranked results with key details.

## Purpose

Quickly find existing backlog items without leaving the CLI. Useful for:
- Checking if a feature request already exists before creating a duplicate
- Finding an issue ID to reference in `/ImplementFeature` or `/FixBug`
- Reviewing the current state of related work

## Usage

```
/LookupFeature [query] [options]
```

## Arguments

- `query`: Search text to match against issue titles and descriptions (optional, prompted if not provided)
- `state`: Filter by state — `all`, `open`, `closed`, `backlog`, `in-progress` (default: `open`)
- `limit`: Maximum number of results to return (default: 10)
- `includeCompleted`: Include completed/done issues in results (default: false)

## Execution

1. **Check configuration**: Load `project-metadata.getLinearContext`
   - If not configured, instruct user to run `/LogNewFeature` or `/SetupProjectMeta` first
   - Abort if no task management platform is set up

2. **Collect query**: If `query` not provided as argument, prompt the user for search text

3. **Search issues**: Invoke `mcp-sync.syncLinear.searchIssues` with:
   - `query`: The search text
   - `team`: From project metadata
   - `project`: From project metadata
   - `limit`: Result limit
   - `includeCompleted`: Whether to include done issues

4. **Display results**: Show matching issues in ranked table format

5. **No results handling**: If no matches found, suggest:
   - Broadening the search query
   - Using `/LogNewFeature` to create a new item

## Prerequisites

- Task management MCP configured and authenticated
- Project metadata configured (run `/SetupProjectMeta` if not)

## Output Format

### With Results
```
### Search Results for "dark mode"

Found 3 matching issues:

| # | ID        | Title                          | State       | Priority | Labels        |
|---|-----------|--------------------------------|-------------|----------|---------------|
| 1 | PROJ-456  | Add dark mode support          | Backlog     | P2       | feature       |
| 2 | PROJ-389  | Dark theme CSS variables       | In Progress | P1       | ui, feature   |
| 3 | PROJ-512  | Fix dark mode toggle on mobile | Done        | P3       | bug, ui       |

Use `/ImplementFeature linearIssue="PROJ-456"` to start working on an issue.
```

### No Results
```
### Search Results for "quantum computing"

No matching issues found.

- Try broadening your search terms
- Use `/LogNewFeature` to create a new backlog item
```

## Examples

### Basic Search
```
/LookupFeature "authentication"
```

### Search Including Completed
```
/LookupFeature query="dark mode" includeCompleted=true
```

### Search with State Filter
```
/LookupFeature query="API" state="backlog" limit=5
```

### Interactive (No Arguments)
```
/LookupFeature
```
Prompts for search query interactively.

## Error Handling

| Error | Action |
|-------|--------|
| Task management not configured | Display setup instructions, suggest `/SetupProjectMeta` |
| MCP connection failed | Display error, suggest `/CheckMCPStatus` |
| Search returned error | Display MCP error details |
| Settings file missing | Instruct user to run `/SetupProjectMeta` |

## Related Commands

- `/LogNewFeature` - Create a new backlog item
- `/SyncLinear` - Full sync operations (fetch, update, create)
- `/ImplementFeature` - Pick up an issue and implement it
- `/FixBug` - Pick up a bug and fix it
- `/SetupProjectMeta` - Configure project integrations

## Tasks Invoked

- `project-metadata.verify`
- `project-metadata.getLinearContext`
- `mcp-sync.syncLinear.searchIssues`
