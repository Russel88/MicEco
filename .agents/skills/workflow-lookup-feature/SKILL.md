---
name: workflow-lookup-feature
description: Search existing backlog items by keyword or phrase Use when the user asks to lookup feature in this repository or refers to the old `/LookupFeature` command.
---

# Lookup Feature

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# `workflow-lookup-feature`

Search the project backlog for existing issues matching a text query. Returns ranked results with key details.

## Purpose

Quickly find existing backlog items without leaving the CLI. Useful for:
- Checking if a feature request already exists before creating a duplicate
- Finding an issue ID to reference in `workflow-implement-feature` or `workflow-fix-bug`
- Reviewing the current state of related work

## Usage

```
`workflow-lookup-feature` [query] [options]
```

## Arguments

- `query`: Search text to match against issue titles and descriptions (optional, prompted if not provided)
- `state`: Filter by state — `all`, `open`, `closed`, `backlog`, `in-progress` (default: `open`)
- `limit`: Maximum number of results to return (default: 10)
- `includeCompleted`: Include completed/done issues in results (default: false)

## Execution

1. **Check configuration**: Load `project-metadata.getLinearContext`
   - If not configured, instruct user to run `workflow-log-new-feature` or `workflow-setup-project-meta` first
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
   - Using `workflow-log-new-feature` to create a new item

## Prerequisites

- Task management MCP configured and authenticated
- Project metadata configured (run `workflow-setup-project-meta` if not)

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

Use ``workflow-implement-feature` linearIssue="PROJ-456"` to start working on an issue.
```

### No Results
```
### Search Results for "quantum computing"

No matching issues found.

- Try broadening your search terms
- Use `workflow-log-new-feature` to create a new backlog item
```

## Examples

### Basic Search
```
`workflow-lookup-feature` "authentication"
```

### Search Including Completed
```
`workflow-lookup-feature` query="dark mode" includeCompleted=true
```

### Search with State Filter
```
`workflow-lookup-feature` query="API" state="backlog" limit=5
```

### Interactive (No Arguments)
```
`workflow-lookup-feature`
```
Prompts for search query interactively.

## Error Handling

| Error | Action |
|-------|--------|
| Task management not configured | Display setup instructions, suggest `workflow-setup-project-meta` |
| MCP connection failed | Display error, suggest `workflow-check-mcp-status` |
| Search returned error | Display MCP error details |
| Settings file missing | Instruct user to run `workflow-setup-project-meta` |

## Related Commands

- `workflow-log-new-feature` - Create a new backlog item
- `workflow-sync-linear` - Full sync operations (fetch, update, create)
- `workflow-implement-feature` - Pick up an issue and implement it
- `workflow-fix-bug` - Pick up a bug and fix it
- `workflow-setup-project-meta` - Configure project integrations

## Tasks Invoked

- `project-metadata.verify`
- `project-metadata.getLinearContext`
- `mcp-sync.syncLinear.searchIssues`
