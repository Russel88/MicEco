---
name: workflow-log-new-feature
description: Quickly log a new feature, bug, or task to the project backlog Use when the user asks to log new feature in this repository or refers to the old `/LogNewFeature` command.
---

# Log New Feature

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# `workflow-log-new-feature`

Quickly create a new backlog item (feature, bug, or task) in the configured task management platform without starting a full session.

## Purpose

Provides a lightweight way to capture ideas, bugs, and tasks directly into the project backlog. On first use, guides the user through task management platform setup if not already configured.

## Usage

```
`workflow-log-new-feature` [title] [options]
```

## Arguments

- `title`: Item title (optional, prompted if not provided)
- `description`: Item description (optional)
- `type`: Item type — `feature`, `bug`, or `task` (optional, prompted if not provided)
- `priority`: Priority level — `urgent`, `high`, `medium`, `low`, `none` (optional)
- `labels`: Comma-separated labels to apply (optional)

## Execution

### First-Time Setup Check

1. Load `.codex-workflow/settings.json`
2. Check for `projectMetadata.linear` configuration
3. If not configured:
   - Inform user that task management is not set up
   - Invoke `project-metadata.setup` with `--linear-only` to configure
   - If setup is declined or fails, abort with instructions to run `workflow-setup-project-meta`

### Item Creation Flow

1. **Determine item type**: If `type` not provided, ask the user:
   - Feature — New functionality or enhancement
   - Bug — Defect or regression
   - Task — Maintenance, chore, or operational work

2. **Collect title**: If `title` not provided as argument, prompt the user for a concise title

3. **Collect description** (optional): Ask if user wants to add a description
   - If yes, prompt for description text
   - If no, proceed without description

4. **Collect priority** (optional): If `priority` not provided, ask the user:
   - Urgent (P0)
   - High (P1)
   - Medium (P2)
   - Low (P3)
   - No priority

5. **Create issue**: Invoke `mcp-sync.syncLinear.createIssue` with:
   - `title`: The item title
   - `description`: The description (if provided)
   - `team`: From `project-metadata.getLinearContext`
   - `labels`: The item type as a label (e.g., `feature`, `bug`, `task`) plus any user-specified labels
   - `priority`: Mapped priority level

6. **Display result**: Show created item details

## Prerequisites

- Task management MCP configured and authenticated
- `.codex-workflow/settings.json` exists (created by `workflow-setup-workflow` or on first run)

## Output Format

```
### Backlog Item Created

| Field       | Value                              |
|-------------|------------------------------------|
| ID          | PROJ-456                           |
| Type        | Feature                            |
| Title       | Add dark mode support              |
| Priority    | P2 (Medium)                        |
| State       | Backlog                            |
| Labels      | feature                            |
| Link        | https://linear.app/team/issue/...   |
```

## Examples

### With Arguments
```
`workflow-log-new-feature` title="Add dark mode support" type="feature" priority="medium"
```

### Interactive (No Arguments)
```
`workflow-log-new-feature`
```
Prompts for type, title, description, and priority interactively.

### Quick Bug Report
```
`workflow-log-new-feature` title="Login button unresponsive on mobile" type="bug" priority="high"
```

### Minimal
```
`workflow-log-new-feature` "Refactor auth middleware"
```
Prompts for type and priority; uses provided string as title.

## Error Handling

| Error | Action |
|-------|--------|
| Task management not configured | Trigger `project-metadata.setup` |
| MCP connection failed | Display error, suggest `workflow-check-mcp-status` |
| Issue creation failed | Display MCP error details, suggest retry |
| Settings file missing | Create default settings, trigger setup |

## Related Commands

- `workflow-lookup-feature` - Search existing backlog items
- `workflow-sync-linear` - Full sync operations (fetch, update, create)
- `workflow-setup-project-meta` - Configure project integrations
- `workflow-implement-feature` - Pick up a backlog item and implement it
- `workflow-check-mcp-status` - Verify MCP connectivity

## Tasks Invoked

- `project-metadata.verify`
- `project-metadata.setup` (first-time only)
- `project-metadata.getLinearContext`
- `mcp-sync.syncLinear.createIssue`
