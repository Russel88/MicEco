---
name: workflow-update-docs
description: Manage project documentation - audit, sync, and update docs Use when the user asks to update docs in this repository or refers to the old `/UpdateDocs` command.
---

# Update Docs

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# UpdateDocs

Manage project documentation including auditing, syncing with external systems, and updating local files.

## Purpose

Provides documentation management capabilities for auditing document health, syncing with external systems, and maintaining documentation organization.

## Arguments

- `action`: Documentation action (required)
  - `audit` - Audit documentation health
  - `sync` - Sync with external documentation (e.g., Coda)
  - `archive` - Archive old session files
  - `map` - Generate implementation map
- `scope`: Audit scope (for audit action)
  - `all` - Full audit
  - `docs` - Only `/docs`
  - `knowledge` - Only `/knowledge`
  - `planning` - Only `/docs/planning`
- `direction`: Sync direction (for sync action)
  - `pull` - External → Local
  - `push` - Local → External

## Execution

1. Invoke `documentation-sync` task with appropriate action
   - For `audit`: Scan and analyze documents
   - For `sync`: Reconcile with external system
   - For `archive`: Move old session files
   - For `map`: Generate implementation map

2. Display:
   - Action results
   - Findings (for audit)
   - Sync operations (for sync)
   - Files moved (for archive)

## Prerequisites

- Documentation folders exist
- Documentation MCP configured (for sync)

## Examples

### Audit Documentation
```
`workflow-update-docs` action="audit" scope="all"
```

### Sync with External System
```
`workflow-update-docs` action="sync" direction="pull"
```

### Archive Old Sessions
```
`workflow-update-docs` action="archive"
```

### Generate Implementation Map
```
`workflow-update-docs` action="map"
```

## Output Format

### Audit
```
### Documentation Audit
Total Documents: 45

#### Findings
- Stale (>30 days): 5
- Orphaned (no references): 3
- Missing specs: 2
- Sync conflicts: 1
- Healthy: 34

#### Recommendations
1. Archive: /docs/old-feature/implementation.md
2. Create spec: Notification System
3. Sync from external: PRD updates
```

### Archive
```
### Session Archive
Threshold: 7 days

#### Archived Files
- session-summary-2026-01-05.md → archive/2026-01/
- session-plan-2026-01-05.md → archive/2026-01/

Summary: 5 files archived, 3 remaining
```

## Output Files

- `/docs/reports/doc-audit-YYYY-MM-DD.md` (for audit)
- `/knowledge/architecture/IMPLEMENTATION-MAP.md` (for map)

## Related

- `workflow-start-session` - Creates session plans
- `workflow-end-session` - Creates session summaries
- `workflow-check-mcp-status` - Verify external system connection

## Tasks Invoked

- `documentation-sync.audit`
- `documentation-sync.reconcile`
- `documentation-sync.archiveSessions`
- `documentation-sync.generateImplementationMap`
