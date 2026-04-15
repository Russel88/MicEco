---
description: Manage project documentation - audit, sync, and update docs
---

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
/UpdateDocs action="audit" scope="all"
```

### Sync with External System
```
/UpdateDocs action="sync" direction="pull"
```

### Archive Old Sessions
```
/UpdateDocs action="archive"
```

### Generate Implementation Map
```
/UpdateDocs action="map"
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

- `/StartSession` - Creates session plans
- `/EndSession` - Creates session summaries
- `/CheckMCPStatus` - Verify external system connection

## Tasks Invoked

- `documentation-sync.audit`
- `documentation-sync.reconcile`
- `documentation-sync.archiveSessions`
- `documentation-sync.generateImplementationMap`
