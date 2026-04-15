---
name: documentation-sync
description: Reconcile local documentation with external systems and identify gaps
---

# Documentation Sync Task

Manages documentation synchronization between local files (`/docs`, `/knowledge`) and external systems. Identifies gaps, staleness, and inconsistencies.

## Operations

### `audit`

Full documentation audit.

**Inputs:**
- `scope`: Audit scope
  - `all` - Full audit
  - `docs` - Only `/docs`
  - `knowledge` - Only `/knowledge`
  - `planning` - Only `/docs/planning`

**Steps:**
1. Scan local documentation:
   - List all files in `/docs` and `/knowledge`
   - Extract metadata (title, last modified, type)
   - Identify document types (PRD, spec, plan, report)
2. Check external systems for corresponding documents:
   - Query documentation MCP
   - Compare timestamps
   - Identify missing local copies
3. Check task management for documentation references:
   - Scan issue descriptions for doc links
   - Identify broken links
   - Find undocumented features
4. Analyze findings:
   - Stale documents (>30 days unmodified)
   - Orphaned documents (no references)
   - Missing documentation (features without specs)
   - Conflicting versions (local vs external)
5. Generate audit report

**Outputs:**
```json
{
  "audit": {
    "scope": "all",
    "timestamp": "2026-01-16T10:00:00Z",
    "totalDocuments": 45
  },
  "findings": {
    "stale": [...],
    "orphaned": [...],
    "missing": [...],
    "conflicts": [...]
  },
  "summary": {
    "stale": 5,
    "orphaned": 3,
    "missing": 2,
    "conflicts": 1,
    "healthy": 34
  }
}
```

### `reconcile`

Sync local docs with external system.

**Inputs:**
- `direction`: Sync direction
  - `pull` - External → Local
  - `push` - Local → External
  - `bidirectional` - Merge changes
- `documents`: Specific documents to sync (optional)

**Steps:**
1. Identify documents to sync
2. Compare versions
3. For `pull`:
   - Fetch content via `mcp-sync`
   - Update local files
4. For `push`:
   - Read local content
   - Update external via `mcp-sync`
5. For `bidirectional`:
   - Identify conflicts
   - Present for manual resolution
6. Log sync operations

**Outputs:**
```json
{
  "sync": {
    "direction": "pull",
    "timestamp": "2026-01-16T10:00:00Z"
  },
  "operations": [...],
  "conflicts": [],
  "summary": {
    "synced": 3,
    "skipped": 0,
    "conflicts": 0
  }
}
```

### `archiveSessions`

Archive old session files.

**Inputs:**
- `olderThan`: Days threshold (default: 7)
- `dryRun`: Preview without moving (default: false)

**Steps:**
1. Scan `/docs/planning/` for session files
2. Identify files older than threshold
3. Create archive folder structure: `/docs/planning/archive/YYYY-MM/`
4. Move old files to archive
5. Update any references

**Outputs:**
```json
{
  "archive": {
    "threshold": 7,
    "dryRun": false
  },
  "moved": [...],
  "summary": {
    "archived": 5,
    "remaining": 3
  }
}
```

### `updateCurrentState`

Update CURRENT-STATE.md from session state.

**Inputs:**
- None (reads from session-state.json)

**Steps:**
1. Load `session-state.json`
2. Generate human-readable summary
3. Write to `/docs/planning/CURRENT-STATE.md`

**Outputs:**
```json
{
  "updated": true,
  "path": "/docs/planning/CURRENT-STATE.md"
}
```

### `generateImplementationMap`

Create or update implementation map.

**Inputs:**
- `features`: Features to map (optional, default: all)

**Steps:**
1. Scan `/knowledge/prd/` for PRDs
2. Find corresponding specs in `/docs/specs/`
3. Find related code locations
4. Find related issues
5. Generate map document

**Outputs:**
```json
{
  "map": {
    "features": [...]
  },
  "path": "/knowledge/architecture/IMPLEMENTATION-MAP.md"
}
```

## Document Types

| Type | Location | Staleness Threshold |
|------|----------|---------------------|
| PRD | `/knowledge/prd/` | 90 days |
| Feature Spec | `/docs/specs/` | 30 days |
| Session Plan | `/docs/planning/` | 7 days (then archive) |
| Session Summary | `/docs/planning/` | 7 days (then archive) |
| Report | `/docs/reports/` | Never stale |
| Architecture | `/knowledge/architecture/` | 90 days |

## Dependencies

- **mcp-sync**: For external system operations
- **session-management**: For state access

## File Operations

| Operation | Files Modified |
|-----------|----------------|
| audit | None (read-only) |
| reconcile | Local docs or external pages |
| archiveSessions | Session files |
| updateCurrentState | CURRENT-STATE.md |
| generateImplementationMap | IMPLEMENTATION-MAP.md |

## Error Handling

| Error | Action |
|-------|--------|
| External system unavailable | Audit local only, flag sync issues |
| File permission denied | Skip file, log error |
| Conflict detected | Flag for manual resolution |
| Invalid markdown | Best-effort parse, log warning |
