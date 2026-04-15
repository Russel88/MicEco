---
name: work-breakdown
description: Break down PRD requirements into sized work items with dependencies
---

# Work Breakdown Task

Transforms enriched and validated PRD requirements into implementable, sized work items with dependency graphs. Produces ready-to-import task lists for cycle planning.

## Operations

### `generate`

Create work items from PRD requirements.

**Steps:**
1. Load enriched, validated PRD
2. Parse each requirement into implementable units
3. Size each unit to target (2-4h, 4-8h, 1-2d) - split if too large
4. Categorize by area (backend, frontend, database, infrastructure, testing)
5. Assign priority based on dependencies and risk
6. Build dependency graph between items
7. Add test items for each feature item if `includeTests=true`
8. Generate work item specs with title, description, acceptance criteria, size estimate
9. Store breakdown document

**Inputs:**
- `prd`: PRD slug or ID (string)
- `targetSize`: Target size for work items (enum: `2-4h` | `4-8h` | `1-2d`)
- `includeTests`: Include test items for each feature item (boolean, default: true)

**Outputs:**
```json
{
  "prd": "user-auth-revamp",
  "items": [
    {
      "id": "WB-001",
      "title": "Add Redis session store adapter",
      "description": "Implement session storage adapter using Redis to replace in-memory store",
      "area": "backend",
      "size": "4-8h",
      "priority": 1,
      "dependsOn": [],
      "acceptanceCriteria": [
        "Sessions persist across server restarts",
        "TTL-based session expiry works correctly",
        "Fallback to in-memory if Redis unavailable"
      ]
    },
    {
      "id": "WB-002",
      "title": "Implement SSO login endpoint",
      "description": "Create OAuth2 callback endpoint for SSO provider integration",
      "area": "backend",
      "size": "4-8h",
      "priority": 2,
      "dependsOn": ["WB-001"],
      "acceptanceCriteria": [
        "OAuth2 authorization code flow completes successfully",
        "User profile mapped from SSO claims",
        "Error handling for denied/expired tokens"
      ]
    },
    {
      "id": "WB-003",
      "title": "Test: Redis session store adapter",
      "description": "Unit and integration tests for Redis session adapter",
      "area": "testing",
      "size": "2-4h",
      "priority": 2,
      "dependsOn": ["WB-001"],
      "acceptanceCriteria": [
        "Unit tests for CRUD operations",
        "Integration test with Redis container",
        "Edge cases: connection failure, TTL expiry"
      ]
    }
  ],
  "totalEstimate": {
    "min": "8 days",
    "max": "12 days"
  },
  "dependencyGraph": "WB-001 --> WB-002 --> WB-004\n  |            |\n  v            v\nWB-003      WB-005"
}
```

### `preview`

Dry run showing what would be created without creating tickets.

**Steps:**
1. Run generate logic
2. Return preview without MCP calls

**Inputs:**
- `prd`: PRD slug or ID (string)
- `targetSize`: Target size for work items (enum: `2-4h` | `4-8h` | `1-2d`)

**Outputs:**
```json
{
  "prd": "user-auth-revamp",
  "dryRun": true,
  "items": [],
  "totalEstimate": { "min": "8 days", "max": "12 days" },
  "dependencyGraph": "..."
}
```

## Configuration

Size targets define the granularity of work items:
- **2-4h**: Fine-grained, ideal for junior developers or tightly tracked sprints
- **4-8h**: Standard, good balance of granularity and overhead
- **1-2d**: Coarse-grained, suitable for experienced teams with high trust

Priority assignment follows dependency-first ordering:
1. Items with no dependencies (can start immediately)
2. Items on the critical path
3. Items that unblock the most downstream work
4. Nice-to-have or low-risk items

## Error Handling

| Error Type | Action |
|------------|--------|
| PRD not validated | Warn and suggest running prd-validation first |
| PRD not enriched | Warn and suggest running prd-enrichment first |
| Requirement too vague to size | Flag as "needs-refinement" with questions |
| Circular dependency detected | Error with cycle description |

## Dependencies

- **prd-enrichment**: For technical context and complexity data
- **prd-validation**: For validated requirements and acceptance criteria
