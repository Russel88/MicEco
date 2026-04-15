---
name: requirement-discovery
description: Find and compile requirements from PRDs, documentation, and task management sources
---

# Requirement Discovery Task

Discovers and compiles requirements from multiple authoritative sources (local PRDs, external documentation, task management) into a unified requirements object for feature implementation.

## Operations

### `discover`

Find requirements for a feature or task.

**Inputs:**
- `query`: Feature name, keyword, or issue ID
- `sources`: Array of sources to search (default: all)
  - `local` - Search `/knowledge/prd/` and `/docs/specs/`
  - `external` - Query external documentation MCP
  - `taskManagement` - Query task management MCP

**Steps:**
1. **Local Search:**
   - Glob `/knowledge/prd/*.md` for matching files
   - Glob `/docs/specs/*.md` for matching specs
   - Extract relevant sections using keyword matching
2. **External Search (if enabled):**
   - Invoke `mcp-sync.syncCoda.fetchPRD` for known PRD pages
   - Search documentation for matching content
3. **Task Management Search (if enabled):**
   - Invoke `mcp-sync.syncLinear.fetchIssues` with query filter
   - Extract requirements from issue descriptions
   - Collect acceptance criteria from comments
4. **Merge and deduplicate:**
   - Combine requirements from all sources
   - Flag conflicts between sources
   - Establish source priority (Task Management > External > Local)
5. Return unified requirements object

**Outputs:**
```json
{
  "feature": "User Profile Integration",
  "sources": {
    "local": {
      "prd": "/knowledge/prd/feature-prd.md",
      "spec": "/docs/specs/feature-spec.md"
    },
    "external": {
      "docId": "doc123",
      "pageId": "PRD Page"
    },
    "taskManagement": {
      "issues": ["PROJ-211", "PROJ-212"]
    }
  },
  "requirements": [
    {
      "id": "REQ-001",
      "description": "Requirement description",
      "source": "taskManagement:PROJ-211",
      "priority": "P1",
      "acceptanceCriteria": [...]
    }
  ],
  "dependencies": [...],
  "conflicts": [],
  "discoveredAt": "2026-01-16T10:00:00Z"
}
```

### `extractAcceptanceCriteria`

Extract acceptance criteria from a specific source.

**Inputs:**
- `content`: Markdown content containing requirements
- `format`: Expected format (`checklist` | `prose` | `gherkin`)

**Steps:**
1. Parse content for acceptance criteria patterns:
   - Checkbox lists (`- [ ]`)
   - Numbered lists
   - "Given/When/Then" statements
   - "Should" statements
2. Normalize to standard format
3. Return criteria array

**Outputs:**
```json
{
  "criteria": [
    {
      "id": "AC-001",
      "description": "User can view profile",
      "testable": true
    }
  ],
  "format": "checklist",
  "extractedFrom": "local"
}
```

### `mapToImplementation`

Map requirements to implementation artifacts.

**Inputs:**
- `requirements`: Requirements object from `discover`
- `codebase`: Code location context from `context-loader`

**Steps:**
1. Analyze requirements for implementation hints:
   - UI components mentioned → map to component files
   - API endpoints mentioned → map to route files
   - Database fields mentioned → map to schema/types
2. Search codebase for existing implementations
3. Identify gaps (requirements without implementation)
4. Return mapping

**Outputs:**
```json
{
  "mappings": [
    {
      "requirement": "REQ-001",
      "implementations": [
        {
          "file": "lib/feature/handler.ts",
          "status": "partial",
          "coverage": 0.6
        }
      ],
      "gaps": ["Feature not fully implemented"]
    }
  ],
  "unmappedRequirements": [],
  "orphanedCode": []
}
```

## Source Priority

When requirements conflict between sources:

1. **Task Management** (highest) - Most current, actively managed
2. **External Documentation** - Strategic documentation
3. **Local PRD** - Baseline requirements
4. **Local Spec** (lowest) - Implementation details

Conflicts are flagged for human resolution.

## Pattern Recognition

The task recognizes these requirement patterns:

| Pattern | Example | Extraction |
|---------|---------|------------|
| User story | "As a user, I want to..." | Full story |
| Acceptance criteria | "- [ ] User can..." | Checkbox item |
| Gherkin | "Given... When... Then..." | Full scenario |
| Should statement | "The system should..." | Full statement |
| Must/Shall | "Must validate..." | Full requirement |

## Dependencies

- **mcp-sync**: For external and task management queries
- **context-loader**: For codebase context

## Error Handling

| Error | Action |
|-------|--------|
| Source unavailable | Skip source, continue with others, flag in output |
| No requirements found | Return empty with search terms attempted |
| Conflicting requirements | Include both, flag conflict for resolution |
| Malformed content | Best-effort extraction, log warning |

## Caching

- Local file results cached for session duration
- External/task management results cached for 5 minutes
- Manual refresh available via `force: true` parameter
