---
name: prd-enrichment
description: Technical enrichment and complexity analysis for validated PRDs
---

# PRD Enrichment Task

Performs technical enrichment of validated PRDs by searching the codebase for relevant patterns, identifying reusable components, and generating complexity estimates to inform planning.

## Operations

### `analyze`

Search codebase for patterns relevant to PRD and add technical context.

**Steps:**
1. Load validated PRD
2. Extract key technical concepts from requirements
3. Search codebase for related patterns (existing implementations, similar features, shared utilities)
4. Identify reusable components and libraries
5. Flag potential conflicts with existing code
6. Map internal and external dependencies
7. Document integration points
8. Update PRD with enrichment section

**Inputs:**
- `prd`: PRD slug or ID (string)
- `depth`: Analysis depth (enum: `quick` | `standard` | `deep`)

**Outputs:**
```json
{
  "prd": "user-auth-revamp",
  "patterns": [
    { "name": "OAuth2 flow", "location": "lib/auth/oauth.ts", "relevance": "high" },
    { "name": "Session middleware", "location": "middleware/session.ts", "relevance": "high" }
  ],
  "reusableComponents": [
    { "name": "AuthProvider", "path": "components/auth/AuthProvider.tsx", "notes": "Can extend for SSO" }
  ],
  "conflicts": [
    { "description": "Current session store uses in-memory cache, new requirements need Redis", "severity": "medium" }
  ],
  "dependencies": {
    "internal": ["lib/auth/", "lib/session/", "middleware/"],
    "external": ["next-auth", "redis", "jose"]
  },
  "integrationPoints": [
    { "system": "User service API", "type": "REST", "impact": "schema change needed" },
    { "system": "Analytics pipeline", "type": "event", "impact": "new events required" }
  ]
}
```

### `estimateComplexity`

Generate complexity estimate for a PRD.

**Steps:**
1. Analyze requirement count and scope
2. Count integration touchpoints
3. Assess novelty vs existing patterns
4. Factor in dependency risks
5. Calculate T-shirt size (XS, S, M, L, XL)
6. Identify key unknowns that affect estimate

**Inputs:**
- `prd`: PRD slug or ID (string)

**Outputs:**
```json
{
  "complexity": "M",
  "confidence": "medium",
  "factors": [
    { "name": "requirement-count", "value": 12, "impact": "medium" },
    { "name": "integration-points", "value": 3, "impact": "high" },
    { "name": "novelty", "value": "partial - extends existing auth", "impact": "low" }
  ],
  "unknowns": [
    "Redis cluster configuration for session store",
    "SSO provider API rate limits"
  ],
  "estimate": {
    "optimistic": "2 weeks",
    "likely": "3 weeks",
    "pessimistic": "5 weeks"
  }
}
```

## Configuration

Complexity sizing thresholds:
- **XS**: < 1 day, single file change
- **S**: 1-3 days, few files, no integrations
- **M**: 1-2 weeks, multiple areas, some integrations
- **L**: 2-4 weeks, cross-cutting concerns, multiple integrations
- **XL**: 4+ weeks, architectural changes, high uncertainty

## Error Handling

| Error Type | Action |
|------------|--------|
| PRD not validated | Warn and suggest running validation first |
| Codebase search timeout | Return partial results with warning |
| No matching patterns found | Return empty patterns, flag as novel work |
| Conflicting dependencies | Flag in conflicts array with severity |

## Dependencies

- **context-loader**: For loading codebase and project context
- **git-workflow**: For accessing repository history and file analysis
