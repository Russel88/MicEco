---
name: priority-matrix
description: Build prioritized task list from issues and project context
---

# Priority Matrix Task

Builds a prioritized task list by combining issues, blockers, dependencies, and business impact scoring into an actionable work queue.

## Operations

### `build`

Build complete priority matrix.

**Inputs:**
- `team`: Team filter (optional)
- `assignee`: Filter by assignee (default: "me")
- `includeBacklog`: Include backlog items (default: false)

**Steps:**
1. Fetch issues via `mcp-sync.syncLinear.fetchIssues`
2. Load blockers from `session-state.json`
3. Analyze dependencies between issues
4. Score each issue:
   - Base priority (P0-P4)
   - Blocker multiplier (+2 if blocking others)
   - Dependency penalty (-1 if blocked)
   - Age factor (+0.5 per week old)
   - Business impact (from labels)
5. Sort by composite score
6. Group into priority tiers
7. Return matrix

**Outputs:**
```json
{
  "matrix": {
    "P0_Critical": [
      {
        "id": "PROJ-211",
        "title": "Critical bug fix",
        "score": 9.5,
        "factors": {
          "basePriority": 1,
          "blocking": ["PROJ-212"],
          "blockedBy": [],
          "ageWeeks": 1,
          "labels": ["bug"]
        },
        "recommendedAction": "Fix immediately - blocking other work"
      }
    ],
    "P1_High": [...],
    "P2_Medium": [...],
    "P3_Low": [],
    "Backlog": []
  },
  "summary": {
    "totalItems": 15,
    "blockedItems": 3,
    "blockingItems": 2,
    "overdueItems": 1
  },
  "generatedAt": "2026-01-16T10:00:00Z"
}
```

### `score`

Score a single issue.

**Inputs:**
- `issue`: Issue object

**Steps:**
1. Apply scoring algorithm:
   ```
   score = basePriority * 2
         + (isBlocking ? 2 : 0)
         - (isBlocked ? 1 : 0)
         + (ageWeeks * 0.5)
         + (hasBugLabel ? 1 : 0)
         + (hasSecurityLabel ? 3 : 0)
   ```
2. Return score with breakdown

**Outputs:**
```json
{
  "issueId": "PROJ-211",
  "score": 9.5,
  "breakdown": {
    "basePriority": 2,
    "blockingBonus": 2,
    "blockedPenalty": 0,
    "ageBonus": 0.5,
    "labelBonus": 1
  }
}
```

### `recommendNext`

Recommend the next task to work on.

**Inputs:**
- `currentWork`: Currently in-progress items (from session state)
- `preferredType`: Optional preference (`bug` | `feature` | `refactor`)

**Steps:**
1. Build priority matrix
2. Filter out in-progress items
3. Filter out blocked items
4. Apply preference filter if specified
5. Return top recommendation with rationale

**Outputs:**
```json
{
  "recommended": {
    "id": "PROJ-211",
    "title": "Critical bug fix",
    "rationale": "Highest priority bug, blocking 2 other tasks"
  },
  "alternatives": [
    {
      "id": "PROJ-213",
      "rationale": "Next highest priority, no dependencies"
    }
  ]
}
```

### `detectBlockers`

Identify blocking relationships and circular dependencies.

**Inputs:**
- `issues`: Array of issues to analyze

**Steps:**
1. Build dependency graph from issue links
2. Detect cycles (circular blocking)
3. Identify critical path (longest chain)
4. Flag unresolved blockers
5. Return analysis

**Outputs:**
```json
{
  "blockingChains": [
    {
      "chain": ["PROJ-211", "PROJ-212", "PROJ-213"],
      "length": 3,
      "criticalPath": true
    }
  ],
  "circularDependencies": [],
  "unresolvedBlockers": [
    {
      "issue": "PROJ-215",
      "blockedBy": "External: Waiting on API key"
    }
  ]
}
```

## Priority Tiers

| Tier | Score Range | Typical Items |
|------|-------------|---------------|
| P0_Critical | 8+ | Security issues, production bugs, blockers |
| P1_High | 5-7.9 | Active sprint items, high-value features |
| P2_Medium | 3-4.9 | Normal priority, planned work |
| P3_Low | 1-2.9 | Nice-to-have, tech debt |
| Backlog | <1 | Future consideration |

## Scoring Weights

| Factor | Weight | Condition |
|--------|--------|-----------|
| Base Priority | x2 | Priority level (1=Urgent, 4=Low) |
| Blocking | +2 | Issue blocks other issues |
| Blocked | -1 | Issue is blocked |
| Age | +0.5/week | Time since creation |
| Bug Label | +1 | Has "bug" label |
| Security Label | +3 | Has "security" label |
| Due Date | +2 | Due within 3 days |

## Dependencies

- **mcp-sync**: For issue fetching
- **session-management**: For current blocker state

## Integration with Session

The priority matrix is:
- Generated on `StartSession`
- Included in session plan
- Updated when issues change
- Referenced for `recommendNext`

## Error Handling

| Error | Action |
|-------|--------|
| Task management unavailable | Use cached matrix if <1h old, warn user |
| No issues found | Return empty matrix, suggest filters |
| Circular dependency | Flag in output, recommend manual resolution |
