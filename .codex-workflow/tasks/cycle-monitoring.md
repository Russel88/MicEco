---
name: cycle-monitoring
description: Cycle progress tracking, status reporting, and retrospective analysis
---

# Cycle Monitoring Task

Tracks cycle progress in real time, generates status reports, produces end-of-cycle summaries, and drives retrospective analysis to improve future planning accuracy.

## Operations

### `status`

Generate current cycle status report.

**Steps:**
1. Fetch current cycle data from Linear
2. Aggregate session data from `/docs/planning/` directory
3. Calculate progress vs plan (items completed, hours burned)
4. Identify active blockers
5. Calculate velocity (actual vs planned)
6. Generate status report

**Inputs:**
- `format`: Report format (enum: `summary` | `detailed` | `blockers` | `risks`)

**Outputs:**
```json
{
  "cycle": "cycle-abc-123",
  "dayN": 5,
  "totalDays": 10,
  "progress": {
    "completed": 6,
    "inProgress": 3,
    "blocked": 1,
    "remaining": 4
  },
  "velocity": {
    "actual": 1.2,
    "planned": 1.4
  },
  "blockers": [
    { "item": "WB-005", "reason": "Waiting on SSO provider sandbox credentials", "daysSinceBlocked": 2 }
  ],
  "risks": [
    { "type": "velocity", "description": "Actual velocity 14% below plan, may not complete all items" }
  ]
}
```

### `updateProgress`

Record session completion and update cycle metrics.

**Steps:**
1. Load current cycle state
2. Update progress with session results
3. Recalculate velocity and projections
4. Store updated state

**Inputs:**
- `sessionId`: Session identifier (string)
- `itemsCompleted`: Array of completed item IDs (array of strings)
- `itemsBlocked`: Array of newly blocked item IDs (array of strings, optional)

**Outputs:**
```json
{
  "updated": true,
  "newProgress": {
    "completed": 8,
    "inProgress": 2,
    "blocked": 1,
    "remaining": 3
  }
}
```

### `summary`

Generate end-of-cycle summary.

**Steps:**
1. Load all cycle data
2. Aggregate all session summaries
3. Calculate final metrics (planned vs actual, completion rate)
4. Generate demo script listing completed features
5. Create stakeholder-friendly summary
6. Draft release notes

**Inputs:**
- `cycleId`: Cycle identifier (string)
- `includeDemo`: Include demo script in output (boolean, default: true)

**Outputs:**
```json
{
  "metrics": {
    "planned": 14,
    "completed": 12,
    "deferred": 2,
    "completionRate": "86%",
    "estimationAccuracy": "91%"
  },
  "demoScript": "1. Show new login flow with SSO button\n2. Demonstrate session persistence across tabs\n3. Show admin user management panel",
  "summary": "Completed 12 of 14 planned items (86%). Auth revamp core is done. SSO integration and session migration deferred to next cycle due to provider API delays.",
  "releaseNotes": "## What's New\n- Redesigned login flow\n- Redis-backed session storage\n- Password strength requirements\n\n## Known Issues\n- SSO integration pending (next cycle)",
  "summaryPath": "/docs/planning/cycle-summary-2026-01-31.md"
}
```

### `retrospective`

Analyze cycle for retrospective insights.

**Steps:**
1. Load cycle data and compare to plan
2. Calculate estimation accuracy per item
3. Identify blocker patterns (types, frequency, resolution time)
4. Calculate velocity trend vs previous cycles
5. Surface themes from session notes
6. Generate improvement suggestions

**Inputs:**
- `cycleId`: Cycle identifier (string)
- `compareCycleId`: Previous cycle to compare against (string, optional)

**Outputs:**
```json
{
  "estimationAccuracy": 85,
  "blockerPatterns": [
    { "type": "external-dependency", "frequency": 3, "avgResolutionDays": 2.5 },
    { "type": "unclear-requirements", "frequency": 1, "avgResolutionDays": 1 }
  ],
  "velocityTrend": "+15%",
  "themes": [
    "External dependencies caused most delays",
    "Testing estimates consistently undershot by 30%",
    "Frontend items completed faster than estimated"
  ],
  "suggestions": [
    "Add 30% buffer to testing estimates",
    "Request external API credentials 1 week before cycle start",
    "Consider pairing on backend items to reduce key-person risk"
  ]
}
```

## Configuration

Status report formats:
- **summary**: High-level progress, velocity, and key risks (default)
- **detailed**: Full item-by-item breakdown with time tracking
- **blockers**: Focused view on active blockers and resolution actions
- **risks**: Risk assessment with probability and impact ratings

## Error Handling

| Error Type | Action |
|------------|--------|
| No active cycle | Return error suggesting cycle-planning.commit first |
| Linear data unavailable | Use local session data only, mark as partial |
| Session data gaps | Interpolate from available data, flag gaps |
| Previous cycle not found for comparison | Skip comparison, return current cycle data only |

## Dependencies

- **mcp-sync**: For fetching cycle and issue data from Linear
- **session-management**: For session summaries and progress data
