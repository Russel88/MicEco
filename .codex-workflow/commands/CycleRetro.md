---
description: Analyze Cycle data for retrospective with estimation accuracy and pattern detection
---

# CycleRetro

Analyze Cycle data to generate retrospective insights including estimation accuracy, blocker patterns, velocity trends, and session data themes.

## Purpose

Provides data-driven retrospective analysis by examining Cycle metrics, comparing against previous cycles, and surfacing patterns to drive continuous improvement.

## Arguments

- `cycle`: Cycle identifier (required)
- `compare`: Previous cycle identifier to compare against (optional)

## Execution

1. Load Cycle data
   - Fetch completed Cycle from Linear
   - Load all Work Items with actual vs. estimated sizes
   - Collect session summaries from the cycle period

2. Analyze estimation accuracy
   - Compare estimated sizes to actual completion times
   - Identify consistently over/under-estimated areas
   - Calculate accuracy percentage and trend

3. Identify blocker patterns
   - Categorize blockers by type (external, technical, process)
   - Calculate average resolution time
   - Identify recurring blocker sources

4. Calculate velocity trends
   - Compare velocity to previous cycles
   - Identify acceleration or deceleration patterns
   - Normalize for team size and capacity changes

5. Surface session data themes
   - Analyze `/EndSession` summaries for recurring topics
   - Identify estimation patterns (which types of work are misjudged)
   - Detect requirement change frequency

## Prerequisites

- Cycle exists in Linear (completed)
- Linear MCP configured and authenticated
- Previous cycle data available (if using `compare`)

## Output Format

```
### Cycle Retrospective: Complete user profile feature
**Cycle:** Jan 20 - Feb 2

#### Estimation Accuracy
- **Overall:** 82% | **Trend:** +3% vs. last cycle
- **Over-estimated:** Infrastructure tasks (avg +30%)
- **Under-estimated:** UI integration tasks (avg -25%)

#### Blockers
- **Total:** 4 blockers
  - External: 2 (auth service, design assets)
  - Technical: 1 (database migration issue)
  - Process: 1 (unclear requirements)
- **Avg Resolution:** 1.8 days

#### Velocity Trend
- **This Cycle:** 5.5 items/week
- **Last Cycle:** 5.0 items/week
- **Change:** +10%

#### Patterns Identified
1. **Estimation:** UI integration consistently underestimated by 25%
   - Recommendation: Add 25% buffer to UI integration estimates
2. **Sessions:** 60% of sessions reported context-switching overhead
   - Recommendation: Batch similar work items in daily planning
3. **Requirements:** 3 items had scope changes mid-cycle
   - Recommendation: Add requirement freeze after day 2
```

## Output Files

- `/docs/planning/cycles/cycle-{id}-retro.md`

## Example

```
/CycleRetro cycle="2025-01-20" compare="2025-01-06"
```

## Related

- `/CycleSummary` - End-of-cycle delivery summary
- `/CycleStatus` - Mid-cycle status checks
- `/CyclePlan` - Plan next cycle with retro insights

## Tasks Invoked

- `cycle-monitoring.retrospective`
