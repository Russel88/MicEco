---
description: Get current Cycle status with progress, blockers, and async check-in data
---

# CycleStatus

Get the current Cycle status including progress metrics, blockers, risks, and aggregated async check-in data.

## Purpose

Provides a comprehensive view of the current Cycle by aggregating Linear progress data with `/EndSession` check-in summaries, calculating progress against plan, and surfacing blockers and risks.

## Arguments

- `format`: Output format (optional)
  - `summary` - High-level overview (default)
  - `detailed` - Full item-by-item status
  - `blockers` - Blockers and risks only
  - `risks` - Risk assessment only
- `output`: Where to post the report (optional)
  - `terminal` - Display in terminal (default)
  - `slack` - Post to Slack channel
  - `both` - Terminal and Slack

## Execution

1. Load current Cycle from Linear
   - Fetch active Cycle and all associated Work Items
   - Get current state of each item

2. Aggregate `/EndSession` data from team
   - Collect session summaries from the current Cycle period
   - Extract progress updates, blockers, and notes

3. Calculate progress vs. plan
   - Compare completed items against milestone targets
   - Calculate velocity (actual vs. planned)
   - Determine on-track/at-risk/off-track status

4. Identify blockers and risks
   - Surface items marked as blocked
   - Identify items at risk of missing milestones
   - Flag capacity concerns

5. Generate status report
   - Format according to requested format

6. Post to Slack if configured
   - Send formatted report to team channel

## Prerequisites

- Active Cycle exists in Linear
- Linear MCP configured and authenticated
- Slack MCP configured (if posting to Slack)

## Output Format

```
### Cycle Status: Complete user profile feature
**Duration:** Jan 20 - Feb 2 | **Day:** 6 of 10
**Progress:** 58% | **Status:** On Track

#### Summary
| Status | Count |
|--------|-------|
| Completed | 7 |
| In Progress | 3 |
| Blocked | 1 |
| Not Started | 1 |

#### Blockers
1. **WI-009** - Auth service API not available (blocked 2 days)
   - Owner: Bob
   - Impact: Blocks WI-010, WI-011
   - Mitigation: Escalated to platform team

#### Velocity
- **Planned:** 6 items/week
- **Actual:** 5.8 items/week
- **Trend:** Stable
```

## Example

```
/CycleStatus format="summary" output="both"
```

## Related

- `/CyclePlan` - Plan the cycle
- `/CycleCommit` - Commit cycle scope
- `/CycleSummary` - End-of-cycle summary
- `/EndSession` - Session check-in data source

## Tasks Invoked

- `cycle-monitoring.status`
- `mcp-sync.getCycleProgress`

## Agents Used

- `status-aggregator`
