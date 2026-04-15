---
name: workflow-cycle-status
description: Get current Cycle status with progress, blockers, and async check-in data Use when the user asks to cycle status in this repository or refers to the old `/CycleStatus` command.
---

# Cycle Status

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# CycleStatus

Get the current Cycle status including progress metrics, blockers, risks, and aggregated async check-in data.

## Purpose

Provides a comprehensive view of the current Cycle by aggregating Linear progress data with `workflow-end-session` check-in summaries, calculating progress against plan, and surfacing blockers and risks.

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

2. Aggregate `workflow-end-session` data from team
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
`workflow-cycle-status` format="summary" output="both"
```

## Related

- `workflow-cycle-plan` - Plan the cycle
- `workflow-cycle-commit` - Commit cycle scope
- `workflow-cycle-summary` - End-of-cycle summary
- `workflow-end-session` - Session check-in data source

## Tasks Invoked

- `cycle-monitoring.status`
- `mcp-sync.getCycleProgress`

## Agents Used

- `status-aggregator`
