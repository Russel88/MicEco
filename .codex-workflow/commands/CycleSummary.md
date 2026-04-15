---
description: Generate end-of-cycle summary with metrics, demo script, and release notes
---

# CycleSummary

Generate a comprehensive end-of-cycle summary for stakeholder review, including metrics, demo script, and release notes draft.

## Purpose

Aggregates all Cycle data into a summary document suitable for cycle review meetings, stakeholder updates, and release coordination. Optionally generates a demo script and release notes.

## Arguments

- `cycle`: Cycle identifier (required)
- `includeDemo`: Generate demo script (default: true)

## Execution

1. Load Cycle data from Linear
   - Fetch all Work Items and their final states
   - Collect milestone completion data
   - Gather team contribution metrics

2. Aggregate all session summaries from Cycle
   - Collect `/EndSession` summaries across the cycle period
   - Extract key accomplishments, decisions, and learnings

3. Calculate final metrics
   - Completion rate (planned vs. delivered)
   - Velocity (points/items per week)
   - Estimation accuracy
   - Blocker frequency and resolution time

4. Generate demo script if requested
   - Create ordered walkthrough of completed features
   - Include setup steps and key scenarios
   - Note any known limitations

5. Create stakeholder summary
   - High-level narrative of what was delivered
   - Business impact assessment
   - Carry-over items for next cycle

## Prerequisites

- Cycle exists in Linear (completed or in progress)
- Linear MCP configured and authenticated

## Output Format

```
### Cycle Summary: Complete user profile feature
**Duration:** Jan 20 - Feb 2 | **Status:** Completed

#### Delivery
- **Planned:** 12 items | **Delivered:** 11 items (92%)
- **Carry-over:** 1 item (WI-009, deferred due to external dependency)

#### Metrics
| Metric | Value | vs. Last Cycle |
|--------|-------|----------------|
| Velocity | 5.5 items/week | +10% |
| Estimation Accuracy | 85% | +5% |
| Blocker Resolution | 1.5 days avg | -0.5 days |

#### Key Accomplishments
1. User profile API complete with full CRUD operations
2. Profile UI with real-time updates
3. Data migration for existing users

#### Demo Script
1. Navigate to /profile
2. Edit profile fields (name, avatar, bio)
3. Show real-time sync across tabs
4. Demonstrate admin view

#### Release Notes Draft
- New: User profile management
- New: Avatar upload with image processing
- Improved: Account settings page layout
```

## Output Files

- `/docs/planning/cycles/cycle-{id}-summary.md`

## Example

```
/CycleSummary cycle="2025-01-20" includeDemo=true
```

## Related

- `/CycleStatus` - Mid-cycle status checks
- `/CycleRetro` - Retrospective analysis
- `/CycleCommit` - Original cycle commitment

## Tasks Invoked

- `cycle-monitoring.summary`

## Agents Used

- `status-aggregator`
- `doc-writer`
