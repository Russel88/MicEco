---
description: Plan a capacity-aware Cycle with scope recommendations and risk assessment
---

# CyclePlan

Plan a capacity-aware development Cycle with scope recommendations, milestones, and risk assessment.

## Purpose

Generates a draft Cycle plan by loading team capacity and prioritized backlog, then recommending scope that fits within available capacity while identifying milestones and risks.

## Arguments

- `duration`: Cycle duration (required)
  - `1w` - One week
  - `2w` - Two weeks
  - `3w` - Three weeks
- `team`: Team members, comma-separated (required)
- `capacity`: Capacity percentage (default: 80%)
- `goal`: Cycle goal statement (required)

## Execution

1. Load team capacity
   - Calculate available hours per team member for the duration
   - Apply capacity percentage to account for meetings, support, etc.

2. Load prioritized backlog from Linear
   - Invoke `mcp-sync` task with action: `getBacklog`
   - Fetch Work Items sorted by priority
   - Include size estimates and dependencies

3. Invoke `cycle-planning` task with action: `plan`
   - Match backlog items to available capacity
   - Respect dependency ordering
   - Balance work across team members

4. Recommend scope based on capacity
   - Flag items that fit within capacity
   - Identify stretch goals
   - Highlight items deferred to next cycle

5. Identify milestones
   - Group Work Items into logical milestones
   - Set milestone target dates within the cycle

6. Generate risk assessment
   - Identify capacity risks (overallocation, single points of failure)
   - Flag dependency risks (external blockers, long chains)
   - Note estimation uncertainty

7. Save draft Cycle plan
   - Write plan document for review before committing

## Prerequisites

- Work Items exist in Linear backlog (from `/Breakdown`)
- Team members defined
- Linear MCP configured and authenticated

## Output Format

```
### Cycle Plan: {goal}
**Duration:** 2w (Jan 20 - Feb 2) | **Capacity:** 80%

#### Team Capacity
| Member | Available | Allocated | Remaining |
|--------|-----------|-----------|-----------|
| Alice  | 48h       | 42h       | 6h        |
| Bob    | 48h       | 44h       | 4h        |

#### Recommended Scope (12 items)
| ID | Title | Assignee | Size | Priority |
|----|-------|----------|------|----------|
| WI-001 | Database schema | Alice | 4-6h | P1 |
| WI-002 | API endpoints | Bob | 6-8h | P1 |

#### Milestones
1. **API Complete** - Jan 27 (WI-001, WI-002, WI-003)
2. **UI Complete** - Jan 31 (WI-006, WI-007)

#### Risk Register
| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| Auth service dependency | High | Medium | Early spike |
```

## Output Files

- `/docs/planning/cycles/cycle-YYYY-MM-DD.md`

## Example

```
/CyclePlan duration="2w" team="Alice,Bob,Carol" capacity="80%" goal="Complete user profile feature"
```

## Related

- `/Breakdown` - Generate Work Items for planning
- `/CycleCommit` - Commit this plan to Linear
- `/CycleStatus` - Monitor cycle progress

## Tasks Invoked

- `cycle-planning.plan`
- `mcp-sync.getBacklog`

## Agents Used

- `cycle-planner`
