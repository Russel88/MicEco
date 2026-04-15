---
description: Convert a validated PRD into sized, prioritized Work Items with dependency mapping
---

# Breakdown

Convert a PRD into actionable Work Items with size estimates, priorities, and dependency relationships.

## Purpose

Takes an enriched, validated PRD and generates a complete work breakdown structure with sized tickets, dependency graphs, and optional automatic creation in Linear.

## Arguments

- `prd`: PRD identifier (required)
- `targetSize`: Target ticket size (optional)
  - `2-4h` - Small tasks
  - `4-8h` - Medium tasks (default)
  - `1-2d` - Large tasks
- `includeTests`: Include test tickets (default: true)
- `dryRun`: Preview without creating tickets (default: false)

## Execution

1. Load enriched, validated PRD
   - Resolve PRD by identifier from `/docs/planning/prds/`
   - Validate PRD has passed enrichment and validation steps

2. Invoke `work-breakdown` task with action: `generate`
   - Decompose requirements into Work Items at target size
   - Assign area labels (Backend, Frontend, Infrastructure, etc.)
   - Estimate size ranges for each item
   - Set priority based on dependency order and business value

3. Generate Work Items from requirements
   - Map each requirement to one or more Work Items
   - Include acceptance criteria from PRD
   - Add test tickets if `includeTests` is enabled

4. Create dependency relationships
   - Identify blocking/blocked-by relationships
   - Generate dependency graph
   - Validate no circular dependencies

5. If not `dryRun`, invoke `mcp-sync` task with action: `createWorkItems`
   - Create issues in Linear
   - Set dependencies between issues
   - Apply labels and priorities

## Prerequisites

- PRD exists and has been validated
- Linear MCP configured (if not dryRun)

## Output Format

```
### Work Breakdown: {prd}
**Items:** 14 | **Test Tickets:** 6 | **Target Size:** 4-8h

#### Backend (5 items)
| ID | Title | Size | Priority |
|----|-------|------|----------|
| WI-001 | Set up database schema | 4-6h | P1 |
| WI-002 | Implement API endpoints | 6-8h | P1 |
| WI-003 | Add validation layer | 4-6h | P2 |

#### Frontend (4 items)
| ID | Title | Size | Priority |
|----|-------|------|----------|
| WI-006 | Create form components | 4-8h | P2 |
| WI-007 | Implement state management | 4-6h | P2 |

#### Dependency Graph
WI-001 ──> WI-002 ──> WI-003
                  └──> WI-006 ──> WI-007

**Total Estimate:** 68-102h (8.5-12.75 days)
```

## Output Files

- `/docs/planning/prds/{prd}-breakdown.md`

## Example

```
/Breakdown prd="user-profile-v2" targetSize="4-8h" includeTests=true dryRun=true
```

## Related

- `/CyclePlan` - Plan a cycle using breakdown output
- `/ImplementFeature` - Implement individual Work Items
- `/SyncLinear` - Manual Linear operations

## Tasks Invoked

- `work-breakdown.generate`
- `mcp-sync.createWorkItems`

## Agents Used

- `cycle-planner`
