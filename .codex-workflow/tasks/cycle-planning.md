---
name: cycle-planning
description: Cycle planning, sequencing, and capacity management
---

# Cycle Planning Task

Manages cycle planning from initiative sequencing through capacity allocation to Linear cycle creation. Produces capacity-aware cycle plans with milestones and risk analysis.

## Operations

### `sequence`

Build initiative dependency graph and recommend implementation order.

**Steps:**
1. Load feasibility data for each initiative
2. Build dependency graph
3. Apply sequencing weights (risk, value, dependencies, balanced)
4. Calculate critical path
5. Generate ordered list with rationale for each position
6. Identify parallelizable streams

**Inputs:**
- `initiatives`: Array of initiative IDs or names (array of strings)
- `weights`: Sequencing strategy (enum: `risk` | `value` | `dependencies` | `balanced`)

**Outputs:**
```json
{
  "sequence": [
    {
      "initiative": "user-auth-revamp",
      "position": 1,
      "rationale": "No dependencies, unblocks SSO integration and admin dashboard",
      "parallelWith": []
    },
    {
      "initiative": "notification-system",
      "position": 2,
      "rationale": "Independent of auth, can run in parallel with SSO",
      "parallelWith": ["sso-integration"]
    },
    {
      "initiative": "sso-integration",
      "position": 2,
      "rationale": "Depends on auth-revamp, parallelizable with notifications",
      "parallelWith": ["notification-system"]
    }
  ],
  "criticalPath": ["user-auth-revamp", "sso-integration", "admin-dashboard"],
  "alternativeSequences": [
    {
      "strategy": "risk-first",
      "sequence": ["sso-integration", "user-auth-revamp", "notification-system"],
      "rationale": "Addresses highest-risk item first to surface unknowns early"
    }
  ]
}
```

### `plan`

Create capacity-aware cycle plan.

**Steps:**
1. Calculate available capacity from team composition and duration
2. Apply capacity percentage (e.g., 80% to account for meetings/interrupts)
3. Load prioritized backlog from Linear
4. Match backlog items to available capacity
5. Define milestones (checkpoints within the cycle)
6. Identify risks (over-commitment, dependency bottlenecks, key-person risk)
7. Generate draft cycle plan

**Inputs:**
- `duration`: Cycle length (enum: `1w` | `2w` | `3w`)
- `team`: Team members and their allocation (array of objects)
- `capacity`: Percentage of time available for planned work (number, default: 80)
- `goal`: High-level goal for the cycle (string)

**Outputs:**
```json
{
  "duration": "2w",
  "team": [
    { "name": "Alice", "role": "backend", "allocation": "100%" },
    { "name": "Bob", "role": "frontend", "allocation": "80%" }
  ],
  "availableHours": 128,
  "scope": [
    { "item": "WB-001", "assignee": "Alice", "hours": 6 },
    { "item": "WB-002", "assignee": "Alice", "hours": 8 },
    { "item": "WB-004", "assignee": "Bob", "hours": 6 }
  ],
  "milestones": [
    { "name": "Auth backend complete", "targetDay": 5, "items": ["WB-001", "WB-002"] },
    { "name": "Integration testing", "targetDay": 8, "items": ["WB-005", "WB-006"] }
  ],
  "risks": [
    { "type": "over-commitment", "description": "Alice at 95% capacity, no buffer for unknowns" },
    { "type": "key-person", "description": "All backend items depend on Alice" }
  ],
  "draftPath": "/docs/planning/cycle-draft-2026-01-20.md"
}
```

### `commit`

Commit a planned cycle to Linear.

**Steps:**
1. Load draft cycle plan
2. Validate all items exist in Linear
3. Create Cycle in Linear via mcp-sync
4. Assign work items to cycle
5. Set milestones
6. Return confirmation with Linear cycle URL

**Inputs:**
- `cycleDraft`: Path to the draft cycle plan (string)
- `adjustments`: Last-minute adjustments to the plan (object, optional)

**Outputs:**
```json
{
  "linearCycleId": "cycle-abc-123",
  "url": "https://linear.app/team/cycle/cycle-abc-123",
  "itemsAssigned": 12,
  "milestonesSet": 3
}
```

## Configuration

Capacity defaults:
- **Effective capacity**: 80% of total hours (accounts for meetings, interrupts, code review)
- **Buffer**: 10-20% of planned capacity held as contingency
- **Max allocation per person**: 90% (prevents single-point-of-failure)

## Error Handling

| Error Type | Action |
|------------|--------|
| Feasibility data missing | Warn and suggest running technical-feasibility first |
| Over-capacity plan | Warn with overflow amount, suggest items to defer |
| Linear items not found | List missing items, skip or create stubs |
| Circular initiative dependencies | Error with cycle description |

## Dependencies

- **mcp-sync**: For Linear cycle creation and issue management
- **technical-feasibility**: For initiative feasibility data and resource estimates
