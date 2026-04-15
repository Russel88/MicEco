---
description: Recommend implementation order for multiple initiatives based on dependencies, risk, and value
---

# PRDSequence

Recommend an optimal implementation sequence for multiple initiatives.

## Purpose

Analyzes a set of initiatives with completed feasibility data and produces a recommended implementation order. Builds a dependency graph, applies configurable priority weights (risk, value, dependencies, or balanced), and generates an ordered plan with rationale and alternative sequences for decision-making.

## Arguments

- `initiatives`: Comma-separated list of Initiative IDs or feature slugs (required, e.g. `"INI-12, INI-15, INI-22"`)
- `weights`: Priority weighting strategy (optional, default: `balanced`)
  - `risk` - Prioritize lowest-risk initiatives first (de-risk early)
  - `value` - Prioritize highest business value first
  - `dependencies` - Prioritize based on dependency resolution order
  - `balanced` - Equal weighting across all factors

## Execution

1. Load feasibility data for each initiative
   - Resolve each identifier from `/docs/planning/prds/`
   - Load feasibility reports (`{prd}-feasibility.md`) for all initiatives
   - Verify all initiatives have status `enriched` or later
   - Flag any initiatives missing feasibility data

2. Build dependency graph
   - Extract inter-initiative dependencies from enrichment and feasibility reports
   - Identify shared infrastructure prerequisites
   - Detect circular dependencies and flag for resolution
   - Map team resource contention across initiatives

3. Apply sequencing algorithm based on weights
   - Score each initiative on: risk (from feasibility), value (from PRD), dependency depth, resource fit
   - Apply selected weighting strategy to compute priority scores
   - Resolve ordering conflicts using dependency constraints
   - Generate primary recommended sequence

4. Generate recommended order with rationale
   - Produce ordered list with per-initiative justification
   - Include dependency graph visualization (ASCII)
   - Generate 1-2 alternative sequences with trade-off explanations
   - Store sequence report

## Prerequisites

- All listed initiatives have completed `/PRDFeasibility` (feasibility reports exist)
- At least 2 initiatives provided for meaningful sequencing
- `.claude/tasks/cycle-planning.md` task definition available

## Output Files

- `/docs/planning/sequence-report-YYYY-MM-DD.md` (sequence recommendation)
- `/docs/planning/session-state.json` (updated with sequence record)

## Output Format

```
### Sequence Report
Strategy: balanced
Initiatives: 4

### Recommended Order
1. **INI-12: User Authentication Upgrade** (Score: 92)
   - Rationale: Unblocks INI-15 and INI-22, lowest risk, foundational
2. **INI-15: Profile Management** (Score: 78)
   - Rationale: Depends on INI-12, high value, medium complexity
3. **INI-22: Search Enhancement** (Score: 71)
   - Rationale: Independent but benefits from INI-12 auth patterns
4. **INI-18: Reporting Dashboard** (Score: 65)
   - Rationale: No dependencies, lowest value-to-effort ratio

### Dependency Graph
INI-12 ──> INI-15
  │
  └──> INI-22
INI-18 (independent)

### Alternative Sequences
**Alt A: Value-first**
INI-15 > INI-12 > INI-22 > INI-18
Trade-off: Higher initial risk, requires auth spike before INI-15

**Alt B: Risk-minimized**
INI-12 > INI-18 > INI-22 > INI-15
Trade-off: Delays highest-value initiative, but all early work is low-risk
```

## Example

```
/PRDSequence initiatives="INI-12, INI-15, INI-22, INI-18" weights="balanced"
```

```
/PRDSequence initiatives="user-onboarding, checkout-redesign, search-v2" weights="risk"
```

## Related

- `/PRDFeasibility` - Feasibility data required for sequencing
- `/PRDEnrich` - Enrichment data used in dependency analysis
- `/PRDValidate` - Validation must precede enrichment
- `/PRDIntake` - Initial PRD import
- `/StartSession` - Start session before sequencing work

## Tasks Invoked

- `cycle-planning.sequence`

## Agents Used

- `cycle-planner` - Dependency analysis and sequence optimization
