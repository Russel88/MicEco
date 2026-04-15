---
description: Generate detailed feasibility analysis with resource estimates and risk assessment for an enriched PRD
---

# PRDFeasibility

Generate a comprehensive technical feasibility analysis for an enriched PRD.

## Purpose

Produces a detailed feasibility assessment that evaluates whether a PRD can be implemented given team composition, timeline constraints, and technical complexity. Calculates resource requirements, assesses integration complexity, identifies prerequisites, and delivers a clear recommendation (PROCEED, DEFER, or REJECT) with supporting rationale.

## Arguments

- `prd`: PRD identifier - feature slug or Linear Initiative ID (required)
- `team`: Team composition string (required, e.g. `"2 BE, 1 FE"`, `"1 BE, 1 FE, 1 QA"`)
- `timeline`: Target timeline (optional, e.g. `"2 weeks"`, `"1 sprint"`, `"Q2 2025"`)

## Execution

1. Load enriched PRD
   - Resolve identifier from `/docs/planning/prds/{prd}.md`
   - Verify PRD status is `enriched`
   - Load enrichment report for complexity estimates and risk flags

2. Invoke `technical-feasibility.assess` task
   - Cross-reference complexity estimate with team composition
   - Model capacity against requirements volume
   - Factor in ramp-up time for unfamiliar technologies

3. Calculate resource requirements
   - Map enrichment complexity (T-shirt size) to person-day estimates
   - Distribute work across team roles (BE, FE, QA, DevOps)
   - Account for code review, testing, and deployment overhead
   - If `timeline` provided, check fit against available capacity

4. Assess integration complexity
   - Score each integration touchpoint from enrichment data
   - Evaluate API contract stability of dependencies
   - Assess data migration effort and rollback complexity
   - Factor in cross-team coordination overhead

5. Identify prerequisites
   - Infrastructure changes needed before implementation
   - Dependent features or services that must be completed first
   - Access or permissions required (API keys, environments)
   - Knowledge gaps requiring spike or research tasks

6. Generate feasibility report
   - Consolidate all assessments into structured report
   - Apply risk color coding (RED/AMBER/GREEN)
   - Produce final recommendation with conditions
   - Store report in `/docs/planning/prds/{prd}-feasibility.md`

## Prerequisites

- PRD enriched via `/PRDEnrich` (status: `enriched`)
- Team composition known
- `.claude/tasks/technical-feasibility.md` task definition available

## Output Files

- `/docs/planning/prds/{prd}-feasibility.md` (feasibility report)
- `/docs/planning/prds/{prd}.md` (updated status in frontmatter)

## Output Format

```
### Feasibility Report: {prd}
Team: 2 BE, 1 FE
Timeline: 2 weeks

### Technical Feasibility: HIGH

### Resource Estimate
- Backend:        8 person-days (2 engineers x 4 days)
- Frontend:       5 person-days (1 engineer x 5 days)
- Testing:        3 person-days
- DevOps:         1 person-day
- **Total:        17 person-days**
- Timeline Fit:   YES (10 working days x 3 engineers = 30 available)

### Integration Touchpoints
1. auth-service API (stable, low risk)
2. OAuth provider (external, medium risk)
3. notification-service (stable, low risk)

### Prerequisites
1. Redis cache expansion (DevOps, 1 day)
2. OAuth provider sandbox credentials (blocked on vendor)
3. Spike: mobile API pagination pattern (1 day)

### Risk Summary
- RED:   OAuth provider sandbox access not yet confirmed
- AMBER: Database migration requires maintenance window
- GREEN: All internal service APIs stable and documented

### Recommendation: PROCEED
Conditions:
- Confirm OAuth sandbox access before sprint start
- Schedule database migration window with SRE team
- Complete pagination pattern spike in first 2 days
```

## Example

```
/PRDFeasibility prd="user-onboarding" team="2 BE, 1 FE" timeline="2 weeks"
```

```
/PRDFeasibility prd="INI-42" team="1 BE, 1 FE, 1 QA"
```

## Related

- `/PRDEnrich` - Enrich PRD before feasibility assessment
- `/PRDValidate` - Validate PRD completeness
- `/PRDSequence` - Sequence initiatives using feasibility data
- `/ImplementFeature` - Implement after feasibility approval

## Tasks Invoked

- `technical-feasibility.assess`

## Agents Used

- `technical-analyst` - Resource modeling and risk assessment
