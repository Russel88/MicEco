---
description: Add technical context, risk analysis, and complexity estimates to a validated PRD
---

# PRDEnrich

Enrich a validated PRD with technical analysis, codebase patterns, risk flags, and complexity estimates.

## Purpose

Augments a validated PRD with engineering-relevant context by analyzing the existing codebase for relevant patterns, identifying technical risks and dependencies, and estimating implementation complexity. Produces actionable technical notes that bridge the gap between product requirements and engineering execution.

## Arguments

- `prd`: PRD identifier - feature slug or Linear Initiative ID (required)
- `depth`: Analysis depth (optional, default: `standard`)
  - `quick` - Surface-level pattern matching, key risks only
  - `standard` - Full codebase analysis, dependency mapping, complexity estimate
  - `deep` - Exhaustive analysis including cross-service impacts, performance modeling, migration planning

## Execution

1. Load validated PRD
   - Resolve identifier from `/docs/planning/prds/{prd}.md`
   - Verify PRD status is `validated`
   - Load validation report for context on any flagged areas

2. Invoke `prd-enrichment.analyze` task
   - Parse requirements into technical implications
   - Map each requirement to affected system components
   - Identify architectural patterns relevant to implementation

3. Analyze codebase for relevant patterns
   - Invoke `context-loader.codebasePatterns` to scan existing code
   - Identify similar implementations in the codebase
   - Surface reusable modules, services, and utilities
   - Flag deprecated patterns that should not be followed

4. Identify risks and dependencies
   - External service dependencies (APIs, databases, third-party)
   - Internal cross-team dependencies
   - Data migration requirements
   - Breaking change potential
   - Performance impact areas

5. Estimate complexity
   - Map requirements to effort categories
   - Generate T-shirt size estimate: XS, S, M, L, XL
   - Break down by component (backend, frontend, infrastructure)
   - Flag uncertainty areas that may affect estimate accuracy

6. Update PRD with enrichment section
   - Append technical context section to PRD document
   - Include architecture notes, pattern references, and risk flags
   - Update PRD status to `enriched` in frontmatter

7. Store enrichment report
   - Write detailed analysis to `/docs/planning/prds/{prd}-enrichment.md`
   - Link report in PRD frontmatter

## Prerequisites

- PRD validated via `/PRDValidate` (status: `validated`)
- Codebase accessible for pattern analysis
- `.claude/tasks/prd-enrichment.md` task definition available

## Output Files

- `/docs/planning/prds/{prd}-enrichment.md` (detailed enrichment report)
- `/docs/planning/prds/{prd}.md` (updated with technical context section)

## Output Format

```
### Enrichment Report: {prd}
Depth: standard

### Technical Notes
- Existing auth module can be extended for SSO requirement
- Search functionality aligns with existing Elasticsearch patterns
- No existing pagination pattern in mobile API - new pattern needed

### Risk Flags
- HIGH: Third-party OAuth provider has 99.5% SLA (below our 99.9% target)
- MEDIUM: Database migration required for new user fields
- LOW: Frontend bundle size increase ~15KB

### Dependency Map
- Internal: auth-service, user-service, notification-service
- External: OAuth provider, email gateway
- Infrastructure: Redis cache expansion needed

### Complexity Estimate
- Overall: L (Large)
- Backend: M | Frontend: L | Infrastructure: S
- Uncertainty: Medium (OAuth integration untested)
```

## Example

```
/PRDEnrich prd="user-onboarding" depth="standard"
```

```
/PRDEnrich prd="INI-42" depth="deep"
```

## Related

- `/PRDValidate` - Validate PRD before enrichment
- `/PRDFeasibility` - Assess feasibility after enrichment
- `/PRDSequence` - Sequence initiatives using enrichment data
- `/ImplementFeature` - Implement using enriched technical context

## Tasks Invoked

- `prd-enrichment.analyze`
- `context-loader.codebasePatterns`

## Agents Used

- `technical-analyst` - Codebase analysis and risk identification
