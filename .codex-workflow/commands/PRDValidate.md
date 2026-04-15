---
description: Check PRD completeness and generate a validation report with scores and recommended actions
---

# PRDValidate

Validate a PRD for completeness, clarity, and implementability.

## Purpose

Performs structured validation of an imported PRD against quality criteria. Generates a scored report covering problem definition, user targeting, acceptance criteria, and security considerations. Surfaces gaps and generates clarifying questions for the Product Owner when the PRD does not meet the validation threshold.

## Arguments

- `prd`: PRD identifier - feature slug or Linear Initiative ID (required)
- `checklist`: Validation focus area (optional, default: `all`)
  - `completeness` - All required sections present and populated
  - `clarity` - Unambiguous language, clear definitions
  - `implementability` - Sufficient detail for engineering execution
  - `all` - Run all validation checks

## Execution

1. Load PRD from `/docs/planning/prds/{prd}.md`
   - Resolve identifier: match by feature slug or Initiative ID in frontmatter
   - Verify PRD status is `imported` or later
   - Load full document content and metadata

2. Invoke `prd-validation.validate` task
   - Run selected checklist(s) against PRD content
   - Score each validation dimension independently
   - Aggregate into overall completeness score

3. Generate validation report with scores
   - Problem Statement: clarity and specificity score (0-100)
   - User Definition: target user identification score (0-100)
   - Acceptance Criteria: testability and coverage score (0-100)
   - Security Considerations: threat model presence score (0-100)
   - Overall Completeness: weighted aggregate percentage

4. If score < threshold (default: 70%), generate questions for PO
   - Identify specific gaps per section
   - Generate targeted clarifying questions
   - Prioritize questions by impact on implementability

5. Store report in `/docs/planning/prds/{prd}-validation.md`
   - Update PRD status to `validated` or `needs-revision`
   - Link validation report in PRD frontmatter

## Prerequisites

- PRD imported via `/PRDIntake` (exists in `/docs/planning/prds/`)
- `.claude/tasks/prd-validation.md` task definition available

## Output Files

- `/docs/planning/prds/{prd}-validation.md` (validation report)
- `/docs/planning/prds/{prd}.md` (updated status in frontmatter)

## Output Format

```
### Validation Report: {prd}
Checklist: all

### Scores
- Problem Statement:      85/100
- User Definition:        60/100
- Acceptance Criteria:    45/100
- Security Considerations: 70/100
- **Overall Completeness: 65%**

### Recommended Actions
1. Add measurable acceptance criteria for search functionality
2. Define edge cases for multi-tenant scenarios
3. Include data retention requirements

### Questions for Product Owner
1. What is the expected response time SLA for the search endpoint?
2. Should deleted users retain access to shared resources?
3. Is GDPR data export in scope for v1?
```

## Example

```
/PRDValidate prd="user-onboarding" checklist="all"
```

```
/PRDValidate prd="INI-42" checklist="implementability"
```

## Related

- `/PRDIntake` - Import PRD before validation
- `/PRDEnrich` - Enrich validated PRD with technical context
- `/PRDFeasibility` - Assess feasibility after validation
- `/ImplementFeature` - Implement once PRD is validated

## Tasks Invoked

- `prd-validation.validate`

## Agents Used

- `prd-validator` - PRD quality analysis and scoring
