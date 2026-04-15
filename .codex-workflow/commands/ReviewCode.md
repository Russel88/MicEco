---
description: Comprehensive code review for architecture, security, and PRD alignment
---

# ReviewCode

Comprehensive code review covering architecture, security, performance, and requirements alignment.

## Purpose

Performs thorough code review to ensure quality, security, and alignment with project standards before code is merged or deployed.

## Arguments

- `files`: Files to review (or `staged` for git staged files)
- `prUrl`: Pull request URL (optional)
- `focus`: Review focus areas (optional, multi-select)
  - `architecture` - Structure and patterns
  - `security` - Vulnerabilities
  - `performance` - Efficiency
  - `prd` - Requirements alignment

## Execution

1. Invoke `quality-gates` task with action: `review`
   - Gathers changed files
   - Routes to reviewer agent
   - Analyzes against standards
   - Categorizes findings
   - Generates report

2. Display:
   - Review summary
   - Critical/High/Medium/Low findings
   - PRD alignment status
   - Recommendations

## Prerequisites

- Files to review exist
- Relevant PRD/spec available (for prd focus)

## Output Format

```
### Review Summary
Files Reviewed: 5
Lines Analyzed: 450

### Findings
- **Critical (1)**: Security vulnerability in input handling
- **High (2)**: Missing error handling
- **Medium (5)**: Code style issues
- **Low (3)**: Optimization opportunities

### PRD Alignment
- Covered: REQ-001, REQ-002
- Missing: REQ-003 (partial implementation)
```

## Example

```
/ReviewCode files=["lib/api/handler.ts"] focus=["security", "prd"]
```

## Finding Categories

| Category | Description | Action Required |
|----------|-------------|-----------------|
| Critical | Security vulnerabilities, data loss risks | Must fix before merge |
| High | Significant issues, missing functionality | Should fix before merge |
| Medium | Code quality, maintainability | Consider fixing |
| Low | Style, optimization | Optional |

## Related

- `/ImplementFeature` - Implement features
- `/FixBug` - Fix identified issues
- `/GenerateTests` - Add test coverage
- `/RunSecurityAudit` - Deep security analysis

## Tasks Invoked

- `quality-gates.review`
- `context-loader.loadFull`

## Agents Used

- `reviewer` - Code analysis
