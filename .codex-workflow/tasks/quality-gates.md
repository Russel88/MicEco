---
name: quality-gates
description: Orchestrate code review, testing, and security audit workflows
---

# Quality Gates Task

Orchestrates quality assurance workflows including code review, test planning, and security audits. Ensures code meets quality standards before delivery.

## Operations

### `review`

Comprehensive code review.

**Inputs:**
- `files`: Files to review (or "staged" for git staged files)
- `prUrl`: Pull request URL (optional)
- `focus`: Review focus areas (optional)
  - `architecture` - Structure and patterns
  - `security` - Vulnerabilities
  - `performance` - Efficiency
  - `prd` - Requirements alignment

**Steps:**
1. Gather context:
   - Load changed files
   - If PR, fetch PR details
   - Load relevant specs from `context-loader`
2. Route to `reviewer` agent with:
   - File contents
   - Architecture rules from CLAUDE.md
   - Relevant PRD sections
3. Analyze for:
   - Architecture alignment
   - Security vulnerabilities (OWASP top 10)
   - Performance issues
   - Code style compliance
   - PRD requirement coverage
4. Categorize findings:
   - Critical (must fix)
   - High (should fix)
   - Medium (consider)
   - Low (optional)
5. Generate review report

**Outputs:**
```json
{
  "review": {
    "filesReviewed": 5,
    "linesAnalyzed": 450
  },
  "findings": {
    "critical": [...],
    "high": [...],
    "medium": [...],
    "low": [...]
  },
  "summary": {
    "critical": 1,
    "high": 2,
    "medium": 5,
    "low": 3,
    "passed": false
  },
  "prdAlignment": {
    "covered": ["REQ-001", "REQ-002"],
    "missing": ["REQ-003"]
  }
}
```

### `generateTests`

Generate test plan and specifications.

**Inputs:**
- `target`: File, function, or feature to test
- `types`: Test types to generate
  - `unit` - Unit tests
  - `integration` - Integration tests
  - `e2e` - End-to-end tests
- `coverage`: Target coverage percentage (default: 80)

**Steps:**
1. Load target context via `context-loader`
2. Route to `test-planner` agent with:
   - Target code
   - Existing tests (if any)
   - Coverage requirements
3. Generate:
   - Test cases with descriptions
   - Edge cases to cover
   - Mock requirements
   - Expected behaviors
4. Return test plan (not code)

**Outputs:**
```json
{
  "target": "lib/api/handler.ts",
  "testPlan": {
    "unit": [...],
    "integration": [...],
    "edgeCases": [...]
  },
  "estimatedCoverage": 85
}
```

### `securityAudit`

Security-focused audit.

**Inputs:**
- `scope`: Audit scope
  - `full` - Entire codebase
  - `changed` - Recently changed files
  - `specific` - Specific files/directories
- `files`: Files to audit (for `specific` scope)

**Steps:**
1. Gather code for audit
2. Route to `researcher` agent for:
   - OWASP vulnerability patterns
   - Dependency vulnerability check
   - Configuration security
   - Authentication/authorization review
3. Analyze findings
4. Generate audit report
5. Optionally create issues for findings

**Outputs:**
```json
{
  "audit": {
    "scope": "changed",
    "filesAudited": 12,
    "timestamp": "2026-01-16T10:00:00Z"
  },
  "vulnerabilities": [...],
  "dependencies": {
    "vulnerable": 0,
    "outdated": 3
  },
  "configuration": {
    "issues": [...]
  },
  "summary": {
    "high": 1,
    "medium": 2,
    "low": 5,
    "passed": false
  }
}
```

### `preCommit`

Quick pre-commit quality check.

**Inputs:**
- None (uses git staged files)

**Steps:**
1. Get staged files
2. Run quick lint/type check
3. Run minimal security scan
4. Return pass/fail with issues

**Outputs:**
```json
{
  "passed": false,
  "checks": {
    "lint": { "passed": true },
    "types": { "passed": false, "errors": 2 },
    "security": { "passed": true }
  },
  "blockers": [...]
}
```

## Quality Thresholds

| Check | Threshold | Blocking |
|-------|-----------|----------|
| Critical findings | 0 | Yes |
| High findings | 0 | Yes (for main branch) |
| Test coverage | 80% | No (warning) |
| Type errors | 0 | Yes |
| Lint errors | 0 | Yes |

## Agent Routing

| Operation | Agent | Purpose |
|-----------|-------|---------|
| review | `reviewer` | Code analysis |
| generateTests | `test-planner` | Test specifications |
| securityAudit | `researcher` | Security patterns |

## Dependencies

- **context-loader**: For code context
- **mcp-sync**: For issue creation

## Agents Used

- **reviewer**: Code review analysis
- **test-planner**: Test case generation
- **researcher**: Security standards

## Report Output

Reports are written to `/docs/reports/`:
- Code review: `code-review-YYYY-MM-DD.md`
- Security audit: `security-audit-YYYY-MM-DD.md`
- Test plan: `test-plan-{target}-YYYY-MM-DD.md`

## Error Handling

| Error | Action |
|-------|--------|
| Files not found | Skip missing, report in output |
| Agent timeout | Return partial results, flag incomplete |
| Parse error | Best-effort analysis, log warning |
