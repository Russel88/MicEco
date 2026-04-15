---
description: Execute full security audit for vulnerabilities, dependencies, and configuration
---

# RunSecurityAudit

Comprehensive security audit covering code vulnerabilities, dependencies, and configuration.

## Purpose

Performs thorough security analysis to identify vulnerabilities, outdated dependencies, and configuration issues before deployment.

## Arguments

- `scope`: Audit scope (required)
  - `full` - Entire codebase
  - `changed` - Recently changed files only
  - `specific` - Specific files/directories
- `files`: Files to audit (for `specific` scope)
- `createIssues`: Create task management issues for findings (default: false)

## Execution

1. Invoke `quality-gates` task with action: `securityAudit`
   - Gathers code for audit
   - Routes to researcher agent
   - Checks OWASP patterns
   - Analyzes dependencies
   - Reviews configuration
   - Generates report

2. Display:
   - Audit summary
   - Vulnerabilities by severity
   - Dependency status
   - Configuration issues
   - Remediation recommendations

## Prerequisites

- Code exists for audit scope

## Output Files

- `/docs/reports/security-audit-YYYY-MM-DD.md`

## Output Format

```
### Security Audit Report
Scope: changed
Files Audited: 12
Timestamp: 2026-01-16T10:00:00Z

### Vulnerabilities
| Severity | Count | Top Issue |
|----------|-------|-----------|
| High | 1 | SQL injection in search |
| Medium | 2 | Missing input validation |
| Low | 5 | Verbose error messages |

### Dependencies
- Vulnerable: 0
- Outdated: 3 (non-critical)

### Configuration
- CSP headers: Missing frame-ancestors
- CORS: Properly configured

### Recommendations
1. [HIGH] Add parameterized queries to search endpoint
2. [MEDIUM] Validate user input in forms
...
```

## Example

```
/RunSecurityAudit scope="changed" createIssues=true
```

## Vulnerability Categories

| Category | Examples | Severity |
|----------|----------|----------|
| Injection | SQL, XSS, Command | High |
| Authentication | Weak tokens, session issues | High |
| Data Exposure | PII leaks, verbose errors | Medium |
| Configuration | Missing headers, CORS | Medium |
| Dependencies | Known CVEs | Varies |

## Related

- `/ReviewCode` - General code review
- `/FixBug` - Fix security issues

## Tasks Invoked

- `quality-gates.securityAudit`
- `mcp-sync.syncLinear` (if createIssues=true)

## Agents Used

- `researcher` - Security standards and patterns
