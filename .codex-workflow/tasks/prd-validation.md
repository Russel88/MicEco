---
name: prd-validation
description: PRD intake, validation, and completeness scoring
---

# PRD Validation Task

Handles PRD intake from multiple sources, validates structure and completeness, and generates quality scores with actionable feedback for product owners.

## Operations

### `intake`

Fetch PRD from URL or file, parse structure, store locally with metadata.

**Steps:**
1. Determine source type (Coda URL, Confluence URL, local file)
2. Fetch content via appropriate MCP or file read
3. Parse PRD structure - identify sections (problem statement, user personas, requirements, acceptance criteria, non-functional requirements, security considerations)
4. Generate feature slug from PRD title
5. Store at `/docs/planning/prds/{feature-slug}.md` with YAML frontmatter (source URL, fetched date, status: draft)
6. Return PRD metadata object

**Inputs:**
- `url`: Source URL for the PRD (string, optional)
- `file`: Local file path for the PRD (string, optional)
- `initiative`: Initiative name or ID to associate (string, optional)

**Outputs:**
```json
{
  "id": "prd-001",
  "slug": "user-auth-revamp",
  "title": "User Authentication Revamp",
  "source": {
    "type": "coda",
    "url": "https://coda.io/d/...",
    "fetchedAt": "2026-01-16T10:00:00Z"
  },
  "sections": [
    "problem-statement",
    "user-personas",
    "requirements",
    "acceptance-criteria",
    "non-functional-requirements",
    "security-considerations"
  ],
  "status": "draft"
}
```

### `validate`

Run completeness and quality checks on a PRD.

**Steps:**
1. Load PRD from `/docs/planning/prds/{prd}.md`
2. Check completeness:
   - Problem statement present and clear
   - User personas defined
   - Functional requirements listed
   - Acceptance criteria for each requirement
   - Non-functional requirements specified
   - Security considerations addressed
   - Out of scope defined
3. Check clarity:
   - Requirements are testable
   - No ambiguous language ("should be fast", "user-friendly")
   - Quantifiable metrics where appropriate
4. Check implementability:
   - Technical constraints mentioned
   - Integration points identified
   - Data model implications considered
5. Score each section 0-100
6. Calculate overall completeness score
7. Generate questions for gaps (things the PO needs to clarify)
8. Store validation report

**Inputs:**
- `prd`: PRD slug or ID (string)
- `checklist`: Validation scope (enum: `completeness` | `clarity` | `implementability` | `all`)

**Outputs:**
```json
{
  "prd": "user-auth-revamp",
  "score": 78,
  "sections": [
    {
      "name": "problem-statement",
      "score": 95,
      "issues": [],
      "questions": []
    },
    {
      "name": "acceptance-criteria",
      "score": 60,
      "issues": ["3 requirements missing acceptance criteria"],
      "questions": ["What is the expected latency for SSO login?"]
    }
  ],
  "overallScore": 78,
  "questionsForPO": [
    "What is the expected latency for SSO login?",
    "Are there specific compliance requirements for password storage?"
  ],
  "reportPath": "/docs/planning/prds/user-auth-revamp-validation.md"
}
```

### `checklistConfig`

Load custom validation rules and override default thresholds.

**Steps:**
1. Load custom rules from `/docs/planning/validation-rules.md` if exists
2. Merge with defaults
3. Return active ruleset

**Inputs:**
- `customRulesPath`: Path to custom validation rules file (string, optional)

**Outputs:**
```json
{
  "rules": [
    { "id": "req-001", "name": "problem-statement-required", "severity": "error" },
    { "id": "req-002", "name": "acceptance-criteria-required", "severity": "error" },
    { "id": "req-003", "name": "security-section-present", "severity": "warning" }
  ],
  "thresholds": {
    "pass": 80,
    "warn": 60
  }
}
```

## Configuration

Validation thresholds can be customized per project via `/docs/planning/validation-rules.md`. Default thresholds:
- **Pass**: 80+ overall score
- **Warn**: 60-79 overall score
- **Fail**: Below 60

## Error Handling

| Error Type | Action |
|------------|--------|
| PRD not found | Return error with available PRD list |
| Source URL unreachable | Retry once, then fail with network error |
| Parse failure | Return partial parse with warnings |
| Missing sections | Score as 0 for missing, include in questions |

## Dependencies

- **mcp-sync**: For fetching PRDs from Coda/Confluence documentation sources
