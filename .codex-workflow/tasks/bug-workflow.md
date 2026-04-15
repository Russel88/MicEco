---
name: bug-workflow
description: Orchestrate bug investigation, root cause analysis, and fix delivery
---

# Bug Workflow Task

Structured workflow for bug investigation and resolution, from triage through fix delivery with proper verification.

## Operations

### `investigate`

Full bug investigation workflow.

**Inputs:**
- `description`: Bug description
- `linearIssue`: Task management issue ID (optional)
- `affectedFiles`: Known affected files (optional)
- `reproSteps`: Reproduction steps (optional)

**Steps:**
1. **Triage Phase:**
   - Categorize severity: Critical / High / Medium / Low
   - Assess impact: Users affected, systems impacted
   - Check for workarounds
   - If issue provided, fetch full details via `mcp-sync`

2. **Investigation Phase:**
   - Invoke `context-loader.loadForBug` to gather context
   - Analyze suspected code paths
   - Review recent changes to affected files (git blame/log)
   - Search for similar patterns in codebase
   - Document minimal reproduction steps

3. **Root Cause Analysis:**
   - Trace issue to source
   - Classify issue type:
     - Logic error
     - Race condition
     - Data corruption
     - Integration failure
     - Configuration issue
   - Assess if architectural problem
   - Document unhandled edge cases

4. **Approval Checkpoint:**
   - Present findings summary
   - Proposed fix approach
   - Risk assessment
   - Wait for user confirmation

5. **Fix Phase (after approval):**
   - Route to `code-writer` agent with:
     - Root cause analysis
     - Fix requirements
     - Affected files
   - Generate fix code
   - Include defensive improvements

6. **Verification Phase:**
   - Outline test cases for the fix
   - Check for regression risks
   - Verify edge cases covered

7. **Delivery Phase:**
   - Update task management issue via `mcp-sync`
   - Update session state
   - Generate fix summary

**Outputs:**
```json
{
  "bug": {
    "description": "API returns null on valid requests",
    "severity": "High",
    "linearIssue": "PROJ-211"
  },
  "investigation": {
    "rootCause": "Missing null check on response",
    "issueType": "logic_error",
    "affectedFiles": [...],
    "reproSteps": [...]
  },
  "fix": {
    "approach": "Add validation and error handling",
    "filesModified": [...],
    "linesChanged": 85
  },
  "verification": {
    "testCases": [...],
    "regressionRisk": "Low"
  },
  "status": "fixed|investigating|blocked"
}
```

### `triage`

Quick triage without full investigation.

**Inputs:**
- `description`: Bug description
- `linearIssue`: Issue ID (optional)

**Steps:**
1. Categorize severity
2. Assess impact
3. Recommend priority
4. Return triage summary

**Outputs:**
```json
{
  "severity": "High",
  "impact": "Affects all users of feature",
  "recommendedPriority": "P1",
  "suggestedAssignee": "self",
  "workaround": "Refresh page to retry"
}
```

### `analyzeOnly`

Root cause analysis without fix generation.

**Inputs:**
- `description`: Bug description
- `affectedFiles`: Files to analyze

**Steps:**
1. Run investigation phase
2. Run root cause analysis
3. Return analysis without fix

**Outputs:**
```json
{
  "rootCause": "...",
  "issueType": "...",
  "suggestedFix": "...",
  "estimatedComplexity": "Low|Medium|High"
}
```

## Severity Classification

| Severity | Criteria | Response |
|----------|----------|----------|
| Critical | Production down, data loss, security breach | Immediate fix |
| High | Major feature broken, significant user impact | Same day |
| Medium | Feature degraded, workaround exists | This sprint |
| Low | Minor issue, cosmetic, edge case | Backlog |

## Issue Type Classification

| Type | Description | Typical Fix |
|------|-------------|-------------|
| `logic_error` | Incorrect business logic | Code change |
| `race_condition` | Timing/concurrency issue | Synchronization |
| `data_corruption` | Invalid data state | Validation + migration |
| `integration_failure` | External service issue | Error handling |
| `configuration_issue` | Wrong settings/env | Config update |
| `type_error` | TypeScript/type mismatch | Type fix |
| `ui_bug` | Visual/interaction issue | Component fix |

## Workflow Diagram

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│   Triage    │────▶│ Investigate │────▶│  Root Cause │
└─────────────┘     └─────────────┘     └─────────────┘
                                              │
                    ┌─────────────────────────┘
                    ▼
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  Approval   │────▶│    Fix      │────▶│   Verify    │
│ Checkpoint  │     │ Generation  │     │  & Deliver  │
└─────────────┘     └─────────────┘     └─────────────┘
```

## Agent Routing

| Phase | Agent | Purpose |
|-------|-------|---------|
| Investigation | None (task handles) | Code analysis |
| Root Cause | None (task handles) | Analysis |
| Fix | `code-writer` | Generate fix code |
| Verification | `test-planner` | Test case generation |

## Dependencies

- **context-loader**: For bug context
- **mcp-sync**: For task management updates
- **session-management**: For state updates

## Agents Used

- **code-writer**: Fix generation

## Approval Checkpoints

Workflow pauses for approval:
1. After root cause analysis (before fix)
2. Before applying fix to critical files

## Error Handling

| Error | Action |
|-------|--------|
| Cannot reproduce | Document attempts, request more info |
| Root cause unclear | Present multiple hypotheses |
| Fix introduces regression | Rollback, revise approach |
| Blocked by external | Mark as blocked, document dependency |

## Git Integration

Bug fixes follow this branch convention:
- Branch: `fix/{issue-id}-{short-description}`
- Commit: `fix({scope}): {description}`
- PR title: `Fix: {description} (#{issue-id})`
