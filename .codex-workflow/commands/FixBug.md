---
description: Diagnose and fix bugs with structured triage, root cause analysis, and verification
---

# FixBug

Structured bug investigation and resolution workflow.

## Purpose

Guides the bug fixing process from triage through fix delivery, ensuring proper root cause analysis, approval checkpoints, and verification.

## Arguments

- `description`: Bug description (required)
- `linearIssue`: Task management issue ID (optional)
- `affectedFiles`: Known affected files (optional)
- `severity`: Initial severity assessment (optional: `critical` | `high` | `medium` | `low`)

## Execution

1. Invoke `bug-workflow` task with action: `investigate`
   - Triage: Categorize severity and impact
   - Investigation: Analyze code paths, recent changes
   - Root Cause: Identify source of issue
   - Approval: Present findings before fix
   - Fix: Generate fix code via code-writer agent
   - Verification: Outline test cases

2. Display:
   - Severity and impact assessment
   - Root cause analysis
   - Fix approach
   - Test verification plan

## Prerequisites

- Bug description or task management issue

## Example

```
/FixBug description="API returns null on valid requests" linearIssue="PROJ-211" affectedFiles=["lib/api/handler.ts"]
```

## Workflow Phases

### 1. Triage
- Categorize: Critical / High / Medium / Low
- Assess impact: Users affected, systems impacted
- Check for workarounds

### 2. Investigation
- Reproduce the issue
- Analyze suspected code paths
- Review recent changes (git blame)

### 3. Root Cause Analysis
- Trace issue to source
- Classify issue type (logic error, race condition, etc.)
- Assess architectural implications

### 4. Approval Checkpoint
**STOP** - Present findings for review before proceeding

### 5. Fix Implementation
- Generate targeted fix
- Include defensive programming
- Add comments explaining fix rationale

### 6. Verification
- Define regression tests
- Outline manual verification steps

## Related

- `/ImplementFeature` - For new features
- `/ReviewCode` - Review the fix
- `/GenerateTests` - Generate test plan

## Tasks Invoked

- `bug-workflow.investigate`
- `context-loader.loadForBug`
- `mcp-sync.syncLinear`

## Agents Used

- `code-writer` - Fix generation
