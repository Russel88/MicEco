---
description: Plan and execute code refactoring with structured approach and safety checks
---

# RefactorCode

Structured code refactoring with planning, safety checks, and verification.

## Purpose

Guides refactoring efforts from analysis through implementation, ensuring code quality improvements without introducing regressions.

## Arguments

- `target`: File, module, or pattern to refactor (required)
- `goal`: Refactoring goal (optional)
  - `simplify` - Reduce complexity
  - `extract` - Extract reusable components
  - `consolidate` - Merge duplicate code
  - `modernize` - Update to current patterns
- `linearIssue`: Task management issue ID (optional)

## Execution

1. Invoke `feature-workflow` task with action: `plan` (refactor mode)
   - Analysis: Assess current code structure
   - Plan: Generate refactoring plan
   - Approval: Present plan for confirmation
   - Implementation: Execute refactoring steps
   - Verification: Ensure no regressions

2. Display:
   - Current code analysis
   - Refactoring plan
   - Risk assessment
   - Implementation steps

## Prerequisites

- Target code exists
- Understanding of current functionality

## Example

```
/RefactorCode target="lib/utils/" goal="consolidate" linearIssue="PROJ-220"
```

## Refactoring Types

| Goal | Description | Risk Level |
|------|-------------|------------|
| `simplify` | Reduce complexity, improve readability | Low |
| `extract` | Extract shared logic to utilities | Medium |
| `consolidate` | Merge duplicate implementations | Medium |
| `modernize` | Update to current framework patterns | High |

## Safety Checks

Before refactoring:
- Identify all usages of target code
- Check for test coverage
- Note external dependencies

After refactoring:
- Verify all usages still work
- Run existing tests
- Manual verification of critical paths

## Related

- `/ImplementFeature` - For new features
- `/ReviewCode` - Review the refactoring
- `/GenerateTests` - Add test coverage first

## Tasks Invoked

- `feature-workflow.plan` (refactor mode)
- `context-loader.loadFull`
- `quality-gates.review`

## Agents Used

- `code-writer` - Refactoring implementation

## Approval Checkpoint

The workflow pauses for approval after the refactoring plan is generated, before any code changes.
