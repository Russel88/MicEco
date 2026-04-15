---
description: Generate comprehensive test plans with unit, integration, and edge case coverage
---

# GenerateTests

Generate test specifications and plans for code coverage.

## Purpose

Creates comprehensive test plans including unit tests, integration tests, and edge case coverage. Outputs test specifications, not test code.

## Arguments

- `target`: File, function, or feature to test (required)
- `types`: Test types to generate (optional, multi-select)
  - `unit` - Unit tests
  - `integration` - Integration tests
  - `e2e` - End-to-end tests
- `coverage`: Target coverage percentage (default: 80)

## Execution

1. Invoke `quality-gates` task with action: `generateTests`
   - Loads target context
   - Routes to test-planner agent
   - Generates test specifications
   - Identifies edge cases
   - Returns test plan

2. Display:
   - Test plan summary
   - Test cases by type
   - Edge cases to cover
   - Mock requirements
   - Estimated coverage

## Prerequisites

- Target code exists

## Output Format

```
### Test Plan: lib/api/handler.ts

#### Unit Tests (8 cases)
1. handleRequest returns data for valid input
2. handleRequest throws ValidationError for invalid input
3. handleRequest handles network timeout
...

#### Integration Tests (3 cases)
1. API flow from request to database
...

#### Edge Cases
- Network timeout after partial response
- Rate limiting (429 response)
- Invalid JSON in response
...

#### Mock Requirements
- fetch: Mock network requests
- database: Mock database calls
```

## Example

```
/GenerateTests target="lib/api/handler.ts" types=["unit", "integration"] coverage=90
```

## Test Types

| Type | Scope | Dependencies |
|------|-------|--------------|
| Unit | Single function/class | Mocked |
| Integration | Multiple modules | Partial mocks |
| E2E | Full user flow | None (real) |

## Related

- `/ReviewCode` - Review test coverage
- `/ImplementFeature` - Implement features
- `/FixBug` - Fix bugs with test verification

## Tasks Invoked

- `quality-gates.generateTests`
- `context-loader.loadFull`

## Agents Used

- `test-planner` - Test specification

## Notes

This command generates test **specifications**, not test **code**. Use the specifications to implement tests manually or with the code-writer agent.
