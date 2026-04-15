---
name: feature-workflow
description: Orchestrate feature implementation from requirements to code delivery
---

# Feature Workflow Task

End-to-end orchestration for feature implementation, coordinating requirement discovery, specification generation, database design, code generation, and task management updates.

## Operations

### `implement`

Full feature implementation workflow.

**Inputs:**
- `feature`: Feature name or description
- `linearIssue`: Optional task management issue ID
- `stack`: Technology requirements (e.g., `["database", "api", "ui"]`)
- `skipSpec`: Skip spec generation if already exists (default: false)

**Steps:**
1. **Discovery Phase:**
   - Invoke `requirement-discovery.discover(query: feature)`
   - Invoke `context-loader.loadForFeature(featureName: feature)`
   - If issue provided, fetch full details

2. **Specification Phase (unless skipSpec):**
   - Check for existing spec in `/docs/specs/`
   - If missing, generate feature specification:
     - Architecture alignment
     - File structure plan
     - API design (if applicable)
     - Database changes (if applicable)
   - Write spec to `/docs/specs/{feature-name}.md`

3. **Approval Checkpoint:**
   - Present specification summary
   - List files to create/modify
   - Wait for user confirmation
   - If rejected, return to specification phase

4. **Implementation Phase:**
   - If `database` in stack:
     - Route to `schema-designer` agent for schema design
     - Generate migration (if approved)
   - Route to `code-writer` agent with:
     - Feature spec
     - Architecture constraints
     - File paths to create/modify
   - Collect generated code

5. **Integration Phase:**
   - Verify code aligns with spec
   - Check for missing implementations
   - Identify edge cases not covered

6. **Update Phase:**
   - Update task management issue with progress via `mcp-sync`
   - Update session state with completed work
   - Generate implementation summary

**Outputs:**
```json
{
  "feature": "User Profile Integration",
  "status": "completed|partial|blocked",
  "specification": {
    "path": "/docs/specs/user-profile.md",
    "generated": true
  },
  "implementation": {
    "filesCreated": [...],
    "filesModified": [...],
    "linesOfCode": 250
  },
  "linearUpdates": {
    "issueId": "PROJ-211",
    "newState": "In Review"
  },
  "nextSteps": [...]
}
```

### `plan`

Generate implementation plan without executing.

**Inputs:**
- `feature`: Feature name or description
- `linearIssue`: Optional issue ID

**Steps:**
1. Run discovery phase only
2. Generate specification
3. Create implementation plan
4. Return plan without executing

**Outputs:**
```json
{
  "feature": "User Profile Integration",
  "plan": {
    "phases": [
      {
        "name": "Database",
        "tasks": [...],
        "agent": "schema-designer"
      },
      {
        "name": "API",
        "tasks": [...],
        "agent": "code-writer"
      }
    ],
    "estimatedFiles": 5,
    "dependencies": [...]
  }
}
```

### `resume`

Resume a partially completed feature.

**Inputs:**
- `feature`: Feature name
- `fromPhase`: Phase to resume from

**Steps:**
1. Load existing spec from `/docs/specs/`
2. Check session state for progress
3. Identify remaining work
4. Continue from specified phase

**Outputs:**
Same as `implement`

## Workflow Diagram

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  Discovery  │────▶│   Spec Gen  │────▶│  Approval   │
└─────────────┘     └─────────────┘     └─────────────┘
                                              │
                    ┌─────────────────────────┘
                    ▼
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  Database   │────▶│    Code     │────▶│   Update    │
│  (if needed)│     │  Generation │     │   Tracking  │
└─────────────┘     └─────────────┘     └─────────────┘
```

## Agent Routing

| Requirement Type | Agent | Input |
|------------------|-------|-------|
| Database schema | `schema-designer` | Requirements, existing schema |
| API endpoints | `code-writer` | Spec, route conventions |
| UI components | `code-writer` | Spec, design system |
| Tests | `test-planner` | Implementation, acceptance criteria |

## Specification Template

Generated specs follow this structure:

```markdown
# Feature: {name}

## Overview
{description}

## Requirements
{from requirement-discovery}

## Architecture
### Files to Create
- path/to/file.ts - purpose

### Files to Modify
- existing/file.ts - changes needed

## Database Changes
{if applicable}

## API Design
{if applicable}

## UI Components
{if applicable}

## Testing Strategy
{outline}

## Acceptance Criteria
{from requirements}
```

## Dependencies

- **requirement-discovery**: For requirements
- **context-loader**: For project context
- **mcp-sync**: For task management updates
- **session-management**: For state updates

## Agents Used

- **schema-designer**: Database design
- **code-writer**: All code generation

## Error Handling

| Error | Action |
|-------|--------|
| Requirements not found | Prompt for manual input |
| Spec generation fails | Save partial, flag for review |
| Agent produces invalid code | Flag errors, request fixes |
| Task management update fails | Continue, retry at end |

## Approval Checkpoints

The workflow pauses for user approval at:
1. After specification generation
2. Before database migrations
3. Before major file creation

Approvals prevent unintended changes and ensure alignment.
