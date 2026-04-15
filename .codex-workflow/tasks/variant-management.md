---
name: variant-management
description: Create, validate, and submit workflow variants via branches and pull requests
---

# Variant Management Task

Handles the lifecycle of workflow variants including creation, validation, and contribution back to the main repository via proper branching and pull request workflows.

## Purpose

Enable users to:
- Create new variants for unsupported project types (e.g., Swift, Flutter, Rust)
- Customize existing variants with additional components
- Contribute variants and changes back to the community via pull requests
- Never push directly to main - always use branches and merge requests

## Operations

### `create`

Create a new variant from scratch.

**Inputs:**
- `name`: Variant name (lowercase, hyphens) e.g., "swift-ios"
- `description`: Human-readable description
- `projectType`: Project type (e.g., "swift", "flutter", "rust", "generic")
- `baseVersion`: Base workflow version to extend (default: current)

**Steps:**
1. Validate variant name is unique and follows conventions
2. Create variant directory structure:
   ```
   variants/{name}/
   ├── manifest.json
   ├── README.md
   ├── agents/
   ├── commands/
   ├── tasks/
   └── templates/
   ```
3. Generate manifest.json with provided metadata
4. Create README.md template
5. Create placeholder CLAUDE.md.variant template
6. Return created file paths

**Outputs:**
```json
{
  "success": true,
  "variantPath": "variants/swift-ios",
  "files": [
    "manifest.json",
    "README.md",
    "templates/CLAUDE.md.variant"
  ],
  "nextSteps": [
    "Add agents to agents/",
    "Add commands to commands/",
    "Add tasks to tasks/",
    "Update manifest.json with component references",
    "Run /SubmitVariant when ready"
  ]
}
```

### `addAgent`

Add a new agent to a variant.

**Inputs:**
- `variantName`: Target variant name
- `agentName`: Agent identifier (lowercase, hyphens)
- `description`: Agent description
- `expertiseAreas`: Array of expertise areas

**Steps:**
1. Validate variant exists
2. Create agent markdown file from template
3. Update manifest.json to include new agent
4. Return file path

**Outputs:**
```json
{
  "success": true,
  "filePath": "variants/swift-ios/agents/swift-specialist.md",
  "manifestUpdated": true
}
```

### `addCommand`

Add a new command to a variant.

**Inputs:**
- `variantName`: Target variant name
- `commandName`: Command name (PascalCase)
- `description`: Command description
- `usage`: Usage pattern

**Steps:**
1. Validate variant exists
2. Create command markdown file from template
3. Update manifest.json to include new command
4. Return file path

### `addTask`

Add a new task to a variant.

**Inputs:**
- `variantName`: Target variant name
- `taskName`: Task identifier (lowercase, hyphens)
- `description`: Task description
- `operations`: Array of operation names

**Steps:**
1. Validate variant exists
2. Create task markdown file from template
3. Update manifest.json to include new task
4. Return file path

### `validate`

Validate a variant structure and manifest.

**Inputs:**
- `variantName`: Variant to validate

**Steps:**
1. Check manifest.json exists and is valid JSON
2. Validate required fields (name, version, description, projectType)
3. Verify referenced files exist
4. Check markdown frontmatter is valid
5. Validate version format
6. Return validation results

**Outputs:**
```json
{
  "valid": true,
  "errors": [],
  "warnings": [
    "README.md is minimal - consider adding more documentation"
  ],
  "components": {
    "agents": 2,
    "commands": 1,
    "tasks": 1,
    "templates": 1
  }
}
```

### `prepareBranch`

Prepare a branch for submitting variant changes.

**Inputs:**
- `variantName`: Variant to submit
- `changeType`: "new" | "update"
- `description`: Brief description of changes

**Steps:**
1. Validate variant passes validation
2. Determine branch name: `variant/{variantName}` (new) or `variant/{variantName}-update-{date}` (update)
3. Check if branch already exists
4. Create new branch from main
5. Stage variant files
6. Create commit with descriptive message
7. Return branch info

**Outputs:**
```json
{
  "success": true,
  "branch": "variant/swift-ios",
  "baseBranch": "main",
  "stagedFiles": [
    "variants/swift-ios/manifest.json",
    "variants/swift-ios/agents/swift-specialist.md",
    "..."
  ],
  "commitMessage": "feat(variant): Add swift-ios variant for iOS development"
}
```

### `submitPullRequest`

Create a pull request for variant contribution.

**Inputs:**
- `branch`: Branch name with variant changes
- `title`: PR title
- `description`: PR description
- `variantName`: Variant being submitted

**Steps:**
1. Verify branch exists and has commits
2. Push branch to origin
3. Create pull request via GitHub CLI:
   - Target: main branch
   - Include variant summary in description
   - Add appropriate labels
4. Return PR URL

**Outputs:**
```json
{
  "success": true,
  "prUrl": "https://github.com/Clinical-Microbiomics/claude-workflow/pull/42",
  "prNumber": 42,
  "branch": "variant/swift-ios"
}
```

### `listLocalVariants`

List variants available locally (including custom ones).

**Steps:**
1. Scan `variants/` directory
2. Load manifest.json from each
3. Identify which are from upstream vs local
4. Return list

**Outputs:**
```json
{
  "variants": [
    {
      "name": "nextjs-development",
      "version": "1.0.0",
      "description": "...",
      "source": "upstream"
    },
    {
      "name": "swift-ios",
      "version": "1.0.0",
      "description": "...",
      "source": "local",
      "submitted": false
    }
  ]
}
```

## Branch Naming Conventions

| Change Type | Branch Pattern | Example |
|-------------|----------------|---------|
| New variant | `variant/{name}` | `variant/swift-ios` |
| Update variant | `variant/{name}-update-{YYYYMMDD}` | `variant/nextjs-development-update-20260117` |
| Bug fix | `variant/{name}-fix-{issue}` | `variant/nextjs-development-fix-123` |

## Commit Message Format

```
feat(variant): Add {variant-name} variant

- Description of what this variant provides
- Key agents/commands/tasks included
- Target project types

Components:
- Agents: {count}
- Commands: {count}
- Tasks: {count}
```

## Safety Rules

1. **Never push to main** - All changes go through branches and PRs
2. **Validate before submit** - Variants must pass validation
3. **No force push** - Preserve history on shared branches
4. **Require PR review** - Changes need approval before merge

## Error Handling

| Error | Action |
|-------|--------|
| Variant name exists | Suggest alternative name or update flow |
| Invalid manifest | Return validation errors, don't proceed |
| Branch already exists | Offer to use existing or create new |
| GitHub CLI not authenticated | Prompt user to run `gh auth login` |
| Push failed | Check permissions, suggest fork workflow |

## Dependencies

- GitHub CLI (`gh`) for PR creation
- Git for branch management

## Usage Examples

```
// Create new variant
invoke variant-management.create(
  name="swift-ios",
  description="iOS development with Swift and SwiftUI",
  projectType="swift"
)

// Add an agent
invoke variant-management.addAgent(
  variantName="swift-ios",
  agentName="swift-specialist",
  description="Expert in Swift and iOS development"
)

// Validate variant
invoke variant-management.validate(variantName="swift-ios")

// Submit via PR
invoke variant-management.prepareBranch(
  variantName="swift-ios",
  changeType="new"
)
invoke variant-management.submitPullRequest(
  branch="variant/swift-ios",
  title="Add Swift iOS variant",
  variantName="swift-ios"
)
```
