---
name: workflow-submit-variant
description: Submit a variant or variant changes to the community via branch and pull request Use when the user asks to submit variant in this repository or refers to the old `/SubmitVariant` command.
---

# Submit Variant

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# `workflow-submit-variant`

Submit a new variant or changes to an existing variant to the claude-workflow repository via a proper branching and pull request workflow. This command ensures contributions follow the repository's contribution guidelines.

## Purpose

Contribute your variants and improvements back to the community:
- Submit new variants you've created
- Propose updates to existing variants
- Fix bugs in variant components
- All contributions go through branches and PRs - never direct push to main

## Usage

```
`workflow-submit-variant` <variant-name> [options]
```

## Arguments

- `variant-name`: The variant to submit (required)

## Options

- `--type`: Submission type: `new` | `update` | `fix` (default: auto-detect)
- `--title`: Custom PR title
- `--dry-run`: Preview what would be submitted without actually creating PR

## Submission Process

### 1. Pre-flight Checks

Before submission, the command verifies:
- Variant exists and passes validation
- GitHub CLI is authenticated (`gh auth status`)
- User has push access or can fork
- No uncommitted changes in variant directory

### 2. Validate Variant

Runs full validation:
```
Validating swift-ios...
  ✓ manifest.json valid
  ✓ All referenced files exist
  ✓ Markdown frontmatter valid
  ✓ Version format correct
  ✓ README.md present
```

### 3. Create Branch

Creates appropriately named branch:

| Submission Type | Branch Name |
|-----------------|-------------|
| New variant | `variant/{name}` |
| Update | `variant/{name}-update-YYYYMMDD` |
| Bug fix | `variant/{name}-fix-{description}` |

### 4. Stage and Commit

- Stages all variant files
- Creates descriptive commit message
- Shows diff summary for confirmation

### 5. Push and Create PR

- Pushes branch to origin (or fork)
- Creates pull request via GitHub CLI
- Includes variant summary in PR description
- Adds appropriate labels

## Examples

```bash
# Submit a new variant
`workflow-submit-variant` swift-ios

# Submit update to existing variant
`workflow-submit-variant` nextjs-development --type update

# Preview without submitting
`workflow-submit-variant` swift-ios --dry-run

# Custom PR title
`workflow-submit-variant` swift-ios --title "Add Swift iOS variant with SwiftUI support"
```

## Output

```
Submitting variant: swift-ios

Pre-flight checks:
  ✓ Variant validation passed
  ✓ GitHub CLI authenticated
  ✓ Repository access confirmed
  ✓ Working directory clean

Creating branch: variant/swift-ios

Files to be committed:
  + variants/swift-ios/manifest.json
  + variants/swift-ios/README.md
  + variants/swift-ios/agents/swift-specialist.md
  + variants/swift-ios/templates/AGENTS.md.variant

Commit message:
  feat(variant): Add swift-ios variant for iOS development

  - Swift and SwiftUI specialist agent
  - iOS-specific project structure
  - Xcode integration patterns

  Components: 1 agent, 0 commands, 0 tasks

Push to origin? [Y/n]: y

Creating pull request...

✓ Pull request created: https://github.com/Clinical-Microbiomics/codex-workflow/pull/42

Next steps:
1. Review your PR at the link above
2. Respond to any review feedback
3. Once approved, maintainers will merge to main
```

## Pull Request Template

The generated PR follows this format:

```markdown
## Summary

Add {variant-name} variant for {project-type} development.

## Components

| Type | Count | Details |
|------|-------|---------|
| Agents | 1 | swift-specialist |
| Commands | 0 | - |
| Tasks | 0 | - |
| Templates | 1 | AGENTS.md.variant |

## Description

{variant-description}

## Checklist

- [x] Variant validation passes
- [x] README.md included
- [x] All components documented
- [ ] Tested with sample project

## Test Instructions

1. Clone this PR branch
2. Run: ``workflow-setup-workflow` /path/to/test --variant {variant-name}`
3. Verify setup completes successfully
```

## Handling Updates to Existing Variants

When submitting changes to variants in the main repository:

1. Ensure you're working with the latest main
2. Make your changes in the variant directory
3. Run ``workflow-submit-variant` variant-name --type update`
4. Branch name includes date to avoid conflicts
5. PR describes what changed

## Fork Workflow

If you don't have push access to the main repository:

1. Command detects lack of push access
2. Offers to fork the repository
3. Creates branch on your fork
4. PR targets the upstream repository

```
You don't have push access to Clinical-Microbiomics/claude-workflow.
Would you like to fork and submit from there? [Y/n]: y

Forking repository...
✓ Fork created: yourusername/claude-workflow

Pushing to fork...
Creating PR to upstream...
```

## Safety Rules

This command enforces safe contribution practices:

1. **Never pushes to main** - Always creates feature branches
2. **Requires validation** - Won't submit invalid variants
3. **Requires clean state** - No uncommitted changes
4. **Descriptive commits** - Auto-generates meaningful commit messages
5. **PR review required** - All changes need maintainer approval

## Error Handling

| Error | Action |
|-------|--------|
| Validation failed | Show errors, abort submission |
| Not authenticated | Prompt to run `gh auth login` |
| Branch exists | Offer to use existing or create new |
| Push rejected | Check permissions, suggest fork |
| PR creation failed | Show error, provide manual instructions |

## Related Commands

- `workflow-create-variant` - Create a new variant
- `workflow-setup-workflow` - Test variant before submitting
- `workflow-check-mcp-status` - Verify GitHub CLI status

## Tasks Invoked

- `variant-management.validate`
- `variant-management.prepareBranch`
- `variant-management.submitPullRequest`

## Prerequisites

- GitHub CLI installed and authenticated (`gh auth login`)
- Git configured with user.name and user.email
- Variant passes validation

## Notes

- First-time contributors may need to fork
- Maintainers review all PRs before merge
- Complex variants may require discussion before submission
- Consider opening an issue first for major new variants
