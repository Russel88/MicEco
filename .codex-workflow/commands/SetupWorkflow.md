---
name: SetupWorkflow
description: Interactively set up or update the Claude Workflow System in a project
user_invocable: true
---

# /SetupWorkflow

Interactively set up the Claude Workflow System in a new or existing project. This command guides you through configuration options and handles file copying/merging.

## Usage

```
/SetupWorkflow [target-path]
```

## Arguments

- `target-path`: Path to the project to set up (optional, will prompt if not provided)

## Interactive Setup Process

When invoked, this command will:

1. **Determine Target Project**
   - If path provided, use that
   - Otherwise, ask user for the target project path
   - Validate the path exists or offer to create it

2. **Detect Existing Setup**
   - Check if `.claude/` directory exists
   - Check if `CLAUDE.md` exists
   - If existing setup found, ask: Update, Replace, or Cancel

3. **Select Variant**
   - List available variants from `/variants/`
   - Show description of each variant
   - Allow "base only" option for generic projects
   - **Create New Variant** option if no suitable variant exists
     - Invokes `/CreateVariant` to build custom variant
     - Returns to setup with new variant selected

4. **Gather Project Information**
   - Project name (default: directory name)
   - Project description
   - Any variant-specific prompts from manifest

5. **Execute Setup**
   - Copy base layer components
   - Apply variant layer if selected
   - Process templates with provided values
   - Set up GitHub workflows

6. **Configure Project Integrations**
   - Detect GitHub repository from git remote (auto-configure)
   - Prompt for Linear team and project selection (via MCP)
   - Prompt for Coda document and page selection (via MCP)
   - Store configuration in `.claude/settings.json` under `projectMetadata`
   - Skip any integration if MCP not available (can configure later)

7. **Post-Setup**
   - Show summary of files created/modified
   - Show configured integrations status
   - Provide next steps guidance

## Setup Modes

### New Project Setup
For empty directories or projects without `.claude/`:
- Creates full directory structure
- Copies all components
- Generates CLAUDE.md from template
- Runs full project integration setup (Linear, Coda, GitHub)

### Update Existing Project
For projects with existing workflow setup:
- Syncs new/updated components
- Preserves local customizations in CLAUDE.md
- Merges settings where possible
- Checks if project metadata is configured
- Offers to configure missing integrations (skips already configured ones)

### Replace Existing Setup
For projects where you want to start fresh:
- Backs up existing `.claude/` to `.claude.backup/`
- Removes old setup
- Performs fresh installation
- Runs full project integration setup

## Workflow Repository Location

This command expects to be run from within the claude-workflows repository, or requires the `CLAUDE_WORKFLOW_REPO` environment variable to be set to the repository path.

## Examples

```bash
# Set up a new project interactively
/SetupWorkflow /path/to/my-new-project

# Set up current directory (when in target project)
/SetupWorkflow .

# Will prompt for all options
/SetupWorkflow
```

## Implementation

This command uses Claude Code's native file operations to:
1. Read variant manifests and templates
2. Copy files to target project
3. Process template variables
4. Create directory structures
5. Configure project integrations via MCP

No external scripts or dependencies required.

## Creating Custom Variants

If your project type isn't supported by existing variants:

1. During variant selection, choose **"Create New Variant"**
2. Follow the `/CreateVariant` prompts to define your variant
3. Setup continues with your new variant
4. Optionally submit to community with `/SubmitVariant`

Example flow:
```
Variant selection:
  1. base (Generic - any project)
  2. nextjs-development (Next.js with React)
  → 3. Create New Variant...

Creating new variant...
  Name: swift-ios
  Description: iOS development with Swift and SwiftUI
  Project type: swift
  ...

Variant created! Continuing setup with swift-ios...
```

## Tasks Invoked

- `project-metadata.setup` - Configure Linear, Coda, GitHub integrations
- `mcp-sync.verifyConnections` - Check MCP availability before integration setup
- `variant-management.create` - Create new variant (if selected)
- `variant-management.validate` - Validate variant before use

## Related Commands

- `/SetupProjectMeta` - Reconfigure project integrations after initial setup
- `/CreateVariant` - Create a new variant for unsupported project types
- `/SubmitVariant` - Submit variant to community via PR
- `/StartSession` - Begin development session (verifies integrations)
- `/CheckMCPStatus` - Verify MCP connection status
