---
name: SetupProjectMeta
description: Reconfigure project integration metadata for Linear, Coda, and GitHub
user_invocable: true
---

# /SetupProjectMeta

Reconfigure the project's integration metadata for Linear (task management), Coda (documentation), and GitHub (source control). Use this command to change integrations after initial setup or to configure integrations that were skipped.

## Purpose

Project integrations are normally configured during `/SetupWorkflow`. Use this command when you need to:
- Change the Linear team or project
- Switch to a different Coda document or page
- Update the GitHub repository association
- Configure integrations that were skipped during initial setup

The workflow system uses this metadata to:
- Create and query issues in the correct Linear project
- Read and update product requirements in the correct Coda page
- Verify GitHub repository access

## Usage

```
/SetupProjectMeta [options]
```

## Arguments

- `--force`: Overwrite existing metadata configuration
- `--verify-only`: Only verify existing configuration without prompting for changes
- `--linear-only`: Only configure Linear integration
- `--coda-only`: Only configure Coda integration
- `--github-only`: Only configure GitHub integration

## Interactive Setup Process

When invoked, this command will:

### 1. Check Current Configuration

- Load `.claude/settings.json`
- Display current metadata status
- If already configured, ask: Update or Keep existing

### 2. Configure Linear Integration

1. Verify Linear MCP connection
2. List available teams via `mcp__linear__list_teams`
3. Prompt user to select team
4. List projects in selected team via `mcp__linear__list_projects`
5. Prompt user to select or create project
6. If creating new project, prompt for project name and details
7. Store team and project IDs in settings

### 3. Configure Coda Integration

1. Verify Coda MCP connection
2. List available documents via `mcp__coda__coda_list_documents`
3. Prompt user to select document
4. List pages in selected document via `mcp__coda__coda_list_pages`
5. Prompt user to select or create project requirements page
6. If creating new page, use PRD template
7. Store document and page IDs in settings

### 4. Configure GitHub Integration

1. Attempt auto-detection from git remote
2. If detected, confirm with user
3. If not detected or user wants different repo, prompt for repo URL
4. Verify access via `gh repo view`
5. Store repository details in settings

### 5. Finalize

- Save updated settings to `.claude/settings.json`
- Display summary of configured integrations
- Provide next steps guidance

## Prerequisites

- Linear MCP server configured and authenticated
- Coda MCP server configured and authenticated (optional)
- GitHub CLI (`gh`) authenticated
- `.claude/settings.json` exists (created by `/SetupWorkflow`)

## Output

Updates `.claude/settings.json` with:

```json
{
  "projectMetadata": {
    "linear": {
      "teamId": "team-uuid",
      "teamName": "Engineering",
      "projectId": "project-uuid",
      "projectName": "My Project"
    },
    "coda": {
      "docId": "document-id",
      "docName": "Project Documentation",
      "pageId": "page-id",
      "pageName": "My Project Requirements"
    },
    "github": {
      "repo": "org/repo-name",
      "owner": "org",
      "name": "repo-name"
    }
  }
}
```

## Examples

```bash
# Full interactive setup
/SetupProjectMeta

# Verify existing configuration
/SetupProjectMeta --verify-only

# Only set up Linear integration
/SetupProjectMeta --linear-only

# Force reconfigure all integrations
/SetupProjectMeta --force
```

## Execution Steps

1. Load current settings from `.claude/settings.json`
2. Check MCP connection status via `mcp-sync.verifyConnections`
3. For each integration (Linear, Coda, GitHub):
   - Check if already configured
   - If not configured or `--force`, run interactive setup
   - Verify configuration is valid
4. Save updated settings
5. Display configuration summary

## Error Handling

| Error | Action |
|-------|--------|
| Settings file missing | Create default settings file first |
| Linear MCP not connected | Skip Linear setup, warn user |
| Coda MCP not connected | Skip Coda setup, warn user |
| GitHub CLI not authenticated | Skip GitHub setup, warn user |
| Team/project not found | Clear invalid config, re-prompt |

## Related Commands

- `/StartSession` - Verifies metadata before starting a session
- `/CheckMCPStatus` - View connection status for all MCPs
- `/SyncLinear` - Manually sync with Linear (uses configured project)
- `/SetupWorkflow` - Initial workflow system setup

## Tasks Invoked

- `project-metadata.setup`
- `project-metadata.verify`
- `mcp-sync.verifyConnections`

## Notes

- This command should be run once when setting up a new project
- Re-run with `--force` if you need to change project associations
- Linear and Coda configurations are optional but recommended
- GitHub is auto-detected but can be overridden
