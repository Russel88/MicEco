---
name: project-metadata
description: Manage project-specific metadata for Linear, Coda, and GitHub integrations
---

# Project Metadata Task

Centralized management of project-specific integration metadata. Handles setup, verification, and retrieval of Linear project, Coda documentation, and GitHub repository configurations.

## Purpose

Every project using the Claude Workflow System needs to be connected to external services:
- **Linear**: Project/team for issue tracking and task management
- **Coda**: Document and page for product requirements and documentation
- **GitHub**: Repository for source control and CI/CD

This task ensures these connections are configured and verified before workflow operations.

## Operations

### `verify`

Check if project metadata is configured and valid.

**Steps:**
1. Load `.claude/settings.json`
2. Check for `projectMetadata` section
3. Validate required fields are present:
   - `linear.teamId` or `linear.projectId`
   - `coda.docId` (optional but recommended)
   - `github.repo` (auto-detected from git remote if not set)
4. Test connections via MCP if configured:
   - Verify Linear team/project exists
   - Verify Coda document is accessible
   - Verify GitHub repo access
5. Return verification status

**Outputs:**
```json
{
  "configured": true,
  "verified": true,
  "linear": {
    "configured": true,
    "teamId": "TEAM-123",
    "teamName": "Engineering",
    "projectId": "PROJ-456",
    "projectName": "My Project",
    "verified": true
  },
  "coda": {
    "configured": true,
    "docId": "abc123",
    "docName": "Project Documentation",
    "pageId": "Requirements",
    "verified": true
  },
  "github": {
    "configured": true,
    "repo": "org/repo-name",
    "verified": true,
    "autoDetected": true
  },
  "missing": []
}
```

### `setup`

Interactive setup of project metadata.

**Steps:**
1. Load existing settings (if any)
2. Detect GitHub repository from git remote
3. Prompt for Linear configuration:
   - List available teams via MCP
   - List projects within selected team
   - Store team and project IDs
4. Prompt for Coda configuration:
   - List available documents via MCP
   - List pages within selected document
   - Optionally create new project page
   - Store document and page IDs
5. Confirm GitHub repository:
   - Auto-detect from git remote
   - Allow override if needed
6. Save to `.claude/settings.json`
7. Return setup summary

**Inputs:**
- `interactive`: Boolean (default: true) - if false, only validate existing config
- `force`: Boolean (default: false) - overwrite existing metadata

**Outputs:**
```json
{
  "success": true,
  "metadata": {
    "linear": { ... },
    "coda": { ... },
    "github": { ... }
  },
  "message": "Project metadata configured successfully"
}
```

### `getLinearContext`

Retrieve Linear project context for issue operations.

**Steps:**
1. Load project metadata
2. Return team and project identifiers

**Outputs:**
```json
{
  "teamId": "TEAM-123",
  "teamName": "Engineering",
  "projectId": "PROJ-456",
  "projectName": "My Project"
}
```

### `getCodaContext`

Retrieve Coda document context for documentation operations.

**Steps:**
1. Load project metadata
2. Return document and page identifiers

**Outputs:**
```json
{
  "docId": "abc123",
  "docName": "Project Documentation",
  "pageId": "page-id",
  "pageName": "Requirements"
}
```

### `getGitHubContext`

Retrieve GitHub repository context.

**Steps:**
1. Load project metadata
2. If not configured, detect from git remote
3. Return repository information

**Outputs:**
```json
{
  "repo": "org/repo-name",
  "owner": "org",
  "name": "repo-name",
  "defaultBranch": "main"
}
```

### `ensureLinearProject`

Ensure Linear project exists, create if needed.

**Steps:**
1. Check if project metadata configured
2. If configured, verify project exists
3. If not exists, offer to create:
   - Use project name from settings
   - Create in configured team
4. Update metadata with project ID
5. Return project details

**Inputs:**
- `createIfMissing`: Boolean (default: true)

**Outputs:**
```json
{
  "exists": true,
  "created": false,
  "project": {
    "id": "PROJ-456",
    "name": "My Project"
  }
}
```

### `ensureCodaPage`

Ensure Coda project page exists, create if needed.

**Steps:**
1. Check if Coda metadata configured
2. If configured, verify page exists
3. If not exists, offer to create:
   - Create page with project name
   - Use PRD template
4. Update metadata with page ID
5. Return page details

**Inputs:**
- `createIfMissing`: Boolean (default: true)
- `template`: String (default: "prd") - template to use for new pages

**Outputs:**
```json
{
  "exists": true,
  "created": false,
  "page": {
    "id": "page-123",
    "name": "My Project Requirements"
  }
}
```

### `ensureGitHubRepo`

Verify GitHub repository access and configuration.

**Steps:**
1. Detect repo from git remote
2. Verify GitHub CLI access (`gh auth status`)
3. Check repository exists and is accessible
4. Return repository details

**Outputs:**
```json
{
  "exists": true,
  "accessible": true,
  "repo": {
    "owner": "org",
    "name": "repo-name",
    "url": "https://github.com/org/repo-name"
  }
}
```

## Configuration Schema

Project metadata is stored in `.claude/settings.json` under the `projectMetadata` key:

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

## Integration with Other Tasks

### Session Management
- `session-management.start` calls `project-metadata.verify`
- Warns if metadata is not configured
- Offers to run setup if missing

### MCP Sync
- `mcp-sync.syncLinear` uses `project-metadata.getLinearContext` for team/project filtering
- `mcp-sync.syncCoda` uses `project-metadata.getCodaContext` for document targeting

### Context Loader
- Includes project metadata in session context
- Provides integration identifiers to agents

## Error Handling

| Error Type | Action |
|------------|--------|
| Settings file missing | Create default settings, prompt for setup |
| Metadata section missing | Prompt for setup |
| Linear connection failed | Mark as unverified, allow offline mode |
| Coda connection failed | Mark as unverified, continue without |
| GitHub not detected | Warn, allow manual configuration |
| Invalid project/team ID | Clear invalid config, prompt for reconfiguration |

## Dependencies

- **mcp-sync**: For verifying connections and listing resources

## Usage Examples

```
// Verify project metadata is configured
invoke project-metadata.verify

// Interactive setup
invoke project-metadata.setup(interactive=true)

// Get Linear context for issue creation
invoke project-metadata.getLinearContext

// Ensure all integrations are ready
invoke project-metadata.ensureLinearProject(createIfMissing=true)
invoke project-metadata.ensureCodaPage(createIfMissing=true)
invoke project-metadata.ensureGitHubRepo
```
