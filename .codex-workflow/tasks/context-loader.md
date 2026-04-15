---
name: context-loader
description: Load and bundle project context for agent consumption
---

# Context Loader Task

Gathers and bundles all relevant project context for use by agents and other tasks. Ensures consistent context loading across all operations.

## Operations

### `loadFull`

Load complete project context.

**Steps:**
1. Load core configuration:
   - Read `CLAUDE.md` for architecture rules
   - Read `.claude/agents/*.md` for agent capabilities
   - Read `.claude/commands/*.md` for available commands
   - Read `.claude/tasks/*.md` for task definitions
2. Load documentation context:
   - Scan `/knowledge/prd/` for PRD files
   - Scan `/docs/specs/` for active specifications
   - Read `/docs/planning/CURRENT-STATE.md` if exists
3. Load git context:
   - Current branch name
   - Recent commits (last 10)
   - Uncommitted changes summary
   - Stash list
4. Load technology context:
   - Parse `package.json` for dependencies
   - Check runtime version
   - Identify framework versions
5. Bundle and return

**Outputs:**
```json
{
  "architecture": {
    "rules": "CLAUDE.md content summary",
    "agents": ["code-writer", "schema-designer", ...],
    "commands": ["StartSession", "EndSession", ...],
    "tasks": ["session-management", "mcp-sync", ...]
  },
  "documentation": {
    "prds": ["/knowledge/prd/consolidated-prd.md", ...],
    "specs": ["/docs/specs/feature-spec.md", ...],
    "currentState": "CURRENT-STATE.md content if exists"
  },
  "git": {
    "branch": "feature/implementation",
    "recentCommits": [...],
    "uncommittedChanges": 5,
    "stashCount": 2
  },
  "technology": {
    "runtime": "node/bun version",
    "framework": "next/react version",
    "typescript": "version"
  },
  "loadedAt": "2026-01-16T10:00:00Z"
}
```

### `loadMinimal`

Load minimal context for quick operations.

**Steps:**
1. Read `CLAUDE.md` summary section only
2. Get current git branch
3. Check for session state
4. Return minimal bundle

**Outputs:**
```json
{
  "branch": "feature/implementation",
  "hasActiveSession": true,
  "sessionId": "session-2026-01-16-1000",
  "loadedAt": "2026-01-16T10:00:00Z"
}
```

### `loadForFeature`

Load context relevant to a specific feature.

**Inputs:**
- `featureName`: Name or identifier of feature
- `includeSpecs`: Whether to include full spec content (default: true)

**Steps:**
1. Search `/knowledge/prd/` for matching PRD
2. Search `/docs/specs/` for matching specifications
3. Identify related files in codebase (by naming convention)
4. Load relevant issues via `mcp-sync`
5. Bundle feature-specific context

**Outputs:**
```json
{
  "feature": "feature-name",
  "prd": {
    "path": "/knowledge/prd/feature-prd.md",
    "relevantSections": [...]
  },
  "specs": [
    {
      "path": "/docs/specs/feature-spec.md",
      "content": "..."
    }
  ],
  "codeLocations": [
    "app/feature/",
    "lib/feature/"
  ],
  "linearIssues": ["PROJ-211", "PROJ-212"],
  "loadedAt": "2026-01-16T10:00:00Z"
}
```

### `loadForBug`

Load context relevant to bug investigation.

**Inputs:**
- `bugDescription`: Description of the bug
- `affectedFiles`: Known affected files (optional)
- `linearIssue`: Issue ID (optional)

**Steps:**
1. If issue provided, fetch issue details via `mcp-sync`
2. Identify potentially affected code paths
3. Load recent git history for affected files
4. Check for related test files
5. Bundle bug investigation context

**Outputs:**
```json
{
  "bug": {
    "description": "...",
    "linearIssue": "PROJ-XXX",
    "severity": "High"
  },
  "affectedFiles": [
    {
      "path": "lib/feature/handler.ts",
      "recentChanges": [...]
    }
  ],
  "relatedTests": [
    "lib/feature/__tests__/handler.test.ts"
  ],
  "loadedAt": "2026-01-16T10:00:00Z"
}
```

### `codebasePatterns`

Search for patterns in the codebase relevant to a PRD.

**Inputs:**
- `keywords`: Array of technical keywords from PRD
- `areas`: Code areas to search (optional)

**Steps:**
1. Search codebase for each keyword
2. Identify existing implementations of similar functionality
3. Find reusable utilities, hooks, components
4. Map file relationships and dependencies
5. Return pattern analysis

**Outputs:**
```json
{
  "patterns": [
    {
      "keyword": "authentication",
      "existingImplementations": ["lib/auth/session.ts", "middleware/auth.ts"],
      "reusableComponents": ["components/LoginForm.tsx"],
      "relatedTests": ["__tests__/auth.test.ts"]
    }
  ],
  "suggestedApproach": "Extend existing auth module",
  "loadedAt": "2026-01-16T10:00:00Z"
}
```

### `loadPRDContext`

Load PRD and all enrichment data for a session.

**Inputs:**
- `prdSlug`: PRD feature slug
- `includeEnrichment`: Include enrichment data (default: true)
- `includeValidation`: Include validation report (default: true)

**Steps:**
1. Load PRD from /docs/planning/prds/{prdSlug}.md
2. Load validation report if exists
3. Load enrichment data if exists
4. Load feasibility analysis if exists
5. Load breakdown if exists
6. Bundle all PRD context

**Outputs:**
```json
{
  "prd": { "slug": "feature-name", "content": "...", "status": "enriched" },
  "validation": { "score": 85, "issues": [] },
  "enrichment": { "complexity": "M", "risks": [] },
  "feasibility": { "rating": "HIGH" },
  "breakdown": { "itemCount": 12 },
  "loadedAt": "2026-01-16T10:00:00Z"
}
```

## Context Sources

| Source | Location | Purpose |
|--------|----------|---------|
| CLAUDE.md | Root | Architecture rules, conventions |
| Agents | `.claude/agents/` | Available agent capabilities |
| Commands | `.claude/commands/` | Available user commands |
| Tasks | `.claude/tasks/` | Available orchestration tasks |
| PRDs | `/knowledge/prd/` | Requirements documentation |
| Specs | `/docs/specs/` | Technical specifications |
| Current State | `/docs/planning/CURRENT-STATE.md` | Active work summary |
| Session State | `/docs/planning/session-state.json` | Machine-readable state |
| package.json | Root | Dependencies and versions |
| Git | Repository | Branch, commits, changes |

## Exclusions

Always exclude from context loading:
- `node_modules/`
- `.next/`
- `dist/`
- `*.log` files
- `.env*` files (security)
- Large binary files

## Caching

Context is cached for the duration of a session. Cache invalidation:
- On file modification in watched directories
- On explicit refresh request
- On session end

## Dependencies

- **mcp-sync**: For issue fetching (optional)

## Error Handling

| Error | Action |
|-------|--------|
| File not found | Skip file, log warning |
| Parse error | Return raw content, log error |
| Permission denied | Skip file, log error |
| Large file | Truncate to first 10KB, flag as truncated |
