---
name: ci-integration
description: Coordinate local Claude Code sessions with GitHub Actions CI/CD
---

# CI Integration Task

Bridges local development workflows with GitHub Actions automation, ensuring consistency between local quality checks and CI pipelines.

## Purpose

- Run local checks that mirror CI workflows before pushing
- Analyze CI results after GitHub Actions run
- Maintain consistency between local commands and CI workflows
- Track CI-related issues in task management

## Operations

### `preCommit`
Run local quality checks that mirror GitHub Actions CI before pushing.

**Steps:**
1. Identify changed files via `git diff --name-only`
2. Run local security audit on changed files (mirrors `security-review.yml`)
3. Run local code review on changed files (mirrors `code-review.yml`)
4. Report any findings that would fail CI

**Output:**
```json
{
  "ready": true|false,
  "securityFindings": [...],
  "reviewFindings": [...],
  "recommendation": "string"
}
```

### `reviewCIResults`
Analyze GitHub Actions results after CI runs.

**Steps:**
1. Fetch GitHub Actions run status via `gh run list`
2. If failed, fetch logs via `gh run view --log`
3. Parse security findings from artifacts
4. Prioritize fixes based on severity
5. Optionally create task management issues for findings

**Inputs:**
- `runId`: GitHub Actions run ID (optional, defaults to latest)
- `createIssues`: Whether to create issues for findings

**Output:**
```json
{
  "status": "success|failure",
  "securityFindings": [...],
  "reviewComments": [...],
  "suggestedFixes": [...]
}
```

### `syncWorkflows`
Ensure local Claude Code commands align with CI workflows.

**Steps:**
1. Read `.github/workflows/*.yml`
2. Compare with `.claude/commands/` and `.claude/tasks/`
3. Identify any drift between local and CI configurations
4. Report recommendations for alignment

### `triggerReview`
Trigger GitHub Actions review workflows programmatically.

**Steps:**
1. Verify GitHub CLI (`gh`) is authenticated
2. Check for open PRs on current branch
3. For each requested review type:
   - `codeReview`: Run `gh workflow run code-review.yml`
   - `securityReview`: Run `gh workflow run security-review.yml`
4. Report workflow run URLs

**Inputs:**
- `codeReview`: Whether to trigger code review workflow (boolean)
- `securityReview`: Whether to trigger security review workflow (boolean)
- `prNumber`: Optional PR number (defaults to current branch's PR)

**Output:**
```json
{
  "triggered": ["code-review", "security-review"],
  "workflowRuns": [
    {"workflow": "code-review", "url": "https://github.com/..."},
    {"workflow": "security-review", "url": "https://github.com/..."}
  ],
  "prNumber": 123
}
```

**Implementation:**
```bash
# Trigger code review workflow
gh workflow run code-review.yml -f pr_number=<PR_NUMBER>

# Trigger security review workflow
gh workflow run security-review.yml -f pr_number=<PR_NUMBER>

# Get workflow run URL
gh run list --workflow=code-review.yml --limit=1 --json url
```

## CI/CD Mapping

| Local Command | GitHub Action | Notes |
|--------------|---------------|-------|
| `/ReviewCode` | `code-review.yml` | Same review criteria |
| `/RunSecurityAudit` | `security-review.yml` | Same vulnerability checks |
| `/EndSession codeReview=true` | `code-review.yml` | Triggers on session end |
| `/EndSession securityReview=true` | `security-review.yml` | Triggers on session end |
| `/GenerateTests` | N/A | Local only (test specs) |
| `/ImplementFeature` | N/A | Local only (implementation) |
| `/FixBug` | N/A | Local only (debugging) |

## GitHub Actions Reference

### Security Review (`security-review.yml`)
- **Action**: `anthropics/claude-code-security-review@main`
- **Model**: Claude Opus (deep analysis)
- **Triggers**: PR to main (code files only)
- **Excludes**: `node_modules`, `.next`, `dist`, `coverage`

### Code Review (`code-review.yml`)
- **Action**: `anthropics/claude-code-action@v1`
- **Model**: Claude Sonnet
- **Triggers**: PR opened/synced, @claude mentions
- **Focus**: Architecture, TypeScript, framework patterns, security, performance

## Environment Requirements

- `ANTHROPIC_API_KEY` in GitHub Secrets
- GitHub CLI (`gh`) for local CI interaction
- Task management MCP for issue tracking (optional)

## Collaboration

- Invoked by `/ReviewCode` command for pre-push checks
- Invoked by `/EndSession` command when `codeReview` or `securityReview` flags are set
- Provides CI context to `quality-gates` task
- Updates task management via `mcp-sync` task when issues found
