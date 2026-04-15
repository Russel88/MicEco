---
name: git-workflow
description: Manage git operations including branch management and commit workflows
---

# Git Workflow Task

Manages git operations for the development workflow, including branch management, commit preparation, and repository state tracking.

## Operations

### `getState`

Get current repository state.

**Inputs:**
- None

**Steps:**
1. Get current branch name
2. Check for uncommitted changes
3. Check for unpushed commits
4. Get stash list
5. Get recent commits (last 10)

**Outputs:**
```json
{
  "branch": "feature/implementation",
  "isClean": false,
  "uncommittedChanges": {
    "staged": 3,
    "unstaged": 5,
    "untracked": 2
  },
  "unpushedCommits": 2,
  "stashes": [
    {
      "index": 0,
      "message": "WIP: feature work",
      "branch": "feature/other"
    }
  ],
  "recentCommits": [
    {
      "hash": "abc1234",
      "message": "feat(feature): Add implementation",
      "author": "Developer",
      "date": "2026-01-15T18:00:00Z"
    }
  ],
  "remoteStatus": "ahead by 2 commits"
}
```

### `createBranch`

Create and switch to a new branch.

**Inputs:**
- `type`: Branch type (`feature` | `fix` | `refactor` | `chore`)
- `name`: Branch name (will be slugified)
- `issueId`: Issue ID (optional)
- `fromBranch`: Base branch (default: main)

**Steps:**
1. Validate current state (warn if uncommitted changes)
2. Generate branch name: `{type}/{issueId}-{name}` or `{type}/{name}`
3. Create branch from base
4. Switch to new branch
5. Return confirmation

**Outputs:**
```json
{
  "branch": "feature/PROJ-211-implementation",
  "createdFrom": "main",
  "switched": true
}
```

### `preparePR`

Prepare for pull request.

**Inputs:**
- `title`: PR title (optional, generated from commits if not provided)
- `targetBranch`: Target branch (default: main)

**Steps:**
1. Get current branch state
2. Verify all changes committed
3. Verify branch pushed to remote
4. Generate PR description from commits
5. Return PR preparation summary

**Outputs:**
```json
{
  "ready": true,
  "branch": "feature/PROJ-211-implementation",
  "target": "main",
  "commits": 5,
  "suggestedTitle": "feat(feature): Add implementation",
  "suggestedBody": "## Summary\n...",
  "warnings": []
}
```

### `stash`

Stash current changes.

**Inputs:**
- `message`: Stash message
- `includeUntracked`: Include untracked files (default: false)

**Steps:**
1. Check for changes to stash
2. Create stash with message
3. Return stash reference

**Outputs:**
```json
{
  "stashed": true,
  "index": 0,
  "message": "WIP: Implementation",
  "filesStashed": 5
}
```

### `applyStash`

Apply a stash.

**Inputs:**
- `index`: Stash index (default: 0)
- `drop`: Drop stash after apply (default: false)

**Steps:**
1. Validate stash exists
2. Check for conflicts
3. Apply stash
4. Optionally drop

**Outputs:**
```json
{
  "applied": true,
  "index": 0,
  "conflicts": false,
  "dropped": true
}
```

### `commitChanges`

Prepare and execute commit.

**Inputs:**
- `message`: Commit message (optional, will prompt if not provided)
- `files`: Specific files to commit (optional, default: all staged)
- `type`: Commit type for conventional commits

**Steps:**
1. Stage specified files (or verify staged)
2. Validate commit message format
3. Execute commit
4. Return confirmation

**Outputs:**
```json
{
  "committed": true,
  "hash": "abc1234",
  "message": "feat(feature): Add implementation",
  "filesCommitted": 5
}
```

## Branch Naming Conventions

| Type | Pattern | Example |
|------|---------|---------|
| Feature | `feature/{issue}-{name}` | `feature/PROJ-211-implementation` |
| Bug Fix | `fix/{issue}-{name}` | `fix/PROJ-215-null-check` |
| Refactor | `refactor/{name}` | `refactor/code-cleanup` |
| Chore | `chore/{name}` | `chore/update-dependencies` |

## Commit Message Format

Following conventional commits:
```
{type}({scope}): {description}

{body}

{footer}
```

Types:
- `feat` - New feature
- `fix` - Bug fix
- `refactor` - Code refactoring
- `docs` - Documentation
- `test` - Tests
- `chore` - Maintenance

## Safety Checks

The task includes safety checks:

| Check | Action |
|-------|--------|
| Uncommitted changes before branch switch | Warn, suggest stash |
| Force push | Reject (require explicit override) |
| Push to main/master | Warn, require confirmation |
| Commit with untracked files | Warn |
| Large commit (>20 files) | Warn, suggest splitting |

## Dependencies

None (uses git CLI directly).

## Error Handling

| Error | Action |
|-------|--------|
| Git not available | Fatal error, report |
| Merge conflict | Report files, suggest resolution |
| Push rejected | Report reason, suggest pull first |
| Branch exists | Suggest alternative name |

## Integration with Session

Git state is captured during:
- `StartSession`: Branch, recent commits, stashes
- `EndSession`: Uncommitted changes, push status

Session state includes:
```json
{
  "git": {
    "branch": "feature/PROJ-211",
    "uncommittedAtStart": 0,
    "commitsThisSession": 3
  }
}
```
