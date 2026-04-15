---
name: workflow-check-mcp-status
description: Verify MCP server connectivity and compare against project requirements Use when the user asks to check mcp status in this repository or refers to the old `/CheckMCPStatus` command.
---

# Check MCP Status

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# CheckMCPStatus

Verify the health and connectivity of all MCP servers required by the project.

## Purpose

Evaluates the actual MCP environment, compares it to expected servers, and reports any issues with connectivity, authentication, or configuration.

## Arguments

- `verbose`: Show detailed connection information (default: false)

## Execution

1. Invoke `mcp-sync` task with action: `verifyConnections`
   - Tests task management MCP connection and auth
   - Tests documentation MCP connection and auth
   - Tests GitHub CLI availability
   - Collects capability schemas

2. Display:
   - Human-readable status summary
   - Connection health for each service
   - Any errors or warnings
   - Recommendations for fixing issues

## Prerequisites

- MCP servers configured in Codex

## Output Format

```
### MCP Status Summary
- Task Management MCP: OK (connected, authenticated)
- Documentation MCP: ERROR (authentication failed)
- GitHub CLI: OK (authenticated)

### Recommendations
- Re-authenticate Documentation MCP: Check API token
```

## Example

```
`workflow-check-mcp-status` verbose=true
```

## Related

- `workflow-start-session` - Runs this check automatically
- `workflow-sync-linear` - Manual task management sync

## Tasks Invoked

- `mcp-sync.verifyConnections`

## Notes

This command is read-only and never modifies MCP configurations. It only reports status and provides recommendations.
