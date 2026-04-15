---
name: workflow-router
description: Entry point for the repo-local Codex workflow. Use when the user asks to use the workflow generally, wants to know which workflow skills exist in this repository, or refers to an old Claude slash command.
---

# Workflow Router

This repository has a local Codex workflow install.

Use the generated `workflow-*` skills in `.agents/skills/` for concrete operations. Supporting contracts live in `.codex-workflow/commands/`, `.codex-workflow/tasks/`, and `.codex-workflow/agents/`.

## Installed Command Skills

- `workflow-breakdown`: Convert a validated PRD into sized, prioritized Work Items with dependency mapping
- `workflow-check-mcp-status`: Verify MCP server connectivity and compare against project requirements
- `workflow-create-variant`: Interactively create a new workflow variant for a project type
- `workflow-cycle-commit`: Commit a draft Cycle plan to Linear and notify the team
- `workflow-cycle-plan`: Plan a capacity-aware Cycle with scope recommendations and risk assessment
- `workflow-cycle-retro`: Analyze Cycle data for retrospective with estimation accuracy and pattern detection
- `workflow-cycle-status`: Get current Cycle status with progress, blockers, and async check-in data
- `workflow-cycle-summary`: Generate end-of-cycle summary with metrics, demo script, and release notes
- `workflow-end-session`: Close engineering session with handoff notes and cleanup
- `workflow-fix-bug`: Diagnose and fix bugs with structured triage, root cause analysis, and verification
- `workflow-generate-tests`: Generate comprehensive test plans with unit, integration, and edge case coverage
- `workflow-implement-feature`: Implement a feature end-to-end with specification, code generation, and task updates
- `workflow-log-new-feature`: Quickly log a new feature, bug, or task to the project backlog
- `workflow-lookup-feature`: Search existing backlog items by keyword or phrase
- `workflow-prd-enrich`: Add technical context, risk analysis, and complexity estimates to a validated PRD
- `workflow-prd-feasibility`: Generate detailed feasibility analysis with resource estimates and risk assessment for an enriched PRD
- `workflow-prd-intake`: Import a PRD from external source (Coda, Confluence, local file) and create Linear Initiative
- `workflow-prd-sequence`: Recommend implementation order for multiple initiatives based on dependencies, risk, and value
- `workflow-prd-validate`: Check PRD completeness and generate a validation report with scores and recommended actions
- `workflow-refactor-code`: Plan and execute code refactoring with structured approach and safety checks
- `workflow-review-code`: Comprehensive code review for architecture, security, and PRD alignment
- `workflow-run-security-audit`: Execute full security audit for vulnerabilities, dependencies, and configuration
- `workflow-session-status`: Check current session state without modification
- `workflow-setup-project-meta`: Reconfigure project integration metadata for Linear, Coda, and GitHub
- `workflow-setup-workflow`: Interactively set up or update the Claude Workflow System in a project
- `workflow-start-session`: Start an engineering session with context and task loading
- `workflow-submit-variant`: Submit a variant or variant changes to the community via branch and pull request
- `workflow-sync-linear`: Manually sync with task management for issue management
- `workflow-update-docs`: Manage project documentation - audit, sync, and update docs
