# AGENTS.md

This repository uses the Codex Workflow System.

## Project Overview

MicEco - Project description

## Workflow Layout

- Repo-local workflow contracts live in `.codex-workflow/`
- Repo-local Codex skills live in `.agents/skills/`
- Planning state lives in `docs/planning/`
- Product and architecture knowledge lives in `knowledge/`

## How To Use The Workflow

Codex does not support Claude-style slash commands. Instead, use the generated
repo-local workflow skills.

Examples:

- "Start a workflow session for this repo"
- "Use the workflow router and tell me what workflow operations are available"
- "Run the workflow implement-feature flow for user profile management"
- "Use workflow-review-code on my current changes"

## Development Guidance

### Code Style

- Follow the established patterns in the codebase
- Keep functions focused and easy to test
- Prefer explicit, readable code over clever shortcuts

### Documentation Rules

- Planning files belong in `docs/planning/`
- Technical specs belong in `docs/specs/`
- Reports belong in `docs/reports/`
- Product requirements belong in `knowledge/prd/`
- Architecture notes belong in `knowledge/architecture/`

### Validation Rules

- Run the smallest relevant test or verification command after edits
- Preserve user changes in existing workflow files unless replacement was requested
- Keep workflow-generated guidance in sync with the actual repository layout

## Environment Notes

### Prerequisites

- Git
- Codex
- Project runtime dependencies

### Installation

```bash
# Project-specific installation commands go here
```

### Development

```bash
# Project-specific dev command goes here
```

### Tests

```bash
# Project-specific test command goes here
```

### Build

```bash
# Project-specific build command goes here
```

## Additional Context
