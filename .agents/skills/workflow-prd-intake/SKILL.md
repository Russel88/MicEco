---
name: workflow-prd-intake
description: Import a PRD from external source (Coda, Confluence, local file) and create Linear Initiative Use when the user asks to prd intake in this repository or refers to the old `/PRDIntake` command.
---

# PRD Intake

This repo-local Codex skill is generated from the installed workflow command contract.
Supporting reference material lives in `.codex-workflow/agents/` and `.codex-workflow/tasks/`.

## Command Contract

# PRDIntake

Import a Product Requirements Document from an external source and set up tracking.

## Purpose

Entry point for bringing a PRD into the workflow system. Fetches content from external sources (Coda, Confluence, local markdown files), parses the document structure, creates or links a Linear Initiative, and stores the PRD locally for downstream processing (validation, enrichment, feasibility).

## Arguments

- `url`: URL to PRD document (Coda, Confluence, or other supported source)
- `file`: Path to local markdown file (alternative to `url`)
- `initiative`: Optional Linear Initiative ID to link to (if omitted, a new Initiative is created)

## Execution

1. Fetch PRD content via MCP
   - For `coda.io` URLs: Use Coda MCP to retrieve document content
   - For Confluence URLs: Use Confluence MCP to retrieve page content
   - For local files: Read markdown file from provided `file` path
   - Normalize content to consistent markdown format

2. Parse PRD structure
   - Identify standard sections: problem statement, target users, requirements, acceptance criteria, constraints
   - Extract metadata: title, author, version, date
   - Flag missing or incomplete sections for later validation

3. Create Initiative in Linear (if `initiative` not provided)
   - Invoke `mcp-sync.createInitiative` to create a new Linear Initiative
   - Set Initiative name from PRD title
   - Link PRD source URL as attachment
   - If `initiative` is provided, verify it exists and link the PRD to it

4. Store PRD reference locally
   - Generate feature slug from PRD title
   - Write parsed PRD to `/docs/planning/prds/{feature-slug}.md`
   - Include frontmatter with source URL, Initiative ID, import timestamp

5. Queue for validation
   - Mark PRD status as `imported` in frontmatter
   - Log intake completion to session state
   - Display next step recommendation: run `workflow-prd-validate`

## Prerequisites

- At least one of `url` or `file` must be provided
- Coda MCP configured and authenticated (for Coda URLs)
- Linear MCP configured and authenticated
- `/docs/planning/prds/` directory exists (created if missing)

## Output Files

- `/docs/planning/prds/{feature-slug}.md` (imported PRD with metadata)
- `/docs/planning/session-state.json` (updated with intake record)

## Example

```
`workflow-prd-intake` url="https://coda.io/d/Project-Docs_d1234/PRD-User-Onboarding_su567"
```

```
`workflow-prd-intake` file="./requirements/checkout-redesign.md" initiative="INI-42"
```

## Related

- `workflow-prd-validate` - Validate the imported PRD for completeness
- `workflow-prd-enrich` - Add technical context to the PRD
- `workflow-prd-feasibility` - Assess implementation feasibility
- `workflow-prd-sequence` - Sequence multiple initiatives
- `workflow-start-session` - Start a session before intake

## Tasks Invoked

- `prd-validation.intake`
- `mcp-sync.createInitiative` (if no `initiative` ID provided)

## Agents Used

- None (orchestration only)
