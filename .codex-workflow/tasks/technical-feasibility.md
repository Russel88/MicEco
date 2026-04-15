---
name: technical-feasibility
description: Technical feasibility assessment and approach comparison
---

# Technical Feasibility Task

Performs comprehensive technical feasibility analysis for enriched PRDs. Evaluates resource requirements, integration complexity, risk factors, and compares implementation approaches to guide go/no-go decisions.

## Operations

### `assess`

Full feasibility analysis for a PRD.

**Steps:**
1. Load enriched PRD with complexity data
2. Evaluate technical feasibility (can we build this with current stack?)
3. Calculate resource requirements (team composition, effort estimates)
4. Score integration complexity (how many services/systems touched)
5. Identify hard prerequisites (what must exist before this can start)
6. Assess risk factors (security review needed, migration complexity, external dependencies)
7. Generate feasibility rating (HIGH/MEDIUM/LOW)
8. Produce recommendation (PROCEED/DEFER/REJECT with conditions)

**Inputs:**
- `prd`: PRD slug or ID (string)
- `team`: Team name or ID for capacity context (string, optional)
- `timeline`: Target timeline constraint (string, optional)

**Outputs:**
```json
{
  "feasibility": "HIGH",
  "resourceEstimate": {
    "roles": [
      { "role": "backend-engineer", "allocation": "100%", "duration": "3 weeks" },
      { "role": "frontend-engineer", "allocation": "50%", "duration": "2 weeks" }
    ],
    "effort": "3 weeks"
  },
  "integrationComplexity": 4,
  "prerequisites": [
    "Redis cluster must be provisioned",
    "SSO provider sandbox account required"
  ],
  "risks": [
    { "severity": "high", "description": "SSO provider API stability unknown" },
    { "severity": "medium", "description": "Session migration requires downtime window" },
    { "severity": "low", "description": "Minor UI component library update needed" }
  ],
  "recommendation": {
    "action": "PROCEED",
    "conditions": [
      "Provision Redis cluster before sprint start",
      "Schedule 30-min downtime window for session migration"
    ]
  }
}
```

### `compareOptions`

Compare multiple implementation approaches.

**Steps:**
1. Define approach options
2. Score each on: effort, risk, maintainability, performance, scalability
3. Generate tradeoff matrix
4. Recommend best option with rationale

**Inputs:**
- `prd`: PRD slug or ID (string)
- `options`: Array of approach descriptions (array of strings)

**Outputs:**
```json
{
  "options": [
    {
      "name": "Option A: Extend existing auth",
      "scores": {
        "effort": 8,
        "risk": 3,
        "maintainability": 7,
        "performance": 6,
        "scalability": 5
      },
      "pros": ["Lower effort", "Familiar codebase", "Minimal migration"],
      "cons": ["Technical debt accumulation", "Limited scalability"]
    },
    {
      "name": "Option B: Full auth rewrite",
      "scores": {
        "effort": 4,
        "risk": 7,
        "maintainability": 9,
        "performance": 9,
        "scalability": 9
      },
      "pros": ["Clean architecture", "Better performance", "Future-proof"],
      "cons": ["Higher effort", "Migration risk", "Longer timeline"]
    }
  ],
  "recommended": "Option A: Extend existing auth",
  "rationale": "Given the 3-week timeline constraint, extending the existing auth system provides the best balance of effort and risk. Option B is the better long-term choice but requires 5+ weeks."
}
```

## Configuration

Feasibility rating thresholds:
- **HIGH**: All prerequisites achievable, risks manageable, team has capacity
- **MEDIUM**: Some prerequisites need work, moderate risks, tight on capacity
- **LOW**: Significant blockers, high risks, insufficient capacity

Scoring dimensions for approach comparison use a 1-10 scale where 10 is best.

## Error Handling

| Error Type | Action |
|------------|--------|
| PRD not enriched | Error with suggestion to run prd-enrichment first |
| Missing complexity data | Run estimateComplexity inline |
| Team data unavailable | Skip capacity analysis, note in output |
| No options provided for compare | Return error requesting at least 2 options |

## Dependencies

- **prd-enrichment**: For enriched PRD data and complexity estimates
- **context-loader**: For codebase and technology stack context
