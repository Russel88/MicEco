---
name: async-checkin
description: Automated daily check-in generation from session data
---

# Async Check-in Task

Generates automated daily check-ins by aggregating session data. Supports posting to Slack and storing check-in history for trend analysis.

## Operations

### `generate`

Create daily async check-in from session data.

**Steps:**
1. Scan `/docs/planning/` for yesterday's session summaries
2. Aggregate completed items, in-progress items, blockers
3. Format daily summary in concise format
4. Include cycle progress if active cycle exists
5. Post to Slack if configured (via webhook or MCP)
6. Store check-in in `/docs/planning/checkins/`

**Inputs:**
- `date`: Date to generate check-in for (string, optional - defaults to yesterday)
- `slackChannel`: Slack channel to post to (string, optional)

**Outputs:**
```json
{
  "date": "2026-01-15",
  "completed": [
    "WB-001: Redis session store adapter",
    "WB-003: Tests for session adapter"
  ],
  "inProgress": [
    "WB-002: SSO login endpoint (70% complete)"
  ],
  "blockers": [
    "Waiting on SSO provider sandbox credentials (day 2)"
  ],
  "cycleProgress": "65%",
  "posted": {
    "slack": true
  },
  "checkinPath": "/docs/planning/checkins/2026-01-15.md"
}
```

### `configure`

Set up async check-in preferences.

**Steps:**
1. Set posting schedule (time of day)
2. Set Slack channel or webhook URL
3. Set summary format (brief | detailed)
4. Save configuration to `/docs/planning/checkin-config.json`

**Inputs:**
- `postTime`: Time of day to generate check-in (string, e.g., "09:00")
- `slackChannel`: Slack channel name or webhook URL (string)
- `format`: Summary format (enum: `brief` | `detailed`)

**Outputs:**
```json
{
  "configured": true,
  "configPath": "/docs/planning/checkin-config.json"
}
```

## Configuration

Check-in formats:
- **brief**: Bullet list of completed, in-progress, and blockers (3-5 lines)
- **detailed**: Full item descriptions with time spent and cycle context

Configuration is stored in `/docs/planning/checkin-config.json`:
```json
{
  "postTime": "09:00",
  "slackChannel": "#engineering-standup",
  "format": "brief",
  "includeWeekends": false
}
```

## Error Handling

| Error Type | Action |
|------------|--------|
| No session data for date | Generate empty check-in with "No sessions recorded" note |
| Slack webhook failure | Store check-in locally, retry on next run |
| Invalid date format | Default to yesterday, warn in output |
| Config file missing | Use defaults (brief format, no Slack posting) |

## Dependencies

- **session-management**: For session summaries and progress data
- **cycle-monitoring**: For active cycle progress context
- **mcp-sync**: For Slack posting via webhook or MCP integration
