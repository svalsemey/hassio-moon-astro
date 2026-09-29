"""Diagnostics support for Moon Astro."""

from __future__ import annotations

from typing import Any

from homeassistant.components.diagnostics import async_redact_data
from homeassistant.core import HomeAssistant

from .const import CONF_LAT, CONF_LON
from .coordinator import MoonAstroConfigEntry

# The observer coordinates identify the user's home.
TO_REDACT = {CONF_LAT, CONF_LON}


async def async_get_config_entry_diagnostics(
    hass: HomeAssistant, entry: MoonAstroConfigEntry
) -> dict[str, Any]:
    """Return diagnostics for a config entry.

    Args:
        hass: Home Assistant instance.
        entry: Loaded config entry.

    Returns:
        A JSON-serializable dictionary describing the entry state.
    """
    runtime = entry.runtime_data
    return {
        "entry": {
            "data": async_redact_data(entry.data, TO_REDACT),
            "options": dict(entry.options),
        },
        "coordinator": {
            "last_update_success": runtime.coordinator.last_update_success,
            "data": runtime.coordinator.data,
        },
        "events_coordinator": {
            "last_update_success": runtime.events_coordinator.last_update_success,
            "next_refresh_utc": runtime.events_coordinator.next_refresh_utc,
            "data": runtime.events_coordinator.data,
        },
    }
