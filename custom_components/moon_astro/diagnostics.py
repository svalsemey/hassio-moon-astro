"""Diagnostics support for Moon Astro.

The config entry diagnostics hold the entry configuration, then for each coordinator
its refresh state, the statistics of its computation times and its latest payload.

The observer location is rounded rather than redacted, so that a computation can be
reproduced without identifying a home. The topocentric angles are rounded alike:
given the time of the computation, which the geocentric values reveal, the precise
altitude and azimuth of the Moon would locate the observer to within meters.
"""

from collections.abc import Mapping
from typing import Any, Final

from homeassistant.core import HomeAssistant

from .const import (
    CONF_ALT,
    CONF_LAT,
    CONF_LON,
    KEY_AZIMUTH,
    KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
    KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
    KEY_ELEVATION,
)
from .coordinator import (
    MoonAstroConfigEntry,
    MoonAstroCoordinator,
    MoonAstroEventsCoordinator,
)

# Decimal places kept for the observer location: 0.1° spans 11 km of latitude and
# at most 11 km of longitude; the elevation is rounded to 100 m.
LOCATION_DECIMALS: Final[Mapping[str, int]] = {CONF_LAT: 1, CONF_LON: 1, CONF_ALT: -2}

# Decimal places kept for the payload angles depending on the observer location:
# 0.1° of Moon altitude or azimuth matches about 10 km on the ground.
TOPOCENTRIC_DECIMALS: Final[Mapping[str, int]] = {
    KEY_AZIMUTH: 1,
    KEY_ELEVATION: 1,
    KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC: 1,
    KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC: 1,
}


def _round_values(
    data: Mapping[str, Any], decimals: Mapping[str, int]
) -> dict[str, Any]:
    """Return a copy of a mapping with the values of the listed keys rounded.

    Args:
        data: Mapping to copy.
        decimals: Decimal places by key; negative values round to tens, hundreds...

    Returns:
        The copy, other values unchanged.
    """
    return {
        key: value if (ndigits := decimals.get(key)) is None else round(value, ndigits)
        for key, value in data.items()
    }


def _milliseconds(seconds: float) -> float:
    """Return a duration in milliseconds rounded to a tenth."""
    return round(seconds * 1000.0, 1)


def _coordinator_diagnostics(
    coordinator: MoonAstroCoordinator | MoonAstroEventsCoordinator,
) -> dict[str, Any]:
    """Describe the refresh state, computation times and payload of a coordinator.

    Args:
        coordinator: Coordinator of the loaded config entry.

    Returns:
        A JSON-serializable dictionary.
    """
    stats = coordinator.computation_stats
    return {
        "last_update_success": coordinator.last_update_success,
        "last_exception": (
            None
            if coordinator.last_exception is None
            else repr(coordinator.last_exception)
        ),
        "update_interval_seconds": (
            None
            if (interval := coordinator.update_interval) is None
            else interval.total_seconds()
        ),
        # The entry is only loaded once both coordinators have computed a payload,
        # so the count is never zero.
        "computation": {
            "reference_time": stats.reference_time,
            "count": stats.count,
            "last_ms": _milliseconds(stats.last),
            "mean_ms": _milliseconds(stats.total / stats.count),
            "min_ms": _milliseconds(stats.minimum),
            "max_ms": _milliseconds(stats.maximum),
        },
        "data": _round_values(coordinator.data, TOPOCENTRIC_DECIMALS),
    }


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
            "data": _round_values(entry.data, LOCATION_DECIMALS),
            "options": dict(entry.options),
        },
        "coordinators": {
            "current": _coordinator_diagnostics(runtime.coordinator),
            "events": {
                **_coordinator_diagnostics(runtime.events_coordinator),
                "next_refresh": runtime.events_coordinator.next_refresh_utc,
            },
        },
    }
