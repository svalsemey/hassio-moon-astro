"""Moon Astro integration setup.

The shared DE440 kernel is obtained (and downloaded when needed) before the two
coordinators are created and refreshed once, so that every entity starts with data.
"""

from __future__ import annotations

from datetime import timedelta

from homeassistant.const import Platform
from homeassistant.core import HomeAssistant
from homeassistant.exceptions import ConfigEntryNotReady

from .const import (
    CONF_EVENTS_REFRESH_FALLBACK,
    CONF_SCAN_INTERVAL,
    DEFAULT_EVENTS_REFRESH_FALLBACK,
    DEFAULT_SCAN_INTERVAL,
    DOMAIN,
)
from .coordinator import (
    MoonAstroConfigEntry,
    MoonAstroCoordinator,
    MoonAstroEventsCoordinator,
    MoonAstroRuntimeData,
)
from .utils import (
    EphemerisError,
    async_discard_ephemeris,
    async_get_ephemeris,
    async_resolve_time_zone,
)

PLATFORMS: list[Platform] = [Platform.BINARY_SENSOR, Platform.SENSOR]


async def async_setup_entry(hass: HomeAssistant, entry: MoonAstroConfigEntry) -> bool:
    """Set up Moon Astro from a config entry.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.

    Returns:
        True once the platforms are set up.

    Raises:
        ConfigEntryNotReady: If the ephemeris cannot be obtained or a first
            computation fails.
    """
    try:
        eph, ts = await async_get_ephemeris(hass)
    except EphemerisError as err:
        raise ConfigEntryNotReady(
            translation_domain=DOMAIN,
            translation_key="ephemeris_unavailable",
            translation_placeholders={"error": str(err)},
        ) from err

    coordinator = MoonAstroCoordinator(
        hass,
        entry,
        eph=eph,
        ts=ts,
        interval=timedelta(
            seconds=entry.options.get(CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL)
        ),
    )
    events_coordinator = MoonAstroEventsCoordinator(
        hass,
        entry,
        eph=eph,
        ts=ts,
        interval=timedelta(
            seconds=entry.options.get(
                CONF_EVENTS_REFRESH_FALLBACK, DEFAULT_EVENTS_REFRESH_FALLBACK
            )
        ),
        tz=await async_resolve_time_zone(entry),
    )
    await coordinator.async_config_entry_first_refresh()
    await events_coordinator.async_config_entry_first_refresh()

    entry.runtime_data = MoonAstroRuntimeData(
        coordinator=coordinator, events_coordinator=events_coordinator
    )
    entry.async_on_unload(entry.add_update_listener(_async_update_options))
    await hass.config_entries.async_forward_entry_setups(entry, PLATFORMS)
    return True


async def async_unload_entry(hass: HomeAssistant, entry: MoonAstroConfigEntry) -> bool:
    """Unload a config entry.

    Both coordinators are shut down by Home Assistant through the unload callbacks
    they registered on the entry. The ephemeris cache is kept for the next load.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.

    Returns:
        True if the platforms were unloaded.
    """
    return await hass.config_entries.async_unload_platforms(entry, PLATFORMS)


async def async_remove_entry(hass: HomeAssistant, entry: MoonAstroConfigEntry) -> None:
    """Delete the cached ephemeris once the single config entry is removed.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry being removed.
    """
    await async_discard_ephemeris(hass)


async def _async_update_options(
    hass: HomeAssistant, entry: MoonAstroConfigEntry
) -> None:
    """Reload the entry so that new options take effect.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.
    """
    await hass.config_entries.async_reload(entry.entry_id)
