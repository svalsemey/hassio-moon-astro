"""Moon Astro integration setup.

The shared DE440 kernel is obtained (and downloaded when needed) before the two
coordinators are created; the first refreshes are deferred so that setup returns
quickly even on low-power devices.
"""

from __future__ import annotations

from datetime import datetime, timedelta
import logging

from homeassistant.const import Platform
from homeassistant.core import HomeAssistant, callback
from homeassistant.exceptions import ConfigEntryNotReady
from homeassistant.helpers.event import async_call_later

from .const import (
    CONF_EVENTS_REFRESH_FALLBACK,
    CONF_SCAN_INTERVAL,
    DEFAULT_EVENTS_REFRESH_FALLBACK,
    DEFAULT_EVENTS_STARTUP_DELAY,
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

_LOGGER = logging.getLogger(__name__)


async def async_setup_entry(hass: HomeAssistant, entry: MoonAstroConfigEntry) -> bool:
    """Set up Moon Astro from a config entry.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.

    Returns:
        True once the platforms are set up.

    Raises:
        ConfigEntryNotReady: If the ephemeris cannot be loaded or downloaded yet.
    """
    try:
        eph, ts = await async_get_ephemeris(hass)
    except EphemerisError as err:
        raise ConfigEntryNotReady(
            translation_domain=DOMAIN,
            translation_key="ephemeris_unavailable",
            translation_placeholders={"error": str(err)},
        ) from err

    entry.runtime_data = MoonAstroRuntimeData(
        coordinator=MoonAstroCoordinator(
            hass,
            entry,
            eph=eph,
            ts=ts,
            interval=timedelta(
                seconds=entry.options.get(CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL)
            ),
        ),
        events_coordinator=MoonAstroEventsCoordinator(
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
        ),
    )

    entry.async_on_unload(entry.add_update_listener(_async_update_options))
    await hass.config_entries.async_forward_entry_setups(entry, PLATFORMS)

    entry.async_create_background_task(
        hass,
        _async_deferred_initial_refresh(hass, entry),
        name=f"{DOMAIN}-{entry.entry_id}-initial_refresh",
    )
    return True


async def _async_deferred_initial_refresh(
    hass: HomeAssistant, entry: MoonAstroConfigEntry
) -> None:
    """Run the initial refresh sequence without blocking the setup path.

    The main coordinator is refreshed first to populate frequently changing sensors.
    The events coordinator refresh is then scheduled after a startup delay through
    the Home Assistant scheduler; the pending timer is cancelled automatically if
    the entry is unloaded before it fires.

    Args:
        hass: Home Assistant instance.
        entry: Loaded config entry.
    """
    runtime = entry.runtime_data

    await runtime.coordinator.async_refresh()
    if not runtime.coordinator.last_update_success:
        _LOGGER.debug(
            "Initial refresh: main coordinator refresh failed (entry_id=%s)",
            entry.entry_id,
        )
        return

    @callback
    def _events_cb(_: datetime) -> None:
        """Start the event-based refresh once the startup delay has elapsed."""
        entry.async_create_background_task(
            hass,
            runtime.events_coordinator.async_refresh(),
            name=f"{DOMAIN}-{entry.entry_id}-events_initial_refresh",
        )

    _LOGGER.debug(
        "Initial refresh: event-based refresh deferred by %s seconds (entry_id=%s)",
        DEFAULT_EVENTS_STARTUP_DELAY,
        entry.entry_id,
    )
    entry.async_on_unload(
        async_call_later(hass, DEFAULT_EVENTS_STARTUP_DELAY, _events_cb)
    )


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
