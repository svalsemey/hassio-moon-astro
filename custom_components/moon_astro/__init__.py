"""Moon Astro integration setup.

This module provides Home Assistant entry setup, unload, and removal handlers.
It ensures the Skyfield ephemeris is available at startup, and triggers a
download when missing or invalid. A global lock is used to avoid concurrent
downloads across flows and entry reloads.
"""

from __future__ import annotations

from datetime import datetime, timedelta
import logging
from pathlib import Path
import time

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
    SHARED_EPHEMERIS_KEY,
    MoonAstroConfigEntry,
    MoonAstroCoordinator,
    MoonAstroEventsCoordinator,
    MoonAstroRuntimeData,
    async_get_shared_ephemeris,
)
from .utils import (
    EPHEMERIS_LOCK_KEY,
    async_resolve_time_zone,
    cleanup_cache_dir,
    ensure_valid_ephemeris,
    get_ephemeris_lock,
    get_ephemeris_path,
    validate_ephemeris_file,
)

PLATFORMS: list[Platform] = [Platform.BINARY_SENSOR, Platform.SENSOR]

_LOGGER = logging.getLogger(__name__)


def _safe_stat_size(path: Path) -> int | None:
    """Return file size in bytes if available.

    Args:
        path: File path.

    Returns:
        File size in bytes, or None if stat fails.
    """
    try:
        return path.stat().st_size
    except OSError:
        return None


async def _async_prepare_ephemeris(hass: HomeAssistant, *, reason: str) -> None:
    """Ensure the ephemeris file is present and valid.

    A global lock is used to avoid concurrent downloads. This function emits
    explicit INFO logs describing the detected state and the performed action.

    Args:
        hass: Home Assistant instance.
        reason: A short string describing why the preparation is called.

    Returns:
        None.

    Raises:
        ConfigEntryNotReady: When the ephemeris could not be prepared.
    """
    lock = get_ephemeris_lock(hass)
    eph_path = get_ephemeris_path(hass)
    cache_dir = eph_path.parent

    started = time.monotonic()
    async with lock:
        _LOGGER.info(
            "Ephemeris check started (%s): cache_dir=%s file=%s",
            reason,
            str(cache_dir),
            str(eph_path),
        )

        await cleanup_cache_dir(hass)

        exists = eph_path.exists()
        size = _safe_stat_size(eph_path) if exists else None

        if exists:
            _LOGGER.info(
                "Ephemeris state (%s): present size=%s bytes",
                reason,
                str(size) if size is not None else "unknown",
            )
        else:
            _LOGGER.info("Ephemeris state (%s): missing", reason)

        valid_before = await validate_ephemeris_file(hass, remove_on_invalid=False)
        _LOGGER.info(
            "Ephemeris validation (%s): %s", reason, "ok" if valid_before else "failed"
        )

        if valid_before:
            _LOGGER.info(
                "Ephemeris check completed (%s): no download needed (elapsed=%.3fs)",
                reason,
                time.monotonic() - started,
            )
            return

        _LOGGER.info(
            "Ephemeris download triggered (%s): file was missing or invalid", reason
        )

        dl_started = time.monotonic()
        download_ok = await ensure_valid_ephemeris(hass)
        # The file on disk may have been replaced: drop any kernel loaded from it.
        hass.data.pop(SHARED_EPHEMERIS_KEY, None)
        _LOGGER.info(
            "Ephemeris download completed (%s): %s (elapsed=%.3fs)",
            reason,
            "success" if download_ok else "failed",
            time.monotonic() - dl_started,
        )

        if not download_ok:
            _LOGGER.info(
                "Ephemeris check completed (%s): download failed (elapsed=%.3fs)",
                reason,
                time.monotonic() - started,
            )
            raise ConfigEntryNotReady("Ephemeris file is missing or invalid")

        valid_after = await validate_ephemeris_file(hass, remove_on_invalid=False)
        _LOGGER.info(
            "Ephemeris validation after download (%s): %s",
            reason,
            "ok" if valid_after else "failed",
        )

        if not valid_after:
            _LOGGER.info(
                "Ephemeris check completed (%s): validation failed after download (elapsed=%.3fs)",
                reason,
                time.monotonic() - started,
            )
            raise ConfigEntryNotReady("Ephemeris file is missing or invalid")

        final_size = _safe_stat_size(eph_path)
        _LOGGER.info(
            "Ephemeris check completed (%s): ready size=%s bytes (elapsed=%.3fs)",
            reason,
            str(final_size) if final_size is not None else "unknown",
            time.monotonic() - started,
        )


async def async_setup_entry(hass: HomeAssistant, entry: MoonAstroConfigEntry) -> bool:
    """Set up Moon Astro from a config entry.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.

    Returns:
        True if setup succeeded.

    Raises:
        ConfigEntryNotReady: If required resources cannot be prepared yet.
    """
    await _async_prepare_ephemeris(hass, reason="startup_or_reload")

    try:
        eph, ts = await async_get_shared_ephemeris(hass)
    except (OSError, ValueError, RuntimeError) as err:
        raise ConfigEntryNotReady("Ephemeris file could not be loaded") from err

    time_zone = await async_resolve_time_zone(entry)
    entry.runtime_data = MoonAstroRuntimeData(
        coordinator=MoonAstroCoordinator(
            hass,
            entry,
            eph=eph,
            ts=ts,
            interval=timedelta(
                seconds=int(entry.options.get(CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL))
            ),
            tz=time_zone,
        ),
        events_coordinator=MoonAstroEventsCoordinator(
            hass,
            entry,
            eph=eph,
            ts=ts,
            interval=timedelta(
                seconds=int(
                    entry.options.get(
                        CONF_EVENTS_REFRESH_FALLBACK, DEFAULT_EVENTS_REFRESH_FALLBACK
                    )
                )
            ),
            tz=time_zone,
        ),
    )

    entry.async_on_unload(entry.add_update_listener(async_update_options))
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

    Returns:
        None.
    """
    runtime = entry.runtime_data

    _LOGGER.debug(
        "Initial refresh: starting main coordinator refresh (entry_id=%s)",
        entry.entry_id,
    )
    await runtime.coordinator.async_refresh()
    if not runtime.coordinator.last_update_success:
        _LOGGER.debug(
            "Initial refresh: main coordinator refresh failed (entry_id=%s)",
            entry.entry_id,
        )
        return

    @callback
    def _events_cb(_: datetime) -> None:
        """Start the event-based refresh once the startup delay has elapsed.

        Args:
            _: The trigger time provided by the scheduler.

        Returns:
            None.
        """
        _LOGGER.debug(
            "Initial refresh: running scheduled events refresh (entry_id=%s)",
            entry.entry_id,
        )
        entry.async_create_background_task(
            hass,
            runtime.events_coordinator.async_refresh(),
            name=f"{DOMAIN}-{entry.entry_id}-events_initial_refresh",
        )

    _LOGGER.info(
        "Deferring event-based sensors initial refresh by %s seconds; a periodic fallback refresh is also enabled via options",
        DEFAULT_EVENTS_STARTUP_DELAY,
    )
    entry.async_on_unload(
        async_call_later(hass, DEFAULT_EVENTS_STARTUP_DELAY, _events_cb)
    )


async def async_unload_entry(hass: HomeAssistant, entry: MoonAstroConfigEntry) -> bool:
    """Unload a config entry.

    Both coordinators are shut down by Home Assistant through the unload callbacks
    they registered on the entry, and the pending startup task and timer are
    cancelled the same way. The ephemeris cache is kept to avoid a new download
    after a reload.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.

    Returns:
        True if unload succeeded.
    """
    if not await hass.config_entries.async_unload_platforms(entry, PLATFORMS):
        return False

    await cleanup_cache_dir(hass, remove_empty_dir=True)
    return True


async def async_remove_entry(hass: HomeAssistant, entry: MoonAstroConfigEntry) -> None:
    """Handle config entry removal.

    This function is called when the config entry is removed from Home Assistant.
    It performs definitive cleanup of cached resources, including the ephemeris file,
    and drops the shared in-memory objects once no other entry uses them.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry being removed.

    Returns:
        None.
    """
    _LOGGER.info(
        "Entry removed: deleting ephemeris cache (file=%s)", get_ephemeris_path(hass)
    )
    await cleanup_cache_dir(hass, remove_ephemeris=True, remove_empty_dir=True)

    # The removed entry is still listed at this point, so only other entries count.
    if all(
        other.entry_id == entry.entry_id
        for other in hass.config_entries.async_entries(DOMAIN)
    ):
        hass.data.pop(SHARED_EPHEMERIS_KEY, None)
        hass.data.pop(EPHEMERIS_LOCK_KEY, None)


async def async_update_options(hass: HomeAssistant, entry: MoonAstroConfigEntry) -> None:
    """Handle options update by reloading the entry.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.

    Returns:
        None.
    """
    await hass.config_entries.async_reload(entry.entry_id)
