"""Moon Astro integration setup."""

from __future__ import annotations

import asyncio
from datetime import timedelta
import importlib

from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant
from homeassistant.exceptions import ConfigEntryNotReady

from .const import CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL, DOMAIN
from .coordinator import MoonAstroCoordinator
from .utils import cleanup_cache_dir, ensure_valid_ephemeris

PLATFORMS: list[str] = ["binary_sensor", "sensor"]

_DATA_REFRESH_TASKS = "refresh_tasks"


def _get_refresh_tasks(hass: HomeAssistant) -> dict[str, asyncio.Task[None]]:
    """Return the internal mapping of entry_id -> initial refresh task.

    Args:
        hass: Home Assistant instance.

    Returns:
        A dictionary mapping config entry IDs to asyncio tasks.
    """
    domain_data = hass.data.setdefault(DOMAIN, {})
    tasks: dict[str, asyncio.Task[None]] = domain_data.setdefault(
        _DATA_REFRESH_TASKS, {}
    )
    return tasks


async def _import_platform(hass: HomeAssistant, platform: str) -> None:
    """Import a platform module in the executor to avoid blocking the event loop.

    Args:
        hass: Home Assistant instance.
        platform: Platform name to import.
    """
    module_path = f"{__package__}.{platform}"
    await hass.async_add_executor_job(importlib.import_module, module_path)


async def async_setup_entry(hass: HomeAssistant, entry: ConfigEntry) -> bool:
    """Set up Moon Astro from a config entry.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.

    Returns:
        True if setup succeeded.

    Raises:
        ConfigEntryNotReady: If required resources cannot be prepared yet.
    """
    # Clean up orphaned temporary files from previous interrupted downloads
    await cleanup_cache_dir(hass)

    # Ensure we have a valid ephemeris file before proceeding
    if not await ensure_valid_ephemeris(hass):
        raise ConfigEntryNotReady("Failed to obtain valid ephemeris file")

    scan_seconds = entry.options.get(CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL)
    coordinator = MoonAstroCoordinator.from_config_entry(
        hass, entry, timedelta(seconds=int(scan_seconds))
    )

    hass.data.setdefault(DOMAIN, {})
    hass.data[DOMAIN][entry.entry_id] = coordinator

    entry.async_on_unload(entry.add_update_listener(async_update_options))

    await hass.config_entries.async_forward_entry_setups(entry, PLATFORMS)

    # Start the first refresh in the background to keep entry setup responsive.
    refresh_task = hass.async_create_task(
        coordinator.async_refresh(),
        name=f"{DOMAIN}-{entry.entry_id}-initial_refresh",
    )
    _get_refresh_tasks(hass)[entry.entry_id] = refresh_task

    return True


async def async_unload_entry(hass: HomeAssistant, entry: ConfigEntry) -> bool:
    """Unload a config entry.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.

    Returns:
        True if unload succeeded.
    """
    unload_ok = await hass.config_entries.async_unload_platforms(entry, PLATFORMS)
    if not unload_ok:
        return False

    tasks = _get_refresh_tasks(hass)
    task = tasks.pop(entry.entry_id, None)
    if task is not None:
        task.cancel()

    domain_data = hass.data.get(DOMAIN, {})
    domain_data.pop(entry.entry_id, None)

    if not tasks:
        domain_data.pop(_DATA_REFRESH_TASKS, None)

    if not domain_data:
        hass.data.pop(DOMAIN, None)

    # Remove ephemeris and cache directory on integration unload.
    await cleanup_cache_dir(
        hass,
        remove_ephemeris=True,
        remove_empty_dir=True,
    )

    return True


async def async_update_options(hass: HomeAssistant, entry: ConfigEntry) -> None:
    """Handle options update by reloading the entry.

    Args:
        hass: Home Assistant instance.
        entry: Configuration entry.
    """
    await hass.config_entries.async_reload(entry.entry_id)
