"""Moon Astro integration setup."""

from __future__ import annotations

from datetime import timedelta
import importlib

from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant

from .const import CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL, DOMAIN
from .coordinator import MoonAstroCoordinator

PLATFORMS: list[str] = ["binary_sensor", "sensor"]


async def _import_platform(hass: HomeAssistant, platform: str) -> None:
    """Import a platform module in the executor to avoid blocking the event loop."""
    module_path = f"{__package__}.{platform}"
    await hass.async_add_executor_job(importlib.import_module, module_path)


async def async_setup_entry(hass: HomeAssistant, entry: ConfigEntry) -> bool:
    """Set up Moon Astro from a config entry."""
    scan_seconds = entry.options.get(CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL)
    coordinator = MoonAstroCoordinator.from_config_entry(
        hass, entry, timedelta(seconds=int(scan_seconds))
    )
    await coordinator.async_config_entry_first_refresh()

    hass.data.setdefault(DOMAIN, {})
    hass.data[DOMAIN][entry.entry_id] = coordinator

    # Pre-import platforms to avoid blocking forward_entry_setups
    for platform in PLATFORMS:
        await _import_platform(hass, platform)

    await hass.config_entries.async_forward_entry_setups(entry, PLATFORMS)
    return True


async def async_unload_entry(hass: HomeAssistant, entry: ConfigEntry) -> bool:
    """Unload a config entry."""
    unload_ok = await hass.config_entries.async_unload_platforms(entry, PLATFORMS)
    if unload_ok:
        hass.data.get(DOMAIN, {}).pop(entry.entry_id, None)
    return unload_ok


async def async_update_options(hass: HomeAssistant, entry: ConfigEntry) -> None:
    """Handle options update by reloading the entry."""
    await hass.config_entries.async_reload(entry.entry_id)
