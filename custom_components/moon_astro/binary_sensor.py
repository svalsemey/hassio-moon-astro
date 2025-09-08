"""Binary sensors for Moon Astro"""

from __future__ import annotations

from homeassistant.components.binary_sensor import BinarySensorEntity
from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant
from homeassistant.helpers.entity import DeviceInfo
from homeassistant.helpers.update_coordinator import CoordinatorEntity
from homeassistant.helpers.entity_platform import AddEntitiesCallback

from .const import DOMAIN, KEY_ABOVE_HORIZON
from .coordinator import MoonAstroCoordinator


async def async_setup_entry(hass: HomeAssistant, entry: ConfigEntry, async_add_entities: AddEntitiesCallback):
    """Set up binary sensor entities."""
    coordinator: MoonAstroCoordinator = hass.data[DOMAIN][entry.entry_id]
    device_info = DeviceInfo(
        identifiers={(DOMAIN, entry.entry_id)},
        manufacturer="Moon Astro",
        model="Skyfield DE440",
        name="Moon Astro",
    )
    async_add_entities([MoonAboveHorizonBinary(coordinator, entry.entry_id, device_info)], True)


class MoonAboveHorizonBinary(CoordinatorEntity[MoonAstroCoordinator], BinarySensorEntity):
    """Binary sensor indicating whether the Moon is above the horizon."""

    def __init__(self, coordinator: MoonAstroCoordinator, entry_id: str, device_info: DeviceInfo) -> None:
        super().__init__(coordinator)
        self._attr_unique_id = f"moon_astro_{entry_id}_above_horizon"
        self._attr_has_entity_name = True
        self._attr_translation_key = "binary_above_horizon"
        self._attr_device_info = device_info
        # Stable slug for entity_id creation
        self._attr_suggested_object_id = "above_horizon"

    @property
    def is_on(self) -> bool:
        data = self.coordinator.data or {}
        return bool(data.get(KEY_ABOVE_HORIZON, False))

    @property
    def icon(self) -> str:
        """Return MDI icon based on state."""
        return "mdi:weather-night" if self.is_on else "mdi:weather-night-off"