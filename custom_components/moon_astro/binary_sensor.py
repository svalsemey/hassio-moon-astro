"""Binary sensor entities for Moon Astro."""

from __future__ import annotations

from homeassistant.components.binary_sensor import (
    BinarySensorEntity,
    BinarySensorEntityDescription,
)
from homeassistant.core import HomeAssistant
from homeassistant.helpers.entity_platform import AddConfigEntryEntitiesCallback

from .const import KEY_ABOVE_HORIZON
from .coordinator import MoonAstroConfigEntry
from .entity import MoonAstroEntity

PARALLEL_UPDATES = 0

ABOVE_HORIZON_DESCRIPTION = BinarySensorEntityDescription(
    key=KEY_ABOVE_HORIZON, translation_key=KEY_ABOVE_HORIZON
)


async def async_setup_entry(
    hass: HomeAssistant,
    entry: MoonAstroConfigEntry,
    async_add_entities: AddConfigEntryEntitiesCallback,
) -> None:
    """Set up binary sensor entities.

    Args:
        hass: Home Assistant instance.
        entry: The config entry being set up.
        async_add_entities: Callback used to register entities.
    """
    async_add_entities(
        [
            MoonAstroAboveHorizonBinarySensor(
                entry.runtime_data.coordinator, entry, ABOVE_HORIZON_DESCRIPTION
            )
        ]
    )


class MoonAstroAboveHorizonBinarySensor(MoonAstroEntity, BinarySensorEntity):
    """Binary sensor indicating whether the Moon is above the horizon."""

    @property
    def is_on(self) -> bool | None:
        """Return True above the horizon, False below, or None without data."""
        if (data := self.coordinator.data) is None:
            return None
        return data.get(KEY_ABOVE_HORIZON)
