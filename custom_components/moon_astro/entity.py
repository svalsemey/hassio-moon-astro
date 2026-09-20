"""Base entity for Moon Astro."""

from __future__ import annotations

from homeassistant.helpers.device_registry import DeviceEntryType, DeviceInfo
from homeassistant.helpers.entity import EntityDescription
from homeassistant.helpers.update_coordinator import CoordinatorEntity

from .const import DOMAIN, MANUFACTURER, MODEL, NAME
from .coordinator import (
    MoonAstroConfigEntry,
    MoonAstroCoordinator,
    MoonAstroEventsCoordinator,
)


class MoonAstroEntity(
    CoordinatorEntity[MoonAstroCoordinator | MoonAstroEventsCoordinator]
):
    """Entity bound to one of the Moon Astro coordinators.

    All entities of a config entry share a single service device. The unique_id is
    derived from the entry id and the description key, which is also the coordinator
    payload key and the translation key.
    """

    _attr_has_entity_name = True

    def __init__(
        self,
        coordinator: MoonAstroCoordinator | MoonAstroEventsCoordinator,
        entry: MoonAstroConfigEntry,
        description: EntityDescription,
    ) -> None:
        """Initialize the entity.

        Args:
            coordinator: Coordinator providing the entity value.
            entry: Config entry owning the entity.
            description: Entity description; its key is the coordinator payload key.
        """
        super().__init__(coordinator)
        self.entity_description = description
        self._attr_unique_id = f"{DOMAIN}_{entry.entry_id}_{description.key}"
        self._attr_device_info = DeviceInfo(
            identifiers={(DOMAIN, entry.entry_id)},
            entry_type=DeviceEntryType.SERVICE,
            manufacturer=MANUFACTURER,
            model=MODEL,
            name=NAME,
        )
