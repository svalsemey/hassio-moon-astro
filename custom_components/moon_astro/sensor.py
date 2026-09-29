"""Sensor entities for Moon Astro.

Each sensor exposes one key of a coordinator payload. Event-based sensors are bound
to the events coordinator, the others to the main coordinator.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from homeassistant.components.sensor import (
    SensorDeviceClass,
    SensorEntity,
    SensorEntityDescription,
    SensorStateClass,
)
from homeassistant.const import DEGREE, PERCENTAGE, UnitOfLength
from homeassistant.core import HomeAssistant
from homeassistant.helpers.entity_platform import AddConfigEntryEntitiesCallback

from .const import (
    ATTR_NEXT_UPDATE,
    FULL_MOON_NAMES,
    KEY_AZIMUTH,
    KEY_DISTANCE,
    KEY_ECLIPTIC_LATITUDE_GEOCENTRIC,
    KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON,
    KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON,
    KEY_ECLIPTIC_LATITUDE_PREVIOUS_FULL_MOON,
    KEY_ECLIPTIC_LATITUDE_PREVIOUS_NEW_MOON,
    KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
    KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC,
    KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON,
    KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON,
    KEY_ECLIPTIC_LONGITUDE_PREVIOUS_FULL_MOON,
    KEY_ECLIPTIC_LONGITUDE_PREVIOUS_NEW_MOON,
    KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
    KEY_ELEVATION,
    KEY_ILLUMINATION,
    KEY_NEXT_APOGEE,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_FULL_MOON_ALT_NAMES,
    KEY_NEXT_FULL_MOON_NAME,
    KEY_NEXT_LAST_QUARTER,
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_PERIGEE,
    KEY_NEXT_RISE,
    KEY_NEXT_SET,
    KEY_PARALLAX,
    KEY_PHASE,
    KEY_PREVIOUS_APOGEE,
    KEY_PREVIOUS_FIRST_QUARTER,
    KEY_PREVIOUS_FULL_MOON,
    KEY_PREVIOUS_FULL_MOON_ALT_NAMES,
    KEY_PREVIOUS_FULL_MOON_NAME,
    KEY_PREVIOUS_LAST_QUARTER,
    KEY_PREVIOUS_NEW_MOON,
    KEY_PREVIOUS_PERIGEE,
    KEY_PREVIOUS_RISE,
    KEY_PREVIOUS_SET,
    KEY_ZODIAC_DEGREE_CURRENT_MOON,
    KEY_ZODIAC_DEGREE_NEXT_FULL_MOON,
    KEY_ZODIAC_DEGREE_NEXT_NEW_MOON,
    KEY_ZODIAC_DEGREE_PREVIOUS_FULL_MOON,
    KEY_ZODIAC_DEGREE_PREVIOUS_NEW_MOON,
    KEY_ZODIAC_SIGN_CURRENT_MOON,
    KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
    KEY_ZODIAC_SIGN_NEXT_NEW_MOON,
    KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON,
    KEY_ZODIAC_SIGN_PREVIOUS_NEW_MOON,
    PHASE_CODES,
    PRECISION_AZIMUTH,
    PRECISION_DISTANCE,
    PRECISION_ECL_GEO,
    PRECISION_ECL_TOPO,
    PRECISION_ELEVATION,
    PRECISION_ILLUMINATION,
    PRECISION_PARALLAX,
    PRECISION_ZODIAC_DEGREE,
    ZODIAC_SIGNS,
)
from .coordinator import MoonAstroConfigEntry, MoonAstroEventsCoordinator
from .entity import MoonAstroEntity

PARALLEL_UPDATES = 0


@dataclass(frozen=True, kw_only=True)
class MoonAstroSensorDescription(SensorEntityDescription):
    """Describe a Moon Astro sensor.

    Attributes:
        is_event_based: True when the value is produced by the events coordinator and
            only changes at astronomical event boundaries.
    """

    is_event_based: bool = False


def _timestamp(key: str, *, event: bool = False) -> MoonAstroSensorDescription:
    """Return the description of a timestamp sensor."""
    return MoonAstroSensorDescription(
        key=key,
        translation_key=key,
        device_class=SensorDeviceClass.TIMESTAMP,
        is_event_based=event,
    )


def _angle(
    key: str,
    precision: int,
    *,
    event: bool = False,
    enabled: bool = True,
    state_class: SensorStateClass | None = None,
) -> MoonAstroSensorDescription:
    """Return the description of an angular sensor expressed in degrees."""
    return MoonAstroSensorDescription(
        key=key,
        translation_key=key,
        native_unit_of_measurement=DEGREE,
        state_class=state_class,
        suggested_display_precision=precision,
        entity_registry_enabled_default=enabled,
        is_event_based=event,
    )


def _enum(
    key: str, options: Sequence[str], *, event: bool = False
) -> MoonAstroSensorDescription:
    """Return the description of a sensor whose state is one of a fixed set of codes."""
    return MoonAstroSensorDescription(
        key=key,
        translation_key=key,
        device_class=SensorDeviceClass.ENUM,
        options=list(options),
        is_event_based=event,
    )


_FULL_MOON_NAME_OPTIONS: tuple[str, ...] = (*FULL_MOON_NAMES, "blue_moon")

# Detailed coordinates whose human-friendly counterpart exists (zodiac signs, geocentric
# coordinates) are registered disabled; users enable them on demand.
SENSOR_DESCRIPTIONS: tuple[MoonAstroSensorDescription, ...] = (
    # Current position and derived values (main coordinator)
    _enum(KEY_PHASE, PHASE_CODES),
    _angle(KEY_AZIMUTH, PRECISION_AZIMUTH),
    _angle(
        KEY_ELEVATION, PRECISION_ELEVATION, state_class=SensorStateClass.MEASUREMENT
    ),
    MoonAstroSensorDescription(
        key=KEY_ILLUMINATION,
        translation_key=KEY_ILLUMINATION,
        native_unit_of_measurement=PERCENTAGE,
        state_class=SensorStateClass.MEASUREMENT,
        suggested_display_precision=PRECISION_ILLUMINATION,
    ),
    MoonAstroSensorDescription(
        key=KEY_DISTANCE,
        translation_key=KEY_DISTANCE,
        device_class=SensorDeviceClass.DISTANCE,
        native_unit_of_measurement=UnitOfLength.KILOMETERS,
        state_class=SensorStateClass.MEASUREMENT,
        suggested_display_precision=PRECISION_DISTANCE,
    ),
    _angle(KEY_PARALLAX, PRECISION_PARALLAX, state_class=SensorStateClass.MEASUREMENT),
    _angle(KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC, PRECISION_ECL_TOPO, enabled=False),
    _angle(KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC, PRECISION_ECL_TOPO, enabled=False),
    _angle(KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC, PRECISION_ECL_GEO),
    _angle(KEY_ECLIPTIC_LATITUDE_GEOCENTRIC, PRECISION_ECL_GEO),
    _timestamp(KEY_NEXT_RISE),
    _timestamp(KEY_NEXT_SET),
    _timestamp(KEY_PREVIOUS_RISE),
    _timestamp(KEY_PREVIOUS_SET),
    _enum(KEY_ZODIAC_SIGN_CURRENT_MOON, ZODIAC_SIGNS),
    _angle(KEY_ZODIAC_DEGREE_CURRENT_MOON, PRECISION_ZODIAC_DEGREE),
    # Lunation, apsis and derived values (events coordinator)
    _timestamp(KEY_NEXT_NEW_MOON, event=True),
    _timestamp(KEY_NEXT_FIRST_QUARTER, event=True),
    _timestamp(KEY_NEXT_FULL_MOON, event=True),
    _timestamp(KEY_NEXT_LAST_QUARTER, event=True),
    _timestamp(KEY_NEXT_APOGEE, event=True),
    _timestamp(KEY_NEXT_PERIGEE, event=True),
    _timestamp(KEY_PREVIOUS_NEW_MOON, event=True),
    _timestamp(KEY_PREVIOUS_FIRST_QUARTER, event=True),
    _timestamp(KEY_PREVIOUS_FULL_MOON, event=True),
    _timestamp(KEY_PREVIOUS_LAST_QUARTER, event=True),
    _timestamp(KEY_PREVIOUS_APOGEE, event=True),
    _timestamp(KEY_PREVIOUS_PERIGEE, event=True),
    _enum(KEY_NEXT_FULL_MOON_NAME, _FULL_MOON_NAME_OPTIONS, event=True),
    MoonAstroSensorDescription(
        key=KEY_NEXT_FULL_MOON_ALT_NAMES,
        translation_key=KEY_NEXT_FULL_MOON_ALT_NAMES,
        is_event_based=True,
    ),
    _enum(KEY_PREVIOUS_FULL_MOON_NAME, _FULL_MOON_NAME_OPTIONS, event=True),
    MoonAstroSensorDescription(
        key=KEY_PREVIOUS_FULL_MOON_ALT_NAMES,
        translation_key=KEY_PREVIOUS_FULL_MOON_ALT_NAMES,
        is_event_based=True,
    ),
    _angle(
        KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON, PRECISION_ECL_GEO, event=True, enabled=False
    ),
    _angle(
        KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON, PRECISION_ECL_GEO, event=True, enabled=False
    ),
    _angle(
        KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON, PRECISION_ECL_GEO, event=True, enabled=False
    ),
    _angle(
        KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON, PRECISION_ECL_GEO, event=True, enabled=False
    ),
    _angle(
        KEY_ECLIPTIC_LONGITUDE_PREVIOUS_NEW_MOON,
        PRECISION_ECL_GEO,
        event=True,
        enabled=False,
    ),
    _angle(
        KEY_ECLIPTIC_LATITUDE_PREVIOUS_NEW_MOON,
        PRECISION_ECL_GEO,
        event=True,
        enabled=False,
    ),
    _angle(
        KEY_ECLIPTIC_LONGITUDE_PREVIOUS_FULL_MOON,
        PRECISION_ECL_GEO,
        event=True,
        enabled=False,
    ),
    _angle(
        KEY_ECLIPTIC_LATITUDE_PREVIOUS_FULL_MOON,
        PRECISION_ECL_GEO,
        event=True,
        enabled=False,
    ),
    _enum(KEY_ZODIAC_SIGN_NEXT_NEW_MOON, ZODIAC_SIGNS, event=True),
    _enum(KEY_ZODIAC_SIGN_NEXT_FULL_MOON, ZODIAC_SIGNS, event=True),
    _enum(KEY_ZODIAC_SIGN_PREVIOUS_NEW_MOON, ZODIAC_SIGNS, event=True),
    _enum(KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON, ZODIAC_SIGNS, event=True),
    _angle(
        KEY_ZODIAC_DEGREE_NEXT_NEW_MOON, PRECISION_ZODIAC_DEGREE, event=True, enabled=False
    ),
    _angle(
        KEY_ZODIAC_DEGREE_NEXT_FULL_MOON,
        PRECISION_ZODIAC_DEGREE,
        event=True,
        enabled=False,
    ),
    _angle(
        KEY_ZODIAC_DEGREE_PREVIOUS_NEW_MOON,
        PRECISION_ZODIAC_DEGREE,
        event=True,
        enabled=False,
    ),
    _angle(
        KEY_ZODIAC_DEGREE_PREVIOUS_FULL_MOON,
        PRECISION_ZODIAC_DEGREE,
        event=True,
        enabled=False,
    ),
)


async def async_setup_entry(
    hass: HomeAssistant,
    entry: MoonAstroConfigEntry,
    async_add_entities: AddConfigEntryEntitiesCallback,
) -> None:
    """Set up sensor entities from a config entry.

    Args:
        hass: Home Assistant instance.
        entry: Config entry for the integration.
        async_add_entities: Callback to add entities.
    """
    runtime = entry.runtime_data
    async_add_entities(
        MoonAstroSensor(
            runtime.events_coordinator
            if description.is_event_based
            else runtime.coordinator,
            entry,
            description,
        )
        for description in SENSOR_DESCRIPTIONS
    )


class MoonAstroSensor(MoonAstroEntity, SensorEntity):
    """Sensor exposing one key of a coordinator payload."""

    entity_description: MoonAstroSensorDescription

    @property
    def native_value(self) -> Any:
        """Return the payload value of the sensor."""
        return self.coordinator.data.get(self.entity_description.key)

    @property
    def extra_state_attributes(self) -> dict[str, Any] | None:
        """Return the next scheduled refresh time of event-based sensors."""
        if isinstance(self.coordinator, MoonAstroEventsCoordinator):
            return {ATTR_NEXT_UPDATE: self.coordinator.next_refresh_utc}
        return None
