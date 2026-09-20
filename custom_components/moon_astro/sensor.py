"""Sensor entities for Moon Astro.

Each sensor exposes one key of a coordinator payload. Event-based sensors are bound
to the events coordinator and keep their last known value across restarts until the
first event computation completes.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from homeassistant.components.sensor import (
    RestoreSensor,
    SensorDeviceClass,
    SensorEntityDescription,
    SensorStateClass,
)
from homeassistant.const import DEGREE, PERCENTAGE, UnitOfLength
from homeassistant.core import HomeAssistant, callback
from homeassistant.helpers.entity_platform import AddConfigEntryEntitiesCallback
from homeassistant.util import dt as dt_util

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
from .coordinator import (
    MoonAstroConfigEntry,
    MoonAstroCoordinator,
    MoonAstroEventsCoordinator,
)
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
    state_class: SensorStateClass | None = None,
) -> MoonAstroSensorDescription:
    """Return the description of an angular sensor expressed in degrees."""
    return MoonAstroSensorDescription(
        key=key,
        translation_key=key,
        native_unit_of_measurement=DEGREE,
        state_class=state_class,
        suggested_display_precision=precision,
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
    _angle(KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC, PRECISION_ECL_TOPO),
    _angle(KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC, PRECISION_ECL_TOPO),
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
    _angle(KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON, PRECISION_ECL_GEO, event=True),
    _angle(KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON, PRECISION_ECL_GEO, event=True),
    _angle(KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON, PRECISION_ECL_GEO, event=True),
    _angle(KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON, PRECISION_ECL_GEO, event=True),
    _angle(KEY_ECLIPTIC_LONGITUDE_PREVIOUS_NEW_MOON, PRECISION_ECL_GEO, event=True),
    _angle(KEY_ECLIPTIC_LATITUDE_PREVIOUS_NEW_MOON, PRECISION_ECL_GEO, event=True),
    _angle(KEY_ECLIPTIC_LONGITUDE_PREVIOUS_FULL_MOON, PRECISION_ECL_GEO, event=True),
    _angle(KEY_ECLIPTIC_LATITUDE_PREVIOUS_FULL_MOON, PRECISION_ECL_GEO, event=True),
    _enum(KEY_ZODIAC_SIGN_NEXT_NEW_MOON, ZODIAC_SIGNS, event=True),
    _enum(KEY_ZODIAC_SIGN_NEXT_FULL_MOON, ZODIAC_SIGNS, event=True),
    _enum(KEY_ZODIAC_SIGN_PREVIOUS_NEW_MOON, ZODIAC_SIGNS, event=True),
    _enum(KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON, ZODIAC_SIGNS, event=True),
    _angle(KEY_ZODIAC_DEGREE_NEXT_NEW_MOON, PRECISION_ZODIAC_DEGREE, event=True),
    _angle(KEY_ZODIAC_DEGREE_NEXT_FULL_MOON, PRECISION_ZODIAC_DEGREE, event=True),
    _angle(KEY_ZODIAC_DEGREE_PREVIOUS_NEW_MOON, PRECISION_ZODIAC_DEGREE, event=True),
    _angle(KEY_ZODIAC_DEGREE_PREVIOUS_FULL_MOON, PRECISION_ZODIAC_DEGREE, event=True),
)


def _values_equal(old: Any, new: Any, tolerance: float | None) -> bool:
    """Return True when two native values are equal for state writing purposes.

    Floats are compared with the absolute tolerance when one is given; NaN never
    compares equal, so it is never silently preserved.

    Args:
        old: Last written native value.
        new: Newly computed native value.
        tolerance: Absolute tolerance for float comparisons, or None for equality.

    Returns:
        True if the values are equivalent.
    """
    if tolerance is not None and isinstance(old, float) and isinstance(new, float):
        return abs(old - new) <= tolerance
    return old == new


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


class MoonAstroSensor(MoonAstroEntity, RestoreSensor):
    """Sensor exposing one key of a coordinator payload.

    State writes are skipped while a float value stays within half of the last digit
    shown at the suggested display precision, which limits recorder churn caused by
    numerical jitter. Availability changes are always written.
    """

    entity_description: MoonAstroSensorDescription

    def __init__(
        self,
        coordinator: MoonAstroCoordinator | MoonAstroEventsCoordinator,
        entry: MoonAstroConfigEntry,
        description: MoonAstroSensorDescription,
    ) -> None:
        """Initialize the sensor.

        Args:
            coordinator: Coordinator providing the value.
            entry: Config entry owning the sensor.
            description: Sensor description.
        """
        super().__init__(coordinator, entry, description)
        self._float_tolerance: float | None = (
            None
            if (precision := description.suggested_display_precision) is None
            else 0.5 * 10.0**-precision
        )
        self._last_written_value: Any = None
        self._last_written_available: bool | None = None

    def _compute_native_value(self) -> Any:
        """Return the value to expose without writing the state.

        Event-based sensors keep their last written (or restored) value while the
        events coordinator has not produced the key yet, so a valid restored state
        is never replaced by an unknown one.

        Returns:
            The native value, parsed to an aware UTC datetime for timestamp sensors.
        """
        data = self.coordinator.data
        value = None if data is None else data.get(self.entity_description.key)
        if isinstance(value, str) and self.device_class == SensorDeviceClass.TIMESTAMP:
            parsed = dt_util.parse_datetime(value)
            value = None if parsed is None else dt_util.as_utc(parsed)
        if value is None and self.entity_description.is_event_based:
            return self._last_written_value
        return value

    @property
    def native_value(self) -> Any:
        """Return the state of the sensor."""
        return self._compute_native_value()

    @property
    def extra_state_attributes(self) -> dict[str, Any] | None:
        """Return the next scheduled refresh time of event-based sensors."""
        if isinstance(self.coordinator, MoonAstroEventsCoordinator):
            return {ATTR_NEXT_UPDATE: self.coordinator.next_refresh_utc}
        return None

    async def async_added_to_hass(self) -> None:
        """Restore the last known value of event-based sensors.

        The platform writes the initial state right after this method returns, so
        the restored value is exposed before the first coordinator refresh completes.
        Values outside the ENUM options of a sensor are ignored.
        """
        await super().async_added_to_hass()
        if not self.entity_description.is_event_based:
            return
        if (last := await self.async_get_last_sensor_data()) is None:
            return
        value = last.native_value
        options = self.entity_description.options
        if options is not None and value not in options:
            return
        self._last_written_value = value

    @callback
    def _handle_coordinator_update(self) -> None:
        """Write the state when the value or the availability changed."""
        value = self._compute_native_value()
        if self.available == self._last_written_available and _values_equal(
            self._last_written_value, value, self._float_tolerance
        ):
            return
        self._last_written_value = value
        self._last_written_available = self.available
        self.async_write_ha_state()
