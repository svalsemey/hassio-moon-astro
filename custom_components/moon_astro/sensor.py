"""Sensor entities for Moon Astro.

This module defines SensorEntity instances backed by the integration coordinator.
Each sensor reads a single computed value from the coordinator data dictionary.
"""

from __future__ import annotations

from datetime import UTC, datetime
import logging
import math
from typing import Any

from homeassistant.components.sensor import SensorDeviceClass, SensorEntity
from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant
from homeassistant.helpers.entity import DeviceInfo
from homeassistant.helpers.entity_platform import AddEntitiesCallback
from homeassistant.helpers.update_coordinator import CoordinatorEntity

from .const import (
    DOMAIN,
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
    PRECISION_AZIMUTH,
    PRECISION_DISTANCE,
    PRECISION_ECL_GEO,
    PRECISION_ECL_TOPO,
    PRECISION_ELEVATION,
    PRECISION_ILLUMINATION,
    PRECISION_PARALLAX,
    PRECISION_ZODIAC_DEGREE,
)
from .coordinator import MoonAstroCoordinator

_LOGGER = logging.getLogger(__name__)


# key, translation_key, unit, device_class, suggested_display_precision
SENSORS: list[tuple[str, str, str | None, SensorDeviceClass | None, int | None]] = [
    (
        KEY_PHASE,
        "sensor_phase",
        None,
        None,
        None,
    ),
    (
        KEY_AZIMUTH,
        "sensor_azimuth",
        "°",
        None,
        PRECISION_AZIMUTH,
    ),
    (
        KEY_ELEVATION,
        "sensor_elevation",
        "°",
        None,
        PRECISION_ELEVATION,
    ),
    (
        KEY_ILLUMINATION,
        "sensor_illumination",
        "%",
        None,
        PRECISION_ILLUMINATION,
    ),
    (
        KEY_DISTANCE,
        "sensor_distance",
        "km",
        None,
        PRECISION_DISTANCE,
    ),
    (
        KEY_PARALLAX,
        "sensor_parallax",
        "°",
        None,
        PRECISION_PARALLAX,
    ),
    # Topocentric ecliptic coordinates (current)
    (
        KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
        "sensor_ecliptic_longitude_topocentric",
        "°",
        None,
        PRECISION_ECL_TOPO,
    ),
    (
        KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
        "sensor_ecliptic_latitude_topocentric",
        "°",
        None,
        PRECISION_ECL_TOPO,
    ),
    # Geocentric ecliptic coordinates (current)
    (
        KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC,
        "sensor_ecliptic_longitude_geocentric",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    (
        KEY_ECLIPTIC_LATITUDE_GEOCENTRIC,
        "sensor_ecliptic_latitude_geocentric",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    # Time sensors (next)
    (
        KEY_NEXT_RISE,
        "sensor_next_rise",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_NEXT_SET,
        "sensor_next_set",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_NEXT_APOGEE,
        "sensor_next_apogee",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_NEXT_PERIGEE,
        "sensor_next_perigee",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_NEXT_FIRST_QUARTER,
        "sensor_next_first_quarter",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_NEXT_FULL_MOON,
        "sensor_next_full_moon",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_NEXT_LAST_QUARTER,
        "sensor_next_last_quarter",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_NEXT_NEW_MOON,
        "sensor_next_new_moon",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    # Time sensors (previous)
    (
        KEY_PREVIOUS_RISE,
        "sensor_previous_rise",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_PREVIOUS_SET,
        "sensor_previous_set",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_PREVIOUS_APOGEE,
        "sensor_previous_apogee",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_PREVIOUS_PERIGEE,
        "sensor_previous_perigee",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_PREVIOUS_FIRST_QUARTER,
        "sensor_previous_first_quarter",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_PREVIOUS_FULL_MOON,
        "sensor_previous_full_moon",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_PREVIOUS_LAST_QUARTER,
        "sensor_previous_last_quarter",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    (
        KEY_PREVIOUS_NEW_MOON,
        "sensor_previous_new_moon",
        None,
        SensorDeviceClass.TIMESTAMP,
        None,
    ),
    # Full moon names
    (
        KEY_PREVIOUS_FULL_MOON_NAME,
        "sensor_previous_full_moon_name",
        None,
        None,
        None,
    ),
    (
        KEY_PREVIOUS_FULL_MOON_ALT_NAMES,
        "sensor_previous_full_moon_alt_names",
        None,
        None,
        None,
    ),
    (
        KEY_NEXT_FULL_MOON_NAME,
        "sensor_next_full_moon_name",
        None,
        None,
        None,
    ),
    (
        KEY_NEXT_FULL_MOON_ALT_NAMES,
        "sensor_next_full_moon_alt_names",
        None,
        None,
        None,
    ),
    # Ecliptic coordinates at next lunations (geocentric, true-of-date)
    (
        KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON,
        "sensor_ecliptic_longitude_next_full_moon",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    (
        KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON,
        "sensor_ecliptic_latitude_next_full_moon",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    (
        KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON,
        "sensor_ecliptic_longitude_next_new_moon",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    (
        KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON,
        "sensor_ecliptic_latitude_next_new_moon",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    # Ecliptic coordinates at previous lunations (geocentric, true-of-date)
    (
        KEY_ECLIPTIC_LONGITUDE_PREVIOUS_FULL_MOON,
        "sensor_ecliptic_longitude_previous_full_moon",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    (
        KEY_ECLIPTIC_LATITUDE_PREVIOUS_FULL_MOON,
        "sensor_ecliptic_latitude_previous_full_moon",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    (
        KEY_ECLIPTIC_LONGITUDE_PREVIOUS_NEW_MOON,
        "sensor_ecliptic_longitude_previous_new_moon",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    (
        KEY_ECLIPTIC_LATITUDE_PREVIOUS_NEW_MOON,
        "sensor_ecliptic_latitude_previous_new_moon",
        "°",
        None,
        PRECISION_ECL_GEO,
    ),
    # Zodiac sign sensors (string states)
    (
        KEY_ZODIAC_SIGN_CURRENT_MOON,
        "sensor_zodiac_sign_current_moon",
        None,
        None,
        None,
    ),
    (
        KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
        "sensor_zodiac_sign_next_full_moon",
        None,
        None,
        None,
    ),
    (
        KEY_ZODIAC_SIGN_NEXT_NEW_MOON,
        "sensor_zodiac_sign_next_new_moon",
        None,
        None,
        None,
    ),
    (
        KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON,
        "sensor_zodiac_sign_previous_full_moon",
        None,
        None,
        None,
    ),
    (
        KEY_ZODIAC_SIGN_PREVIOUS_NEW_MOON,
        "sensor_zodiac_sign_previous_new_moon",
        None,
        None,
        None,
    ),
    # Zodiac degree sensors (floats 0..30)
    (
        KEY_ZODIAC_DEGREE_CURRENT_MOON,
        "sensor_zodiac_degree_current_moon",
        "°",
        None,
        PRECISION_ZODIAC_DEGREE,
    ),
    (
        KEY_ZODIAC_DEGREE_NEXT_FULL_MOON,
        "sensor_zodiac_degree_next_full_moon",
        "°",
        None,
        PRECISION_ZODIAC_DEGREE,
    ),
    (
        KEY_ZODIAC_DEGREE_NEXT_NEW_MOON,
        "sensor_zodiac_degree_next_new_moon",
        "°",
        None,
        PRECISION_ZODIAC_DEGREE,
    ),
    (
        KEY_ZODIAC_DEGREE_PREVIOUS_NEW_MOON,
        "sensor_zodiac_degree_previous_new_moon",
        "°",
        None,
        PRECISION_ZODIAC_DEGREE,
    ),
    (
        KEY_ZODIAC_DEGREE_PREVIOUS_FULL_MOON,
        "sensor_zodiac_degree_previous_full_moon",
        "°",
        None,
        PRECISION_ZODIAC_DEGREE,
    ),
]

# Stable, non-localized slugs for initial entity_id creation
SUGGESTED_SLUGS: dict[str, str] = {
    KEY_PHASE: "phase",
    KEY_AZIMUTH: "azimuth",
    KEY_ELEVATION: "elevation",
    KEY_ILLUMINATION: "illumination",
    KEY_DISTANCE: "distance",
    KEY_PARALLAX: "parallax",
    KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC: "ecliptic_longitude_topocentric",
    KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC: "ecliptic_latitude_topocentric",
    KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC: "ecliptic_longitude_geocentric",
    KEY_ECLIPTIC_LATITUDE_GEOCENTRIC: "ecliptic_latitude_geocentric",
    KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON: "ecliptic_longitude_next_full_moon",
    KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON: "ecliptic_latitude_next_full_moon",
    KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON: "ecliptic_longitude_next_new_moon",
    KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON: "ecliptic_latitude_next_new_moon",
    KEY_ECLIPTIC_LONGITUDE_PREVIOUS_FULL_MOON: "ecliptic_longitude_previous_full_moon",
    KEY_ECLIPTIC_LATITUDE_PREVIOUS_FULL_MOON: "ecliptic_latitude_previous_full_moon",
    KEY_ECLIPTIC_LONGITUDE_PREVIOUS_NEW_MOON: "ecliptic_longitude_previous_new_moon",
    KEY_ECLIPTIC_LATITUDE_PREVIOUS_NEW_MOON: "ecliptic_latitude_previous_new_moon",
    KEY_NEXT_APOGEE: "next_apogee",
    KEY_NEXT_FIRST_QUARTER: "next_first_quarter",
    KEY_NEXT_FULL_MOON: "next_full_moon",
    KEY_NEXT_FULL_MOON_NAME: "next_full_moon_name",
    KEY_NEXT_FULL_MOON_ALT_NAMES: "next_full_moon_alt_names",
    KEY_NEXT_LAST_QUARTER: "next_last_quarter",
    KEY_NEXT_NEW_MOON: "next_new_moon",
    KEY_NEXT_PERIGEE: "next_perigee",
    KEY_NEXT_RISE: "next_rise",
    KEY_NEXT_SET: "next_set",
    KEY_PREVIOUS_APOGEE: "previous_apogee",
    KEY_PREVIOUS_FIRST_QUARTER: "previous_first_quarter",
    KEY_PREVIOUS_FULL_MOON: "previous_full_moon",
    KEY_PREVIOUS_FULL_MOON_NAME: "previous_full_moon_name",
    KEY_PREVIOUS_FULL_MOON_ALT_NAMES: "previous_full_moon_alt_names",
    KEY_PREVIOUS_LAST_QUARTER: "previous_last_quarter",
    KEY_PREVIOUS_NEW_MOON: "previous_new_moon",
    KEY_PREVIOUS_PERIGEE: "previous_perigee",
    KEY_PREVIOUS_RISE: "previous_rise",
    KEY_PREVIOUS_SET: "previous_set",
    KEY_ZODIAC_SIGN_CURRENT_MOON: "zodiac_sign_current_moon",
    KEY_ZODIAC_SIGN_NEXT_FULL_MOON: "zodiac_sign_next_full_moon",
    KEY_ZODIAC_SIGN_NEXT_NEW_MOON: "zodiac_sign_next_new_moon",
    KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON: "zodiac_sign_previous_full_moon",
    KEY_ZODIAC_SIGN_PREVIOUS_NEW_MOON: "zodiac_sign_previous_new_moon",
    KEY_ZODIAC_DEGREE_CURRENT_MOON: "zodiac_degree_current_moon",
    KEY_ZODIAC_DEGREE_NEXT_NEW_MOON: "zodiac_degree_next_new_moon",
    KEY_ZODIAC_DEGREE_NEXT_FULL_MOON: "zodiac_degree_next_full_moon",
    KEY_ZODIAC_DEGREE_PREVIOUS_FULL_MOON: "zodiac_degree_previous_full_moon",
    KEY_ZODIAC_DEGREE_PREVIOUS_NEW_MOON: "zodiac_degree_previous_new_moon",
}


def _parse_timestamp_to_utc(value: Any) -> datetime | None:
    """Parse an ISO timestamp into a timezone-aware UTC datetime.

    Args:
        value: Raw value returned by the coordinator (expected ISO string or datetime).

    Returns:
        A timezone-aware datetime in UTC, or None if parsing fails.
    """
    if value is None:
        return None

    if isinstance(value, datetime):
        dt = value
    else:
        try:
            dt = datetime.fromisoformat(str(value))
        except (ValueError, TypeError):
            return None

    if dt.tzinfo is None:
        return dt.replace(tzinfo=UTC)
    return dt.astimezone(UTC)


def _values_equal(old: Any, new: Any, *, tol: float | None = None) -> bool:
    """Return True if old and new values should be considered equal for state writing.

    This function supports:
    - exact comparisons for non-numeric types
    - tolerant comparisons for floats when a tolerance is provided

    Args:
        old: Previous value as stored by the entity.
        new: New computed value.
        tol: Absolute tolerance used for float comparisons. If None, floats are compared
            using strict equality.

    Returns:
        True if values are equivalent, False otherwise.
    """
    if old is None and new is None:
        return True
    if old is None or new is None:
        return False

    if tol is not None and isinstance(old, float) and isinstance(new, float):
        # NaN should never be treated as equal to anything to avoid masking issues.
        if math.isnan(old) or math.isnan(new):
            return False
        return abs(old - new) <= tol

    return old == new


async def async_setup_entry(
    hass: HomeAssistant, entry: ConfigEntry, async_add_entities: AddEntitiesCallback
) -> None:
    """Set up sensor entities from a config entry.

    Args:
        hass: Home Assistant instance.
        entry: Config entry for the integration.
        async_add_entities: Callback to add entities.

    Returns:
        None.
    """
    coordinator: MoonAstroCoordinator = hass.data[DOMAIN][entry.entry_id]
    device_info = DeviceInfo(
        identifiers={(DOMAIN, entry.entry_id)},
        manufacturer="Moon Astro",
        model="Skyfield DE440",
        name="Moon Astro",
    )

    entities: list[MoonAstroSensor] = []
    for key, name_key, unit, device_class, precision in SENSORS:
        suggested_slug = SUGGESTED_SLUGS.get(key, key)
        entities.append(
            MoonAstroSensor(
                coordinator=coordinator,
                entry_id=entry.entry_id,
                key=key,
                name_key=name_key,
                unit=unit,
                device_class=device_class,
                suggested_display_precision=precision,
                device_info=device_info,
                suggested_object_id=suggested_slug,
            )
        )

    async_add_entities(entities, update_before_add=False)


class MoonAstroSensor(CoordinatorEntity[MoonAstroCoordinator], SensorEntity):
    """Generic sensor bound to a coordinator value.

    This entity avoids unnecessary state writes by only writing when the computed value
    changes compared to the last written value.
    """

    _attr_has_entity_name = True

    def __init__(
        self,
        coordinator: MoonAstroCoordinator,
        entry_id: str,
        key: str,
        name_key: str,
        unit: str | None,
        device_class: SensorDeviceClass | None,
        suggested_display_precision: int | None,
        device_info: DeviceInfo,
        suggested_object_id: str,
    ) -> None:
        """Initialize the sensor.

        Args:
            coordinator: Integration data coordinator.
            entry_id: Config entry identifier.
            key: Coordinator dictionary key.
            name_key: Translation key used by the frontend.
            unit: Unit of measurement if applicable.
            device_class: Sensor device class if applicable.
            suggested_display_precision: Suggested display precision for UI rendering.
            device_info: Home Assistant device information.
            suggested_object_id: Stable suggested entity_id suffix.

        Returns:
            None.
        """
        super().__init__(coordinator)
        self._key = key
        self._attr_unique_id = f"moon_astro_{entry_id}_{key}"
        self._attr_translation_key = name_key
        self._attr_native_unit_of_measurement = unit
        self._attr_device_class = device_class
        self._attr_device_info = device_info
        self._attr_suggested_display_precision = suggested_display_precision
        self._attr_suggested_object_id = suggested_object_id

        self._last_written_native_value: Any = object()

    @property
    def translation_key(self) -> str | None:
        """Expose explicit translation_key so frontend can translate state values.

        Returns:
            The translation key for this entity, or None.
        """
        return self._attr_translation_key

    def _float_write_tolerance(self) -> float | None:
        """Return the absolute tolerance used to decide whether to write a new float state.

        The goal is to avoid recorder/history spam caused by tiny numeric jitter while
        keeping full precision in the coordinator for downstream computations.

        Strategy:
        - If a suggested display precision is known, use half of the last displayed unit.

        Example:
            - precision=2 -> tol=0.005
            - precision=0 -> tol=0.5
        - If no precision is configured, do not apply a tolerance.

        Returns:
            Absolute tolerance for float comparisons, or None when no tolerance applies.
        """
        prec = self._attr_suggested_display_precision
        if prec is None:
            return None

        try:
            p = int(prec)
        except (TypeError, ValueError):
            return None

        # half of the unit step at the displayed precision
        # p=2 => 10^-2 / 2 => 0.005
        return 0.5 * (10.0 ** (-p))

    def _compute_native_value(self) -> Any:
        """Compute the current native value without triggering a state write.

        Returns:
            The computed native value, already normalized for the entity device class.
        """
        data = self.coordinator.data or {}
        value = data.get(self._key)

        if self.device_class == SensorDeviceClass.TIMESTAMP:
            return _parse_timestamp_to_utc(value)

        return value

    @property
    def native_value(self) -> Any:
        """Return the current state of the sensor.

        Returns:
            The current sensor value in its native type.
        """
        return self._compute_native_value()

    def _handle_coordinator_update(self) -> None:
        """Handle updated data from the coordinator.

        This method avoids updating the entity state if the computed value is unchanged,
        using a tolerant comparison for floats based on the entity display precision.
        """
        new_value = self._compute_native_value()

        tol: float | None = None
        if isinstance(new_value, float) and isinstance(
            self._last_written_native_value, float
        ):
            tol = self._float_write_tolerance()

        if _values_equal(self._last_written_native_value, new_value, tol=tol):
            return

        self._last_written_native_value = new_value
        self.async_write_ha_state()

    def _icon_for_phase(self) -> str | None:
        """Return an icon for the phase sensor, based on the current phase code.

        Returns:
            An MDI icon string if the entity is the phase sensor, otherwise None.
        """
        if self._key != KEY_PHASE:
            return None

        data = self.coordinator.data or {}
        phase = (data.get(KEY_PHASE) or "unknown") or "unknown"

        return {
            "new_moon": "mdi:moon-new",
            "full_moon": "mdi:moon-full",
            "first_quarter": "mdi:moon-first-quarter",
            "last_quarter": "mdi:moon-last-quarter",
            "waning_crescent": "mdi:moon-waning-crescent",
            "waning_gibbous": "mdi:moon-waning-gibbous",
            "waxing_crescent": "mdi:moon-waxing-crescent",
            "waxing_gibbous": "mdi:moon-waxing-gibbous",
        }.get(str(phase))

    def _icon_for_specific_keys(self) -> str | None:
        """Return an icon for keys with a fixed, well-defined icon.

        Returns:
            An MDI icon string if the key matches, otherwise None.
        """
        if self._key in (KEY_NEXT_NEW_MOON, KEY_PREVIOUS_NEW_MOON):
            return "mdi:moon-new"

        if self._key in (KEY_NEXT_FULL_MOON, KEY_PREVIOUS_FULL_MOON):
            return "mdi:moon-full"

        if self._key in (
            KEY_ELEVATION,
            KEY_ZODIAC_DEGREE_CURRENT_MOON,
            KEY_ZODIAC_DEGREE_NEXT_FULL_MOON,
            KEY_ZODIAC_DEGREE_NEXT_NEW_MOON,
            KEY_ZODIAC_DEGREE_PREVIOUS_FULL_MOON,
            KEY_ZODIAC_DEGREE_PREVIOUS_NEW_MOON,
        ):
            return "mdi:angle-acute"

        if self._key == KEY_AZIMUTH:
            return "mdi:angle-obtuse"

        if self._key in (
            KEY_NEXT_FULL_MOON_NAME,
            KEY_NEXT_FULL_MOON_ALT_NAMES,
            KEY_PREVIOUS_FULL_MOON_NAME,
            KEY_PREVIOUS_FULL_MOON_ALT_NAMES,
        ):
            return "mdi:calendar-badge"

        if self._key in (KEY_NEXT_FIRST_QUARTER, KEY_PREVIOUS_FIRST_QUARTER):
            return "mdi:moon-first-quarter"

        if self._key in (KEY_NEXT_LAST_QUARTER, KEY_PREVIOUS_LAST_QUARTER):
            return "mdi:moon-last-quarter"

        if self._key == KEY_ILLUMINATION:
            return "mdi:weather-night"

        if self._key == KEY_DISTANCE:
            return "mdi:ruler"

        if self._key == KEY_PARALLAX:
            return "mdi:circle-multiple"

        if self._key in (KEY_NEXT_RISE, KEY_PREVIOUS_RISE):
            return "mdi:arrow-up-circle"

        if self._key in (KEY_NEXT_SET, KEY_PREVIOUS_SET):
            return "mdi:arrow-down-circle"

        if self._key in (
            KEY_NEXT_APOGEE,
            KEY_NEXT_PERIGEE,
            KEY_PREVIOUS_APOGEE,
            KEY_PREVIOUS_PERIGEE,
        ):
            return "mdi:orbit"

        return None

    def _icon_for_ecliptic_coords(self) -> str | None:
        """Return an icon for ecliptic longitude/latitude sensors.

        Returns:
            An MDI icon string if the key matches, otherwise None.
        """
        if self._key in (
            KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
            KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC,
            KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON,
            KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON,
            KEY_ECLIPTIC_LONGITUDE_PREVIOUS_FULL_MOON,
            KEY_ECLIPTIC_LONGITUDE_PREVIOUS_NEW_MOON,
        ):
            return "mdi:longitude"

        if self._key in (
            KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
            KEY_ECLIPTIC_LATITUDE_GEOCENTRIC,
            KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON,
            KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON,
            KEY_ECLIPTIC_LATITUDE_PREVIOUS_FULL_MOON,
            KEY_ECLIPTIC_LATITUDE_PREVIOUS_NEW_MOON,
        ):
            return "mdi:latitude"

        return None

    def _icon_for_zodiac_sign(self) -> str | None:
        """Return an icon for zodiac sign sensors, based on the current sign value.

        Returns:
            An MDI icon string if the entity is a zodiac sign sensor, otherwise None.
        """
        if self._key not in (
            KEY_ZODIAC_SIGN_CURRENT_MOON,
            KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
            KEY_ZODIAC_SIGN_NEXT_NEW_MOON,
            KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON,
            KEY_ZODIAC_SIGN_PREVIOUS_NEW_MOON,
        ):
            return None

        data = self.coordinator.data or {}
        raw = data.get(self._key)
        sign = str(raw).strip().lower() if raw is not None else ""

        zodiac_icons: dict[str, str] = {
            "aries": "mdi:zodiac-aries",
            "taurus": "mdi:zodiac-taurus",
            "gemini": "mdi:zodiac-gemini",
            "cancer": "mdi:zodiac-cancer",
            "leo": "mdi:zodiac-leo",
            "virgo": "mdi:zodiac-virgo",
            "libra": "mdi:zodiac-libra",
            "scorpio": "mdi:zodiac-scorpio",
            "sagittarius": "mdi:zodiac-sagittarius",
            "capricorn": "mdi:zodiac-capricorn",
            "aquarius": "mdi:zodiac-aquarius",
            "pisces": "mdi:zodiac-pisces",
        }

        return zodiac_icons.get(sign, "mdi:zodiac-aquarius")

    @property
    def icon(self) -> str | None:
        """Return an MDI icon by sensor type.

        Returns:
            An MDI icon string, or None when no icon is defined.
        """
        return (
            self._icon_for_phase()
            or self._icon_for_specific_keys()
            or self._icon_for_ecliptic_coords()
            or self._icon_for_zodiac_sign()
        )
