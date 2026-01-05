"""Sensor entities for Moon Astro."""

from __future__ import annotations

from datetime import UTC, datetime
import logging
from typing import Any

from homeassistant.components.sensor import SensorDeviceClass, SensorEntity
from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant
from homeassistant.helpers.entity import DeviceInfo
from homeassistant.helpers.entity_platform import AddEntitiesCallback
from homeassistant.helpers.update_coordinator import CoordinatorEntity

from .const import (
    DOMAIN,
    KEY_AZ,
    KEY_DISTANCE,
    KEY_ECLIPTIC_LATITUDE_GEOCENTRIC,
    KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON,
    KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON,
    KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
    KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC,
    KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON,
    KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON,
    KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
    KEY_EL,
    KEY_ILLUM,
    KEY_NEXT_APOGEE,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_LAST_QUARTER,
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_PERIGEE,
    KEY_NEXT_RISE,
    KEY_NEXT_SET,
    KEY_PARALLAX,
    KEY_PHASE,
    KEY_ZODIAC_DEGREE_CURRENT_MOON,
    KEY_ZODIAC_DEGREE_NEXT_FULL_MOON,
    KEY_ZODIAC_DEGREE_NEXT_NEW_MOON,
    KEY_ZODIAC_SIGN_CURRENT_MOON,
    KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
    KEY_ZODIAC_SIGN_NEXT_NEW_MOON,
    PRECISION_AZ,
    PRECISION_DISTANCE,
    PRECISION_ECL_GEO,
    PRECISION_ECL_TOPO,
    PRECISION_EL,
    PRECISION_ILLUM,
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
        KEY_AZ,
        "sensor_azimuth",
        "°",
        None,
        PRECISION_AZ,
    ),
    (
        KEY_EL,
        "sensor_elevation",
        "°",
        None,
        PRECISION_EL,
    ),
    (
        KEY_ILLUM,
        "sensor_illumination",
        "%",
        None,
        PRECISION_ILLUM,
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
    # Time sensors
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
        KEY_NEXT_NEW_MOON,
        "sensor_next_new_moon",
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
    # Zodiac sign sensors (string states)
    (
        KEY_ZODIAC_SIGN_CURRENT_MOON,
        "sensor_zodiac_sign_current_moon",
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
        KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
        "sensor_zodiac_sign_next_full_moon",
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
        KEY_ZODIAC_DEGREE_NEXT_NEW_MOON,
        "sensor_zodiac_degree_next_new_moon",
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
]

# Stable, non-localized slugs for initial entity_id creation
SUGGESTED_SLUGS: dict[str, str] = {
    KEY_PHASE: "phase",
    KEY_AZ: "azimuth",
    KEY_EL: "elevation",
    KEY_ILLUM: "illumination",
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
    KEY_NEXT_RISE: "next_rise",
    KEY_NEXT_SET: "next_set",
    KEY_NEXT_APOGEE: "next_apogee",
    KEY_NEXT_PERIGEE: "next_perigee",
    KEY_NEXT_NEW_MOON: "next_new_moon",
    KEY_NEXT_FIRST_QUARTER: "next_first_quarter",
    KEY_NEXT_FULL_MOON: "next_full_moon",
    KEY_NEXT_LAST_QUARTER: "next_last_quarter",
    KEY_ZODIAC_SIGN_CURRENT_MOON: "zodiac_sign_current_moon",
    KEY_ZODIAC_SIGN_NEXT_NEW_MOON: "zodiac_sign_next_new_moon",
    KEY_ZODIAC_SIGN_NEXT_FULL_MOON: "zodiac_sign_next_full_moon",
    KEY_ZODIAC_DEGREE_CURRENT_MOON: "zodiac_degree_current_moon",
    KEY_ZODIAC_DEGREE_NEXT_NEW_MOON: "zodiac_degree_next_new_moon",
    KEY_ZODIAC_DEGREE_NEXT_FULL_MOON: "zodiac_degree_next_full_moon",
}


async def async_setup_entry(
    hass: HomeAssistant, entry: ConfigEntry, async_add_entities: AddEntitiesCallback
) -> None:
    """Set up sensor entities from a config entry."""
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
    """Generic sensor bound to a coordinator value."""

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
        """Initialize the sensor."""
        super().__init__(coordinator)
        self._key = key
        self._attr_unique_id = f"moon_astro_{entry_id}_{key}"
        self._attr_translation_key = name_key
        self._attr_native_unit_of_measurement = unit
        self._attr_device_class = device_class
        self._attr_device_info = device_info
        self._attr_suggested_display_precision = suggested_display_precision
        # Stable ASCII object_id slug
        self._attr_suggested_object_id = suggested_object_id

    @property
    def translation_key(self) -> str | None:
        """Expose explicit translation_key so frontend can translate state values."""
        return self._attr_translation_key

    @property
    def native_value(self) -> Any:
        """Return the state of the sensor."""
        data = self.coordinator.data or {}
        value = data.get(self._key)

        # Convert ISO string (local or aware) -> aware UTC datetime for timestamp sensors
        if self.device_class == SensorDeviceClass.TIMESTAMP:
            if value is None:
                return None
            try:
                dt = datetime.fromisoformat(str(value))
            except Exception:  # noqa: BLE001
                return None
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=UTC)
            else:
                dt = dt.astimezone(UTC)
            return dt

        _LOGGER.debug(
            "Sensor %s (%s) raw value: %r", self._key, self._attr_translation_key, value
        )
        return value

    @property
    def icon(self) -> str | None:
        """Return an MDI icon by sensor type."""

        # Phase sensor: dynamic based on phase, illumination, and waxing/waning
        if self._key == KEY_PHASE:
            data = self.coordinator.data or {}
            phase = (data.get(KEY_PHASE) or "unknown") or "unknown"

            if phase == "new_moon":
                return "mdi:moon-new"
            if phase == "full_moon":
                return "mdi:moon-full"
            if phase == "first_quarter":
                return "mdi:moon-first-quarter"
            if phase == "last_quarter":
                return "mdi:moon-last-quarter"
            if phase == "waning_crescent":
                return "mdi:moon-waning-crescent"
            if phase == "waning_gibbous":
                return "mdi:moon-waning-gibbous"
            if phase == "waxing_crescent":
                return "mdi:moon-waxing-crescent"
            if phase == "waxing_gibbous":
                return "mdi:moon-waxing-gibbous"

        # Zodiac degree current moon icon (canonical waning crescent)
        if self._key == KEY_ZODIAC_DEGREE_CURRENT_MOON:
            return "mdi:moon-waning-crescent"

        # Canonical phase-related icons for specific next events
        if self._key in (KEY_NEXT_NEW_MOON, KEY_ZODIAC_DEGREE_NEXT_NEW_MOON):
            return "mdi:moon-new"
        if self._key in (KEY_NEXT_FULL_MOON, KEY_ZODIAC_DEGREE_NEXT_FULL_MOON):
            return "mdi:moon-full"
        if self._key == KEY_NEXT_FIRST_QUARTER:
            return "mdi:moon-first-quarter"
        if self._key == KEY_NEXT_LAST_QUARTER:
            return "mdi:moon-last-quarter"

        # Instant angles and measurements
        if self._key == KEY_AZ:
            return "mdi:angle-obtuse"
        if self._key == KEY_EL:
            return "mdi:angle-acute"
        if self._key == KEY_ILLUM:
            return "mdi:weather-night"
        if self._key == KEY_DISTANCE:
            return "mdi:ruler"
        if self._key == KEY_PARALLAX:
            return "mdi:circle-multiple"

        # Ecliptic longitudes
        if self._key in (
            KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
            KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC,
            KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON,
            KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON,
        ):
            return "mdi:longitude"

        # Ecliptic latitudes
        if self._key in (
            KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
            KEY_ECLIPTIC_LATITUDE_GEOCENTRIC,
            KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON,
            KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON,
        ):
            return "mdi:latitude"

        # Time-based events
        if self._key == KEY_NEXT_RISE:
            return "mdi:arrow-up-circle"
        if self._key == KEY_NEXT_SET:
            return "mdi:arrow-down-circle"
        if self._key in (KEY_NEXT_APOGEE, KEY_NEXT_PERIGEE):
            return "mdi:orbit"

        # Zodiac sensors
        if self._key in (
            KEY_ZODIAC_SIGN_CURRENT_MOON,
            KEY_ZODIAC_SIGN_NEXT_NEW_MOON,
            KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
        ):
            data = self.coordinator.data or {}
            raw = data.get(self._key)
            sign = str(raw).strip().lower() if raw is not None else ""

            zodiac_icons = {
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

        return None
