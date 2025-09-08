"""Sensor entities for Moon Astro."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

import logging
_LOGGER = logging.getLogger(__name__)

from homeassistant.components.sensor import SensorEntity, SensorDeviceClass
from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant
from homeassistant.helpers.entity_platform import AddEntitiesCallback
from homeassistant.helpers.entity import DeviceInfo
from homeassistant.helpers.update_coordinator import CoordinatorEntity

from .const import (
    DOMAIN,
    KEY_PHASE,
    KEY_AZ,
    KEY_EL,
    KEY_ILLUM,
    KEY_DISTANCE,
    KEY_PARALLAX,
    KEY_ECLIPTIC_LON_TOPO,
    KEY_ECLIPTIC_LAT_TOPO,
    KEY_ECLIPTIC_LON_GEO,
    KEY_ECLIPTIC_LAT_GEO,
    KEY_ECLIPTIC_LON_NEXT_FULL,
    KEY_ECLIPTIC_LAT_NEXT_FULL,
    KEY_ECLIPTIC_LON_NEXT_NEW,
    KEY_ECLIPTIC_LAT_NEXT_NEW,
    KEY_NEXT_RISE,
    KEY_NEXT_SET,
    KEY_NEXT_APOGEE,
    KEY_NEXT_PERIGEE,
    KEY_NEXT_NEW,
    KEY_NEXT_FIRST,
    KEY_NEXT_FULL,
    KEY_NEXT_LAST,
    KEY_WAXING,
    KEY_ZODIAC_SIGN_NEXT_NEW,
    KEY_ZODIAC_SIGN_NEXT_FULL,
    KEY_ZODIAC_DEGREE_NEXT_NEW,
    KEY_ZODIAC_DEGREE_NEXT_FULL,
    PRECISION_AZ,
    PRECISION_EL,
    PRECISION_ILLUM,
    PRECISION_DISTANCE,
    PRECISION_PARALLAX,
    PRECISION_ECL_TOPO,
    PRECISION_ECL_GEO,
    PRECISION_ZODIAC_DEGREE,
)

from .coordinator import MoonAstroCoordinator


SENSORS = [
    # key, translation key, unit, device_class, suggested_display_precision
    (KEY_PHASE, "sensor_phase", None, None, None),

    (KEY_AZ, "sensor_azimuth", "°", None, PRECISION_AZ),
    (KEY_EL, "sensor_elevation", "°", None, PRECISION_EL),
    (KEY_ILLUM, "sensor_illumination", "%", None, PRECISION_ILLUM),
    (KEY_DISTANCE, "sensor_distance", "km", None, PRECISION_DISTANCE),
    (KEY_PARALLAX, "sensor_parallax", "°", None, PRECISION_PARALLAX),

    # Topocentric ecliptic coordinates (renamed)
    (KEY_ECLIPTIC_LON_TOPO, "sensor_ecliptic_longitude_topocentric", "°", None, PRECISION_ECL_TOPO),
    (KEY_ECLIPTIC_LAT_TOPO, "sensor_ecliptic_latitude_topocentric", "°", None, PRECISION_ECL_TOPO),

    # Geocentric ecliptic coordinates (current moment)
    (KEY_ECLIPTIC_LON_GEO, "sensor_ecliptic_longitude_geocentric", "°", None, PRECISION_ECL_GEO),
    (KEY_ECLIPTIC_LAT_GEO, "sensor_ecliptic_latitude_geocentric", "°", None, PRECISION_ECL_GEO),

    # Time sensors
    (KEY_NEXT_RISE, "sensor_next_rise", None, SensorDeviceClass.TIMESTAMP, None),
    (KEY_NEXT_SET, "sensor_next_set", None, SensorDeviceClass.TIMESTAMP, None),
    (KEY_NEXT_APOGEE, "sensor_next_apogee", None, SensorDeviceClass.TIMESTAMP, None),
    (KEY_NEXT_PERIGEE, "sensor_next_perigee", None, SensorDeviceClass.TIMESTAMP, None),
    (KEY_NEXT_NEW, "sensor_next_new_moon", None, SensorDeviceClass.TIMESTAMP, None),
    (KEY_NEXT_FIRST, "sensor_next_first_quarter", None, SensorDeviceClass.TIMESTAMP, None),
    (KEY_NEXT_FULL, "sensor_next_full_moon", None, SensorDeviceClass.TIMESTAMP, None),
    (KEY_NEXT_LAST, "sensor_next_last_quarter", None, SensorDeviceClass.TIMESTAMP, None),

    # Ecliptic coordinates at next lunations (geocentric, true-of-date)
    (KEY_ECLIPTIC_LON_NEXT_FULL, "sensor_ecliptic_longitude_next_full_moon", "°", None, PRECISION_ECL_GEO),
    (KEY_ECLIPTIC_LAT_NEXT_FULL, "sensor_ecliptic_latitude_next_full_moon", "°", None, PRECISION_ECL_GEO),
    (KEY_ECLIPTIC_LON_NEXT_NEW, "sensor_ecliptic_longitude_next_new_moon", "°", None, PRECISION_ECL_GEO),
    (KEY_ECLIPTIC_LAT_NEXT_NEW, "sensor_ecliptic_latitude_next_new_moon", "°", None, PRECISION_ECL_GEO),

    # Zodiac sign sensors (strings)
    (KEY_ZODIAC_SIGN_NEXT_NEW, "sensor_zodiac_sign_next_new_moon", None, None, None),
    (KEY_ZODIAC_SIGN_NEXT_FULL, "sensor_zodiac_sign_next_full_moon", None, None, None),

    # Zodiac degree sensors (floats 0..30)
    (KEY_ZODIAC_DEGREE_NEXT_NEW, "sensor_zodiac_degree_next_new_moon", "°", None, PRECISION_ZODIAC_DEGREE),
    (KEY_ZODIAC_DEGREE_NEXT_FULL, "sensor_zodiac_degree_next_full_moon", "°", None, PRECISION_ZODIAC_DEGREE),
]

# Stable, non-localized slugs for initial entity_id creation
SUGGESTED_SLUGS = {
    KEY_PHASE: "phase",
    KEY_AZ: "azimuth",
    KEY_EL: "elevation",
    KEY_ILLUM: "illumination",
    KEY_DISTANCE: "distance",
    KEY_PARALLAX: "parallax",
    KEY_ECLIPTIC_LON_TOPO: "ecliptic_longitude_topocentric",
    KEY_ECLIPTIC_LAT_TOPO: "ecliptic_latitude_topocentric",
    KEY_ECLIPTIC_LON_GEO: "ecliptic_longitude_geocentric",
    KEY_ECLIPTIC_LAT_GEO: "ecliptic_latitude_geocentric",
    KEY_ECLIPTIC_LON_NEXT_FULL: "ecliptic_longitude_next_full_moon",
    KEY_ECLIPTIC_LAT_NEXT_FULL: "ecliptic_latitude_next_full_moon",
    KEY_ECLIPTIC_LON_NEXT_NEW: "ecliptic_longitude_next_new_moon",
    KEY_ECLIPTIC_LAT_NEXT_NEW: "ecliptic_latitude_next_new_moon",
    KEY_NEXT_RISE: "next_rise",
    KEY_NEXT_SET: "next_set",
    KEY_NEXT_APOGEE: "next_apogee",
    KEY_NEXT_PERIGEE: "next_perigee",
    KEY_NEXT_NEW: "next_new_moon",
    KEY_NEXT_FIRST: "next_first_quarter",
    KEY_NEXT_FULL: "next_full_moon",
    KEY_NEXT_LAST: "next_last_quarter",
    KEY_ZODIAC_SIGN_NEXT_NEW: "zodiac_sign_next_new_moon",
    KEY_ZODIAC_SIGN_NEXT_FULL: "zodiac_sign_next_full_moon",
    KEY_ZODIAC_DEGREE_NEXT_NEW: "zodiac_degree_next_new_moon",
    KEY_ZODIAC_DEGREE_NEXT_FULL: "zodiac_degree_next_full_moon",
}


async def async_setup_entry(hass: HomeAssistant, entry: ConfigEntry, async_add_entities: AddEntitiesCallback):
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
                coordinator, entry.entry_id, key, name_key, unit, device_class, precision, device_info, suggested_slug
            )
        )

    async_add_entities(entities, True)


class MoonAstroSensor(CoordinatorEntity[MoonAstroCoordinator], SensorEntity):
    """Generic sensor bound to a coordinator value."""

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
        super().__init__(coordinator)
        self._key = key
        self._attr_unique_id = f"moon_astro_{entry_id}_{key}"
        self._attr_has_entity_name = True
        self._attr_translation_key = name_key
        self._attr_native_unit_of_measurement = unit
        self._attr_device_class = device_class
        self._attr_device_info = device_info
        self._attr_suggested_display_precision = suggested_display_precision
        # Lock initial entity_id slug (stable, ASCII, non-localized)
        self._attr_suggested_object_id = suggested_object_id

    @property
    def translation_key(self) -> str | None:
        # Expose explicit translation_key so frontend can translate state values
        return self._attr_translation_key

    @property
    def native_value(self) -> Any:
        data = self.coordinator.data or {}
        value = data.get(self._key)

        # Convert ISO string (local) -> aware UTC datetime for timestamp sensors
        if self.device_class == SensorDeviceClass.TIMESTAMP:
            if value is None:
                return None
            try:
                dt = datetime.fromisoformat(value)
            except Exception:
                return None
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            else:
                dt = dt.astimezone(timezone.utc)
            return dt

        # Normalize zodiac sign values to match translation keys
        if self._key in (KEY_ZODIAC_SIGN_NEXT_NEW, KEY_ZODIAC_SIGN_NEXT_FULL):
            if value is None:
                return None
            try:
                value_str = str(value).strip().lower()
            except Exception:
                return None
            aliases = {
                "aries": "aries",
                "taurus": "taurus",
                "gemini": "gemini",
                "cancer": "cancer",
                "leo": "leo",
                "virgo": "virgo",
                "libra": "libra",
                "scorpio": "scorpio",
                "sagittarius": "sagittarius",
                "capricorn": "capricorn",
                "aquarius": "aquarius",
                "pisces": "pisces",
            }
            resolved = aliases.get(value_str, value_str)
            _LOGGER.debug("Zodiac state for %s -> %r (from %r)", self._attr_translation_key, resolved, value)
            return resolved

        _LOGGER.debug("Sensor %s (%s) raw value: %r", self._key, self._attr_translation_key, value)
        return value

    @property
    def icon(self) -> str | None:
        """Return an MDI icon by sensor type.
        - Phase: dynamic (kept, with thresholds)
        - Angles/measures: static icons
        - Ecliptic coordinates: longitude/latitude icon
        - Time events: static time/orbit/calendar icons
        - Zodiac sensors: static degree icon or dynamic sign → MDI mapping
        """

        # Phase sensor: keep dynamic logic based on phase, illumination, and waxing/waning
        if self._key == KEY_PHASE:
            data = self.coordinator.data or {}
            phase = data.get(KEY_PHASE, "unknown") or "unknown"
            illum = float(data.get(KEY_ILLUM, 0.0) or 0.0)
            waxing = bool(data.get(KEY_WAXING, False))

            # Canonical phase icons
            if phase == "new_moon":
                return "mdi:moon-new"
            if phase == "full_moon":
                return "mdi:moon-full"
            if phase == "first_quarter":
                return "mdi:moon-first-quarter"
            if phase == "last_quarter":
                return "mdi:moon-last-quarter"

            # Continuous thresholds with waxing/waning
            if illum < 5.0:
                return "mdi:moon-new"
            if illum <= 45.0:
                return "mdi:moon-waxing-crescent" if waxing else "mdi:moon-waning-crescent"
            if 45.0 < illum < 55.0:
                return "mdi:moon-first-quarter" if waxing else "mdi:moon-last-quarter"
            if illum < 99.0:
                return "mdi:moon-waxing-gibbous" if waxing else "mdi:moon-waning-gibbous"
            return "mdi:moon-full"

        # Canonical phase icons
        if self._key in (KEY_NEXT_NEW, KEY_ZODIAC_DEGREE_NEXT_NEW):
            return "mdi:moon-new"
        if self._key in (KEY_NEXT_FULL, KEY_ZODIAC_DEGREE_NEXT_FULL):
            return "mdi:moon-full"
        if self._key == KEY_NEXT_FIRST:
            return "mdi:moon-first-quarter"
        if self._key == KEY_NEXT_LAST:
            return "mdi:moon-last-quarter"

        # Instant angles and measurements: use simple, stable icons
        if self._key == KEY_AZ:
            # Azimuth (angle along horizon)
            return "mdi:angle-obtuse"
        if self._key == KEY_EL:
            # Elevation (vertical angle)
            return "mdi:angle-acute"
        if self._key == KEY_ILLUM:
            # Moon illumination percentage
            return "mdi:weather-night"
        if self._key == KEY_DISTANCE:
            # Earth–Moon distance
            return "mdi:ruler"
        if self._key == KEY_PARALLAX:
            # Parallax angle
            return "mdi:circle-multiple"

        # Ecliptic longitude (topocentric / geocentric, current or at lunations)
        if self._key in (
            KEY_ECLIPTIC_LON_TOPO,
            KEY_ECLIPTIC_LON_GEO,
            KEY_ECLIPTIC_LON_NEXT_FULL,
            KEY_ECLIPTIC_LON_NEXT_NEW,
        ):
            return "mdi:longitude"

        # Ecliptic latitude (topocentric / geocentric, current or at lunations)
        if self._key in (
            KEY_ECLIPTIC_LAT_TOPO,
            KEY_ECLIPTIC_LAT_GEO,
            KEY_ECLIPTIC_LAT_NEXT_FULL,
            KEY_ECLIPTIC_LAT_NEXT_NEW,
        ):
            return "mdi:latitude"

        # Time-based events (timestamps): use time/orbit/calendar cues
        if self._key == KEY_NEXT_RISE:
            # Next moonrise
            return "mdi:chevron-up-circle"
        if self._key == KEY_NEXT_SET:
            # Next moonset
            return "mdi:chevron-down-circle"
        if self._key in (KEY_NEXT_APOGEE, KEY_NEXT_PERIGEE):
            # Next apogee / perigee (orbital extrema)
            return "mdi:orbit"

        # Zodiac sensors
        if self._key in (KEY_ZODIAC_SIGN_NEXT_NEW, KEY_ZODIAC_SIGN_NEXT_FULL):
            # Dynamic mapping: zodiac sign -> MDI icon
            data = self.coordinator.data or {}
            sign = data.get(self._key)
            sign = (str(sign).strip().lower() if sign is not None else "")

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
            # Fallback to a neutral zodiac icon if sign unknown
            return zodiac_icons.get(sign, "mdi:zodiac-aquarius")

        # Default: no explicit icon (HA may show a generic one)
        return None