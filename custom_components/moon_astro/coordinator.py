"""Data coordinator for Moon Astro.

This module implements the data coordinator responsible for computing high-precision
Moon position and related ephemerides for Home Assistant.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
import logging
import math
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

import numpy as np
from skyfield import almanac
from skyfield.api import Loader, wgs84
from skyfield.timelib import Time
from timezonefinder import TimezoneFinder

from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant
from homeassistant.helpers.update_coordinator import DataUpdateCoordinator, UpdateFailed

from .const import (
    CACHE_DIR_NAME,
    CONF_ALT,
    CONF_LAT,
    CONF_LON,
    CONF_USE_HA_TZ,
    DARK_MOON,
    DE440_FILE,
    FIRST_QUARTER,
    FULL_MOON,
    FULL_MOON_STRICT_PCT,
    KEY_ABOVE_HORIZON,
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
    KEY_WAXING,
    KEY_ZODIAC_DEGREE_NEXT_FULL_MOON,
    KEY_ZODIAC_DEGREE_NEXT_NEW_MOON,
    KEY_ZODIAC_ICON_NEXT_FULL_MOON,
    KEY_ZODIAC_ICON_NEXT_NEW_MOON,
    KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
    KEY_ZODIAC_SIGN_NEXT_NEW_MOON,
    LAST_QUARTER,
    NEW_MOON_STRICT_PCT,
    QUARTER_TOL_PCT,
)

# Type aliases to improve readability where the 3rd-party library does not expose stable typing.
type Ephemeris = Any
type Timescale = Any
type Apparent = Any

# Centralize recoverable exception sets to avoid broad-except patterns.
_RECOVERABLE_TZ_ERRORS: tuple[type[Exception], ...] = (ZoneInfoNotFoundError,)
_RECOVERABLE_SKYFIELD_ERRORS: tuple[type[Exception], ...] = (
    ValueError,
    RuntimeError,
)
_RECOVERABLE_NUMERIC_ERRORS: tuple[type[Exception], ...] = (
    ArithmeticError,
    FloatingPointError,
    OverflowError,
    ZeroDivisionError,
    ValueError,
    RuntimeError,
)
# Top-level recoverable errors when wrapping update computations into UpdateFailed
_RECOVERABLE_UPDATE_ERRORS: tuple[type[Exception], ...] = (
    OSError,
    ValueError,
    ArithmeticError,
    RuntimeError,
    KeyError,
    TypeError,
)

_LOGGER = logging.getLogger(__name__)


def _detect_timezone(lat: float, lon: float) -> ZoneInfo:
    """Return a best-effort ZoneInfo for given coordinates.

    If timezone cannot be found, UTC is returned.

    Args:
    lat: Latitude in decimal degrees.
    lon: Longitude in decimal degrees.

    Returns:
    A ZoneInfo instance representing the local timezone or UTC as fallback.
    """
    tzname: str | None = None
    tzname = TimezoneFinder().timezone_at(lat=lat, lng=lon)

    try:
        return ZoneInfo(tzname) if tzname else ZoneInfo("UTC")
    except ZoneInfoNotFoundError:
        return ZoneInfo("UTC")


def _tz_for_hass(hass: HomeAssistant) -> ZoneInfo:
    """Return ZoneInfo from Home Assistant configuration with UTC fallback.

    Args:
    hass: Home Assistant instance.

    Returns:
    The configured ZoneInfo or UTC if missing/invalid.
    """
    tzname = hass.config.time_zone
    try:
        return ZoneInfo(tzname) if tzname else ZoneInfo("UTC")
    except _RECOVERABLE_TZ_ERRORS:
        return ZoneInfo("UTC")


def _to_local_iso(dt_utc: datetime | None, tz: ZoneInfo | None) -> str | None:
    """Convert a UTC datetime to ISO 8601 string in provided timezone.

    Args:
    dt_utc: Input datetime expected to be UTC or naive (assumed UTC).
    tz: Target timezone; UTC is used if None.

    Returns:
    The ISO 8601 formatted string or None when input is None.
    """
    if dt_utc is None:
        return None
    tz = tz or ZoneInfo("UTC")
    if dt_utc.tzinfo is None:
        dt_utc = dt_utc.replace(tzinfo=UTC)
    return dt_utc.astimezone(tz).isoformat()


def _safe_time_iso(t_obj: Time | None, tz: ZoneInfo | None) -> str | None:
    """Convert a Skyfield Time to localized ISO 8601 safely.

    Args:
        t_obj: Skyfield Time or None.
        tz: Target timezone; UTC is used if None.

    Returns:
        Localized ISO 8601 string or None.
    """
    if t_obj is None:
        return None

    # Extract the datetime; .utc_datetime() can return ndarray for array Times
    dt_utc_raw = t_obj.utc_datetime()

    # Handle both scalar datetime and ndarray from Skyfield
    if isinstance(dt_utc_raw, datetime):
        # Already a scalar datetime
        dt_utc = dt_utc_raw
    else:
        # Assume it's a numpy array-like object; extract the first element
        try:
            dt_utc = dt_utc_raw.item() if hasattr(dt_utc_raw, "item") else dt_utc_raw[0]
        except (IndexError, TypeError, AttributeError):
            # Fallback: try to convert to datetime
            dt_utc = datetime.fromisoformat(str(dt_utc_raw))

    return _to_local_iso(dt_utc.replace(tzinfo=UTC), tz)


def _moon_illumination_percentage(eph: Ephemeris, t: Time) -> float:
    """Compute Moon illumination in percent at given time.

    This uses the apparent Moon from Earth's geocenter and the ICRF position
    method. The Sun must be provided as a vector function (e.g., eph["sun"]),
    not as an Apparent/Astrometric, so the method can evaluate it at the same
    instant as the Moon position.

    Args:
        eph: Loaded ephemeris.
        t: Skyfield Time.

    Returns:
        Moon illuminated fraction in percent (0..100).
    """
    earth = eph["earth"]
    moon_app = earth.at(t).observe(eph["moon"]).apparent()
    frac = moon_app.fraction_illuminated(eph["sun"])
    return float(frac * 100.0)


def _moon_phase_name(eph: Ephemeris, t: Time, ts: Timescale | None = None) -> str:
    """Return a readable moon phase name using events proximity and illumination.

    Strategy:
    - Compute illumination at t and t+6h to infer waxing.
    - Query discrete moon phase events around t (±30 days).
    - If the nearest event is within a small window, return precise labels.
    - Otherwise, fallback to a generic label based on illumination and waxing.

    Args:
        eph: Loaded ephemeris.
        t: Skyfield Time for which the phase name is requested.
        ts: Skyfield Timescale used to compute t+6h (optional).

    Returns:
        A phase name among: new_moon, waxing_crescent, first_quarter, waxing_gibbous,
        full_moon, waning_gibbous, last_quarter, waning_crescent.
    """
    illum_now = _moon_illumination_percentage(eph, t)

    if ts is not None:
        t_future = ts.tt_jd(t.tt + 6.0 / 24.0)

        def _illum_future(ttime: Time) -> float:
            """Return the illumination percentage at a future time."""
            return _moon_illumination_percentage(eph, ttime)

        illum_future = _illum_future(t_future)
        waxing = illum_future > illum_now + 1e-6
    else:
        waxing = True

    f = almanac.moon_phases(eph)
    t0 = t - 30.0
    t1 = t + 30.0
    times, phases = almanac.find_discrete(t0, t1, f)

    # Identify surrounding phase events
    last_ev: tuple[Time, int] | None = None
    next_ev: tuple[Time, int] | None = None
    for ti, pv in zip(times, phases, strict=False):
        if ti.tt <= t.tt:
            last_ev = (ti, int(pv))
        elif next_ev is None and ti.tt > t.tt:
            next_ev = (ti, int(pv))

    def close_to(pct: float, target: float, tol: float = QUARTER_TOL_PCT) -> bool:
        """Return True when pct is within ±tol of target."""
        return abs(pct - target) <= tol

    # Windows (hours) for "exact" labels
    WIN_FULL_H = 6.0
    WIN_QUARTER_H = 6.0
    WIN_NEW_H = 6.0

    # Helper to test time proximity
    def _within_hours(t_event: Time | None, hours_window: float) -> bool:
        """Return True if t_event is within hours_window of t."""
        if t_event is None:
            return False
        return abs((t_event.tt - t.tt) * 24.0) <= hours_window

    # Find nearest known event types around t
    last_time = last_ev[0] if last_ev else None
    next_time = next_ev[0] if next_ev else None
    last_type = last_ev[1] if last_ev else None
    next_type = next_ev[1] if next_ev else None

    # Exact event labeling using proximity windows
    if (last_type == DARK_MOON and _within_hours(last_time, WIN_NEW_H)) or (
        next_type == DARK_MOON and _within_hours(next_time, WIN_NEW_H)
    ):
        return "new_moon"
    if (last_type == FIRST_QUARTER and _within_hours(last_time, WIN_QUARTER_H)) or (
        next_type == FIRST_QUARTER and _within_hours(next_time, WIN_QUARTER_H)
    ):
        return "first_quarter" if waxing else "last_quarter"
    if (last_type == FULL_MOON and _within_hours(last_time, WIN_FULL_H)) or (
        next_type == FULL_MOON and _within_hours(next_time, WIN_FULL_H)
    ):
        return "full_moon"
    if (last_type == LAST_QUARTER and _within_hours(last_time, WIN_QUARTER_H)) or (
        next_type == LAST_QUARTER and _within_hours(next_time, WIN_QUARTER_H)
    ):
        return "last_quarter" if not waxing else "first_quarter"

    # Hints based on last canonical event and illumination
    if last_ev is not None:
        _, phase_value = last_ev
        if phase_value == DARK_MOON:
            if waxing and illum_now <= NEW_MOON_STRICT_PCT:
                return "new_moon"
        elif phase_value == FIRST_QUARTER:
            if close_to(illum_now, 50.0) and waxing:
                return "first_quarter"
        elif phase_value == LAST_QUARTER:
            if close_to(illum_now, 50.0) and (not waxing):
                return "last_quarter"

    # Generic mapping by illumination and waxing
    if illum_now <= NEW_MOON_STRICT_PCT:
        return "new_moon" if waxing else "waning_crescent"
    if illum_now <= 45.0:
        return "waxing_crescent" if waxing else "waning_crescent"
    if 45.0 < illum_now < 55.0:
        return "first_quarter" if waxing else "last_quarter"
    if illum_now < FULL_MOON_STRICT_PCT:
        return "waxing_gibbous" if waxing else "waning_gibbous"
    # Close to 100% but not in the proximity window: choose gibbous based on waxing state
    return "waxing_gibbous" if waxing else "waning_gibbous"


def _topocentric_vectors(
    eph: Ephemeris, t: Time, lat: float, lon: float, alt_m: float
) -> tuple[Apparent, float, float, float]:
    """Compute topocentric apparent vector and basic quantities.

    Args:
        eph: Loaded ephemeris.
        t: Skyfield Time.
        lat: Latitude in degrees.
        lon: Longitude in degrees.
        alt_m: Elevation in meters.

    Returns:
        A tuple of (apparent vector, azimuth in degrees, elevation in degrees, distance in km).
    """
    earth = eph["earth"]
    observer = wgs84.latlon(
        latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt_m
    )
    ap = earth + observer
    apparent = ap.at(t).observe(eph["moon"]).apparent()
    alt, az, distance = apparent.altaz()
    return apparent, float(az.degrees), float(alt.degrees), float(distance.km)


def _geocentric_vector(eph: Ephemeris, t: Time) -> Apparent:
    """Return geocentric apparent vector of the Moon at time t.

    Args:
        eph: Loaded ephemeris.
        t: Skyfield Time.

    Returns:
        Apparent vector seen from geocenter.
    """
    earth = eph["earth"]
    return earth.at(t).observe(eph["moon"]).apparent()


# ---------- High-accuracy nutation (IAU 1980: full 106-term) and ecliptic-of-date ----------


def _julian_centuries_TT_from_tt(tt: float) -> float:
    """Convert TT Julian Date to Julian centuries since J2000.0.

    Args:
        tt: Julian Date (Terrestrial Time).

    Returns:
        Julian centuries from J2000.0.
    """
    return (tt - 2451545.0) / 36525.0


def _deg_to_rad(x: float) -> float:
    """Convert degrees to radians.

    Args:
        x: Angle in degrees.

    Returns:
        Angle in radians.
    """
    return x * math.pi / 180.0


def _arcsec_to_rad(x: float) -> float:
    """Convert arcseconds to radians.

    Args:
        x: Angle in arcseconds.

    Returns:
        Angle in radians.
    """
    return _deg_to_rad(x / 3600.0)


def _mean_obliquity_arcsec(T: float) -> float:
    """Compute mean obliquity (IAU 2006) in arcseconds.

    IAU 2006 polynomial for mean obliquity; accurate for centuries near J2000

    Args:
        T: Julian centuries since J2000.0.

    Returns:
        Mean obliquity in arcseconds.
    """
    U = T / 100.0
    return (
        84381.406
        - 4680.93 * U
        - 1.55 * U**2
        + 1999.25 * U**3
        - 51.38 * U**4
        - 249.67 * U**5
        - 39.05 * U**6
        + 7.12 * U**7
        + 27.87 * U**8
        + 5.79 * U**9
        + 2.45 * U**10
    )


def _fundamental_arguments_deg(T: float) -> tuple[float, float, float, float, float]:
    """Return fundamental Delaunay arguments (degrees) for IAU 1980 nutation.

    Args:
        T: Julian centuries since J2000.0.

    Returns:
        Tuple of (M' lunar anomaly, M solar anomaly, F, D, Ω) in degrees.
    """
    # Mean anomaly of the Moon (M'), of the Sun (M), Moon's argument of latitude (F),
    # Moon's elongation from the Sun (D), and longitude of the ascending node (Ω).
    Lm = (
        134.96298139
        + (1325.0 * 360.0 + 198.8673981) * T
        + 0.0086972 * T**2
        + T**3 / 56250.0
    )
    Ls = (
        357.52772333
        + (99.0 * 360.0 + 359.0503400) * T
        - 0.0001603 * T**2
        - T**3 / 300000.0
    )
    F = (
        93.27191028
        + (1342.0 * 360.0 + 82.0175381) * T
        - 0.0036825 * T**2
        + T**3 / 327270.0
    )
    D = (
        297.85036306
        + (1236.0 * 360.0 + 307.1114800) * T
        - 0.0019142 * T**2
        + T**3 / 189474.0
    )
    Om = (
        125.04452222
        - (5.0 * 360.0 + 134.1362608) * T
        + 0.0020708 * T**2
        + T**3 / 450000.0
    )

    def norm(x: float) -> float:
        """Normalize an angle to [0, 360)."""
        return (x % 360.0 + 360.0) % 360.0

    return norm(Lm), norm(Ls), norm(F), norm(D), norm(Om)


# Full IAU 1980 106-term nutation series
# Format: (D, M, M', F, Ω, Δψ_arcsec, Δψ_t_arcsec_per_century, Δε_arcsec, Δε_t_arcsec_per_century)
_IAU1980_TERMS: list[tuple[int, int, int, int, int, float, float, float, float]] = [
    (0, 0, 0, 0, 1, -171996.0, -174.2, 92025.0, 8.9),
    (0, 0, 2, -2, 2, -13187.0, -1.6, 5736.0, -3.1),
    (0, 0, 2, 0, 2, -2274.0, -0.2, 977.0, -0.5),
    (0, 0, 0, 0, 2, 2062.0, 0.2, -895.0, 0.5),
    (0, 1, 0, 0, 0, 1426.0, -3.4, 54.0, -0.1),
    (1, 0, 0, 0, 0, 712.0, 0.1, -7.0, 0.0),
    (0, 1, 2, -2, 2, -517.0, 1.2, 224.0, -0.6),
    (0, 0, 2, 0, 1, -386.0, -0.4, 200.0, 0.0),
    (1, 0, 2, 0, 2, -301.0, 0.0, 129.0, -0.1),
    (0, -1, 2, -2, 2, 217.0, -0.5, -95.0, 0.3),
    (1, 0, 0, -2, 0, -158.0, 0.0, 0.0, 0.0),
    (0, 0, 2, -2, 1, 129.0, 0.1, -70.0, 0.0),
    (-1, 0, 2, 0, 2, 123.0, 0.0, -53.0, 0.0),
    (0, 0, 0, 2, 0, 63.0, 0.0, 0.0, 0.0),
    (1, 0, 2, -2, 2, 63.0, 0.1, -33.0, 0.0),
    (-1, 0, 0, 2, 0, -58.0, -0.1, 0.0, 0.0),
    (-1, 0, 2, 2, 2, -51.0, 0.0, 27.0, 0.0),
    (1, 0, 2, 0, 1, 48.0, 0.0, -24.0, 0.0),
    (0, 0, 2, 2, 2, -38.0, 0.0, 16.0, 0.0),
    (2, 0, 2, 0, 2, -31.0, 0.0, 13.0, 0.0),
    (2, 0, 0, 0, 0, 29.0, 0.0, 0.0, 0.0),
    (0, 0, 2, 0, 0, 29.0, 0.0, 0.0, 0.0),
    (0, 0, 2, -2, 0, 26.0, 0.0, 0.0, 0.0),
    (-1, 0, 2, 0, 1, 21.0, 0.0, -10.0, 0.0),
    (0, 2, 0, 0, 0, -16.0, 0.0, 0.0, 0.0),
    (1, 0, 0, 0, 1, 16.0, 0.0, -8.0, 0.0),
    (0, 0, 0, 0, 3, -15.0, 0.0, 9.0, 0.0),
    (1, 0, 2, -2, 1, -13.0, 0.0, 7.0, 0.0),
    (0, 1, 0, 0, 1, -12.0, 0.0, 6.0, 0.0),
    (-1, 0, 0, 0, 1, 11.0, 0.0, -5.0, 0.0),
    (0, 1, 2, 0, 2, -10.0, 0.0, 5.0, 0.0),
    (0, -1, 2, 0, 2, -8.0, 0.0, 3.0, 0.0),
    (2, 0, 2, -2, 2, -7.0, 0.0, 3.0, 0.0),
    (1, 1, 0, 0, 0, -7.0, 0.0, 0.0, 0.0),
    (-1, 1, 0, 0, 0, -7.0, 0.0, 0.0, 0.0),
    (0, 1, 2, -2, 2, -7.0, 0.0, 3.0, 0.0),
    (0, 0, 0, 2, 1, 6.0, 0.0, -3.0, 0.0),
    (1, 0, 2, 2, 2, -6.0, 0.0, 3.0, 0.0),
    (1, 0, 0, 2, 0, 6.0, 0.0, 0.0, 0.0),
    (2, 0, 2, 0, 1, -6.0, 0.0, 3.0, 0.0),
    (0, 0, 0, 2, 0, -5.0, 0.0, 0.0, 0.0),
    (0, -1, 2, -2, 2, 5.0, 0.0, -3.0, 0.0),
    (2, 0, 2, -2, 1, 5.0, 0.0, -3.0, 0.0),
    (0, 1, 0, 0, 2, -5.0, 0.0, 3.0, 0.0),
    (1, 0, 2, -2, 0, -4.0, 0.0, 0.0, 0.0),
    (0, 0, 0, 1, 0, -4.0, 0.0, 0.0, 0.0),
    (1, 1, 0, 0, 1, -4.0, 0.0, 2.0, 0.0),
    (1, -1, 0, 0, 1, -4.0, 0.0, 2.0, 0.0),
    (1, 0, 0, -1, 0, -4.0, 0.0, 0.0, 0.0),
    (0, 0, 2, 1, 2, -4.0, 0.0, 2.0, 0.0),
    (1, 0, 0, 1, 0, 3.0, 0.0, 0.0, 0.0),
    (1, -1, 2, 0, 2, -3.0, 0.0, 1.0, 0.0),
    (0, -1, 2, 0, 1, -3.0, 0.0, 1.0, 0.0),
    (1, 1, 2, 0, 2, -3.0, 0.0, 1.0, 0.0),
    (-1, 1, 2, 0, 2, -3.0, 0.0, 1.0, 0.0),
    (3, 0, 2, 0, 2, -3.0, 0.0, 1.0, 0.0),
    (0, 0, 0, 0, 1, 3.0, 0.0, -1.0, 0.0),
    (-1, 0, 2, 2, 1, -3.0, 0.0, 1.0, 0.0),
    (0, 0, 2, 2, 1, -3.0, 0.0, 1.0, 0.0),
    (1, 0, 2, 2, 1, -3.0, 0.0, 1.0, 0.0),
    (-1, 0, 2, -2, 1, -2.0, 0.0, 1.0, 0.0),
    (2, 0, 0, 0, 1, 2.0, 0.0, -1.0, 0.0),
    (1, 0, 0, 0, 2, -2.0, 0.0, 1.0, 0.0),
    (2, 0, 2, -2, 2, -2.0, 0.0, 1.0, 0.0),
    (0, -1, 2, -2, 1, -2.0, 0.0, 1.0, 0.0),
    (0, 1, 2, -2, 1, -2.0, 0.0, 1.0, 0.0),
    (0, -2, 0, 2, 0, -2.0, 0.0, 0.0, 0.0),
    (2, 0, 0, -2, 1, 2.0, 0.0, -1.0, 0.0),
    (-2, 0, 2, 0, 1, 2.0, 0.0, -1.0, 0.0),
    (0, 0, 2, 0, 1, 2.0, 0.0, -1.0, 0.0),
    (2, 0, 2, 0, 1, 2.0, 0.0, -1.0, 0.0),
    (0, 0, 0, 2, 1, 2.0, 0.0, -1.0, 0.0),
    (1, 0, 2, -2, 2, -1.0, 0.0, 0.0, 0.0),
    (1, 0, 0, 0, 0, -1.0, 0.0, 0.0, 0.0),
    (-1, 0, 0, 0, 2, 1.0, 0.0, 0.0, 0.0),
    (1, 0, 0, -2, 1, 1.0, 0.0, 0.0, 0.0),
    (0, 0, 0, 2, 2, -1.0, 0.0, 0.0, 0.0),
    (0, 0, 2, 2, 2, -1.0, 0.0, 0.0, 0.0),
    (1, 0, 2, 0, 2, -1.0, 0.0, 0.0, 0.0),
    (0, 0, 2, 0, 2, -1.0, 0.0, 0.0, 0.0),
    (1, 0, 0, 2, 0, 1.0, 0.0, 0.0, 0.0),
    (0, 0, 0, 2, 1, -1.0, 0.0, 0.0, 0.0),
    (1, 0, 2, -2, 1, -1.0, 0.0, 0.0, 0.0),
    (1, 1, 0, -2, 0, -1.0, 0.0, 0.0, 0.0),
    (1, -1, 0, -2, 0, -1.0, 0.0, 0.0, 0.0),
    (2, 0, 0, 0, 0, -1.0, 0.0, 0.0, 0.0),
    (0, 1, 2, 0, 1, -1.0, 0.0, 0.0, 0.0),
    (-1, 0, 2, 2, 2, -1.0, 0.0, 0.0, 0.0),
    (0, -1, 2, 2, 2, -1.0, 0.0, 0.0, 0.0),
    (1, -1, 2, 0, 1, -1.0, 0.0, 0.0, 0.0),
    (0, 0, 2, -1, 2, -1.0, 0.0, 0.0, 0.0),
    (1, 0, 0, 0, 1, -1.0, 0.0, 0.0, 0.0),
    (1, 0, 0, -1, 1, -1.0, 0.0, 0.0, 0.0),
    (0, 1, 0, 1, 0, -1.0, 0.0, 0.0, 0.0),
    (0, -1, 0, 1, 0, -1.0, 0.0, 0.0, 0.0),
    # Additional very small terms with time rates (kept for canonical completeness)
    (0, 0, 1, 0, 1, 0.0, 0.1, 0.0, 0.0),
    (0, 0, 1, 0, 1, 0.0, -0.1, 0.0, 0.0),
    (0, 2, 2, -2, 2, 0.0, 0.1, 0.0, 0.0),
    (0, -2, 2, -2, 2, 0.0, 0.1, 0.0, 0.0),
    (2, 0, 0, -2, 0, 0.0, 0.1, 0.0, 0.0),
    (2, 0, 2, -2, 1, 0.0, 0.1, 0.0, 0.0),
    (2, 0, 2, -2, 1, 0.0, -0.1, 0.0, 0.0),
]


def _nutation_iau1980(
    T: float, Lm_deg: float, Ls_deg: float, F_deg: float, D_deg: float, Om_deg: float
) -> tuple[float, float]:
    """Compute nutation in longitude and obliquity using the IAU 1980 106-term series.

    Args:
        T: Julian centuries since J2000.0.
        Lm_deg: Mean anomaly of the Moon (degrees).
        Ls_deg: Mean anomaly of the Sun (degrees).
        F_deg: Moon's argument of latitude (degrees).
        D_deg: Moon's elongation from the Sun (degrees).
        Om_deg: Longitude of the ascending node (degrees).

    Returns:
        Tuple (Δψ, Δε) in arcseconds.
    """
    # Convert arguments to radians
    Lm = _deg_to_rad(Lm_deg)
    Ls = _deg_to_rad(Ls_deg)
    F = _deg_to_rad(F_deg)
    D = _deg_to_rad(D_deg)
    Om = _deg_to_rad(Om_deg)

    dpsi_as = 0.0
    deps_as = 0.0
    for cD, cM, cMp, cF, cOm, ps, ps_t, pe, pe_t in _IAU1980_TERMS:
        arg = cD * D + cM * Ls + cMp * Lm + cF * F + cOm * Om
        s = math.sin(arg)
        c = math.cos(arg)
        # Linear time dependence per century (IAU 1980 convention)
        dpsi_as += ps * s + ps_t * s * T
        deps_as += pe * c + pe_t * c * T
    return dpsi_as, deps_as  # arcseconds


def _true_obliquity_rad(tt: float) -> float:
    """Return true obliquity of the ecliptic at given TT in radians.

    Args:
        tt: Julian Date (TT).

    Returns:
        True obliquity in radians (mean + nutation in obliquity).
    """
    T = _julian_centuries_TT_from_tt(tt)
    Lm, Ls, F, D, Om = _fundamental_arguments_deg(T)
    _dpsi_as, deps_as = _nutation_iau1980(T, Lm, Ls, F, D, Om)
    eps0_as = _mean_obliquity_arcsec(T)
    return _arcsec_to_rad(eps0_as + deps_as)


def _ecliptic_lon_lat_deg_of_date(apparent_vector: Apparent) -> tuple[float, float]:
    """Compute apparent ecliptic-of-date longitude and latitude in degrees.

    Args:
        apparent_vector: Apparent equatorial position of the Moon.

    Returns:
        Tuple of (longitude_deg, latitude_deg) in the true ecliptic of date.
    """
    eps = _true_obliquity_rad(apparent_vector.t.tt)

    # Apparent equatorial Cartesian (km)
    x, y, z = apparent_vector.position.km

    # Rotate around X by +eps (equatorial -> ecliptic of date)
    cos_e = math.cos(eps)
    sin_e = math.sin(eps)
    x_ecl = x
    y_ecl = y * cos_e + z * sin_e
    z_ecl = -y * sin_e + z * cos_e

    r_xy = math.hypot(x_ecl, y_ecl)
    lon = math.degrees(math.atan2(y_ecl, x_ecl))
    lat = math.degrees(math.atan2(z_ecl, r_xy))
    lon = (lon % 360.0 + 360.0) % 360.0
    return float(lon), float(lat)


# ----------------------------------------------------------------------


def _moon_parallax_angle_deg(distance_km: float) -> float:
    """Approximate equatorial horizontal parallax (degrees) from distance.

    Args:
        distance_km: Topocentric distance to the Moon in kilometers.

    Returns:
        Horizontal parallax angle in degrees.
    """
    Re_km = 6378.137
    x = Re_km / max(distance_km, 1e-6)
    x = min(1.0, max(0.0, x))
    return math.degrees(math.asin(x))


def _next_rise_set(
    eph: Ephemeris, lat: float, lon: float, alt_m: float, t_start: Time
) -> tuple[Time | None, Time | None]:
    """Compute next rise and set times for the Moon at observer location.

    Args:
        eph: Loaded ephemeris.
        lat: Latitude in degrees.
        lon: Longitude in degrees.
        alt_m: Elevation in meters.
        t_start: Start time for the search.

    Returns:
        A tuple (next_rise, next_set) as Skyfield Times or None if not found.
    """
    observer = wgs84.latlon(
        latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt_m
    )
    f = almanac.risings_and_settings(eph, eph["moon"], observer)
    t0 = t_start
    t1 = t_start + 7.0
    times, events = almanac.find_discrete(t0, t1, f)
    next_rise: Time | None = None
    next_set: Time | None = None
    for ti, ei in zip(times, events, strict=False):
        if ei and next_rise is None and ti.tt > t_start.tt:
            next_rise = ti
        if (not ei) and next_set is None and ti.tt > t_start.tt:
            next_set = ti
        if next_rise is not None and next_set is not None:
            break
    return next_rise, next_set


@dataclass
class _BrentResult:
    """Container for Brent extremum search results."""

    tt: float
    fval: float
    iterations: int


def _brent_extremum(
    f: Callable[[float], float],
    a: float,
    b: float,
    is_min: bool = True,
    tol: float = 1e-6,
    max_iter: int = 100,
) -> _BrentResult:
    """Generic Brent method to find an extremum of a univariate function.

    Args:
        f: Function mapping a scalar to a scalar.
        a: Left bound in the independent variable.
        b: Right bound in the independent variable.
        is_min: If True, search for minimum; otherwise maximum.
        tol: Absolute tolerance on the abscissa.
        max_iter: Maximum number of iterations.

    Returns:
        A _BrentResult containing location (tt), function value and iteration count.
    """
    g: Callable[[float], float] = (lambda x: -f(x)) if not is_min else f

    invphi = (math.sqrt(5) - 1) / 2
    invphi2 = (3 - math.sqrt(5)) / 2

    x = w = v = a + invphi2 * (b - a)
    fx = fw = fv = g(x)
    d = e = 0.0

    for it in range(max_iter):
        m = 0.5 * (a + b)
        tol1 = tol * abs(x) + 1e-12
        tol2 = 2.0 * tol1

        if abs(x - m) <= tol2 - 0.5 * (b - a):
            return _BrentResult(tt=x, fval=(fx if is_min else -fx), iterations=it)

        p = q = r = 0.0
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if q > 0:
                p = -p
            q = abs(q)
            parabolic_ok = (
                (abs(p) < abs(0.5 * q * e)) and (p > q * (a - x)) and (p < q * (b - x))
            )
            if parabolic_ok:
                d = p / q
                u = x + d
                if (u - a) < tol2 or (b - u) < tol2:
                    d = tol1 if (x < m) else -tol1
            else:
                e = (b - a) if (x < m) else (a - b)
                d = invphi * e
        else:
            e = (b - a) if (x < m) else (a - b)
            d = invphi * e

        u = x + (d if abs(d) >= tol1 else (tol1 if d > 0 else -tol1))
        fu = g(u)

        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v in (x, w):
                v, fv = u, fu

    return _BrentResult(tt=x, fval=(fx if is_min else -fx), iterations=max_iter)


def _refine_bracket(
    t_list: list[Time], y_list: list[float], kind: str = "max"
) -> tuple[float, float] | None:
    """Select a [a,b] bracket in tt for a local extremum based on slope changes.

    Args:
        t_list: Monotonic list of Time samples.
        y_list: Corresponding function values.
        kind: Either "max" or "min".

    Returns:
        A tuple (tt_a, tt_b) in TT days for refinement, or None if not found.
    """
    y = np.array(y_list)
    slopes = np.diff(y)
    if kind == "max":
        idxs = np.where((slopes[:-1] > 0) & (slopes[1:] < 0))[0] + 1
    else:
        idxs = np.where((slopes[:-1] < 0) & (slopes[1:] > 0))[0] + 1
    if len(idxs) == 0:
        return None
    idx = int(idxs[0])
    i0 = max(0, idx - 1)
    i2 = min(len(t_list) - 1, idx + 1)
    return (t_list[i0].tt, t_list[i2].tt)


def _find_next_apogee(eph: Ephemeris, ts: Timescale, t_start: Time) -> Time | None:
    """Find the next apogee after t_start using distance maximization.

    Args:
        eph: Loaded ephemeris.
        ts: Skyfield Timescale.
        t_start: Start time for the search.

    Returns:
        Skyfield Time at apogee, or None if not found in the window.
    """
    earth, moon = eph["earth"], eph["moon"]

    def geocentric_distance_km(t: Time) -> float:
        """Return geocentric distance Earth-Moon at time t in kilometers."""
        return earth.at(t).observe(moon).distance().km

    days_window = 40.0
    step_hours = 2.0
    steps = int(days_window * 24 / step_hours) + 1
    t_list: list[Time] = []
    d_list: list[float] = []
    for i in range(steps):
        dt_days = i * (step_hours / 24.0)
        t_i = t_start + dt_days
        t_list.append(t_i)
        d_list.append(geocentric_distance_km(t_i))

    bracket = _refine_bracket(t_list, d_list, kind="max")
    if bracket is None:
        return None
    tt_a, tt_b = bracket

    def f_tt(tt: float) -> float:
        """Distance in km as a function of TT (days)."""
        return geocentric_distance_km(ts.tt_jd(tt))

    res = _brent_extremum(f_tt, tt_a, tt_b, is_min=False, tol=1e-7, max_iter=200)
    t_ext = ts.tt_jd(res.tt)
    return t_ext if t_ext.tt > t_start.tt else None


def _find_next_perigee(eph: Ephemeris, ts: Timescale, t_start: Time) -> Time | None:
    """Find the next perigee after t_start using distance minimization.

    Args:
        eph: Loaded ephemeris.
        ts: Skyfield Timescale.
        t_start: Start time for the search.

    Returns:
        Skyfield Time at perigee, or None if not found in the window.
    """
    earth, moon = eph["earth"], eph["moon"]

    def geocentric_distance_km(t: Time) -> float:
        """Return geocentric distance Earth-Moon at time t in kilometers."""
        return earth.at(t).observe(moon).distance().km

    days_window = 40.0
    step_hours = 2.0
    steps = int(days_window * 24 / step_hours) + 1
    t_list: list[Time] = []
    d_list: list[float] = []
    for i in range(steps):
        dt_days = i * (step_hours / 24.0)
        t_i = t_start + dt_days
        t_list.append(t_i)
        d_list.append(geocentric_distance_km(t_i))

    bracket = _refine_bracket(t_list, d_list, kind="min")
    if bracket is None:
        return None
    tt_a, tt_b = bracket

    def f_tt(tt: float) -> float:
        """Distance in km as a function of TT (days)."""
        return geocentric_distance_km(ts.tt_jd(tt))

    res = _brent_extremum(f_tt, tt_a, tt_b, is_min=True, tol=1e-7, max_iter=200)
    t_ext = ts.tt_jd(res.tt)
    return t_ext if t_ext.tt > t_start.tt else None


def _find_phase_next(
    ts: Timescale, f: Any, t_start: Time, phase_value: int
) -> Time | None:
    """Find the next occurrence of a given moon phase after t_start.

    Args:
        ts: Skyfield Timescale.
        f: Almanac function for moon phases.
        t_start: Start time for the search.
        phase_value: Target phase identifier.

    Returns:
        The Skyfield Time of the next matching phase, or None if not found.
    """
    t0 = t_start
    t1 = t_start + 40.0
    times, phases = almanac.find_discrete(t0, t1, f)
    for ti, pv in zip(times, phases, strict=False):
        if int(pv) == int(phase_value) and ti.tt > t_start.tt:
            return ti
    return None


def _zodiac_sign_from_longitude_deg(lon_deg: float) -> str | None:
    """Map ecliptic longitude to a zodiac sign name.

    Args:
        lon_deg: Ecliptic longitude in degrees.

    Returns:
        Lowercase zodiac sign name or None for NaN inputs.
    """
    idx = int((lon_deg % 360.0) // 30) % 12
    names = [
        "aries",
        "taurus",
        "gemini",
        "cancer",
        "leo",
        "virgo",
        "libra",
        "scorpio",
        "sagittarius",
        "capricorn",
        "aquarius",
        "pisces",
    ]
    return names[idx]


def _degree_within_sign(lon_deg: float) -> float:
    """Return the degree within the current zodiac sign.

    Args:
        lon_deg: Ecliptic longitude in degrees.

    Returns:
        Degree within the sign in [0, 30).
    """
    return float((lon_deg % 360.0) % 30.0)


def _zodiac_icon(sign: str | None) -> str | None:
    """Return an MDI icon name for a zodiac sign.

    Args:
        sign: Lowercase zodiac sign name.

    Returns:
        Corresponding icon string or None if unknown.
    """
    if not sign:
        return None
    return {
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
    }.get(sign)


class MoonAstroCoordinator(DataUpdateCoordinator[dict[str, Any]]):
    """Coordinator computing lunar ephemerides and derived quantities for HA sensors."""

    def __init__(
        self,
        hass: HomeAssistant,
        lat: float,
        lon: float,
        elev: float,
        interval: timedelta,
    ) -> None:
        """Initialize the coordinator with observer and scheduling settings.

        Args:
            hass: Home Assistant instance.
            lat: Observer latitude in degrees.
            lon: Observer longitude in degrees.
            elev: Observer elevation in meters.
            interval: Update interval for the coordinator.
        """
        super().__init__(
            hass,
            logger=logging.getLogger(__name__),
            name="Moon Astro",
            update_interval=interval,
        )
        self._lat: float = float(lat)
        self._lon: float = float(lon)
        self._elev: float = float(elev)
        self._eph: Ephemeris | None = None
        self._ts: Timescale | None = None
        self._tz: ZoneInfo | None = None
        self._use_ha_tz: bool = False
        self._hass: HomeAssistant = hass

    @classmethod
    def from_config_entry(
        cls, hass: HomeAssistant, entry: ConfigEntry, interval: timedelta
    ) -> MoonAstroCoordinator:
        """Build the coordinator from a ConfigEntry.

        Args:
            hass: Home Assistant instance.
            entry: Config entry containing coordinates and options.
            interval: Update interval for the coordinator.

        Returns:
            A fully configured MoonAstroCoordinator instance.
        """
        data = entry.data
        c = cls(
            hass=hass,
            lat=float(data.get(CONF_LAT, hass.config.latitude)),
            lon=float(data.get(CONF_LON, hass.config.longitude)),
            elev=float(data.get(CONF_ALT, hass.config.elevation or 0)),
            interval=interval,
        )
        c._use_ha_tz = entry.options.get(CONF_USE_HA_TZ, True)
        c._tz = _tz_for_hass(hass) if c._use_ha_tz else _detect_timezone(c._lat, c._lon)
        return c

    async def _async_load_ephemeris(self) -> tuple[Ephemeris, Timescale]:
        """Load ephemerides and timescale asynchronously with caching.

        Returns:
            A tuple (ephemeris, timescale).
        """

        def _load() -> tuple[Ephemeris, Timescale]:
            """Blocking loader executed in the executor."""
            cache_dir = self._hass.config.path(CACHE_DIR_NAME)
            Path(cache_dir).mkdir(parents=True, exist_ok=True)
            load = Loader(cache_dir)
            eph: Ephemeris = load(DE440_FILE)
            ts: Timescale = load.timescale()
            return eph, ts

        return await self._hass.async_add_executor_job(_load)

    async def _async_update_data(self) -> dict[str, Any]:
        """Compute current lunar data and upcoming events for sensors.

        Returns:
            A dictionary with all computed keys ready to be exposed by entities.

        Raises:
            UpdateFailed: If an unexpected error occurs during calculations.
        """
        try:
            if self._eph is None or self._ts is None:
                self._eph, self._ts = await self._async_load_ephemeris()

            # Use local variables to satisfy type checkers after conditional load
            eph = self._eph
            ts = self._ts

            # At this point, eph and ts are guaranteed non-None by the conditional above
            assert eph is not None, "Ephemeris must be loaded"
            assert ts is not None, "Timescale must be loaded"

            now_utc = datetime.now(UTC)
            t = ts.from_datetime(now_utc)
            t_future = ts.from_datetime(now_utc + timedelta(hours=6))

            def _calc() -> dict[str, Any]:
                """Heavy computation executed in the executor thread."""
                topo_vec, az_deg, el_deg, dist_km = _topocentric_vectors(
                    eph, t, self._lat, self._lon, self._elev
                )
                geo_vec = _geocentric_vector(eph, t)

                ecl_lon_topo, ecl_lat_topo = _ecliptic_lon_lat_deg_of_date(topo_vec)
                ecl_lon_geo, ecl_lat_geo = _ecliptic_lon_lat_deg_of_date(geo_vec)

                illum_now = _moon_illumination_percentage(eph, t)
                parallax_deg = _moon_parallax_angle_deg(dist_km)

                phase = _moon_phase_name(eph, t, ts)

                illum_future = _moon_illumination_percentage(eph, t_future)
                waxing = illum_future > illum_now + 1e-6

                try:
                    t_rise, t_set = _next_rise_set(
                        eph, self._lat, self._lon, self._elev, t
                    )
                except _RECOVERABLE_SKYFIELD_ERRORS:
                    t_rise, t_set = None, None

                try:
                    f = almanac.moon_phases(eph)
                    t_new = _find_phase_next(ts, f, t, DARK_MOON)
                    t_first = _find_phase_next(ts, f, t, FIRST_QUARTER)
                    t_full = _find_phase_next(ts, f, t, FULL_MOON)
                    t_last = _find_phase_next(ts, f, t, LAST_QUARTER)
                except _RECOVERABLE_SKYFIELD_ERRORS:
                    t_new = t_first = t_full = t_last = None

                try:
                    t_apogee = _find_next_apogee(eph, ts, t)
                except _RECOVERABLE_NUMERIC_ERRORS:
                    t_apogee = None
                try:
                    t_perigee = _find_next_perigee(eph, ts, t)
                except _RECOVERABLE_NUMERIC_ERRORS:
                    t_perigee = None

                # Compute ecliptic lon/lat at the lunation instants (geocentric, true-of-date)
                lon_new: float | None = None
                lat_new: float | None = None
                lon_full: float | None = None
                lat_full: float | None = None
                try:
                    if t_new is not None:
                        v_new = _geocentric_vector(eph, t_new)
                        lon_new, lat_new = _ecliptic_lon_lat_deg_of_date(v_new)
                except _RECOVERABLE_NUMERIC_ERRORS as exc:
                    _LOGGER.debug("Failed ecliptic lon/lat at new moon: %r", exc)
                    lon_new = lat_new = None
                try:
                    if t_full is not None:
                        v_full = _geocentric_vector(eph, t_full)
                        lon_full, lat_full = _ecliptic_lon_lat_deg_of_date(v_full)
                except _RECOVERABLE_NUMERIC_ERRORS as exc:
                    _LOGGER.debug("Failed ecliptic lon/lat at full moon: %r", exc)
                    lon_full = lat_full = None

                zodiac_new: str | None = None
                zodiac_full: str | None = None
                zodiac_degree_new: float | None = None
                zodiac_degree_full: float | None = None
                try:
                    if lon_new is not None:
                        zodiac_new = _zodiac_sign_from_longitude_deg(lon_new)
                        zodiac_degree_new = round(_degree_within_sign(lon_new), 4)
                        _LOGGER.debug(
                            "Next new moon t=%s lon_moon_geo=%.6f° lat_moon_geo=%.6f° sign=%s deg_in_sign=%s",
                            _safe_time_iso(t_new, ZoneInfo("UTC")),
                            lon_new,
                            lat_new,
                            zodiac_new,
                            zodiac_degree_new,
                        )
                except (ValueError, ArithmeticError) as exc:
                    _LOGGER.debug("Failed zodiac_new computation: %r", exc)
                    zodiac_new = None
                    zodiac_degree_new = None
                try:
                    if lon_full is not None:
                        zodiac_full = _zodiac_sign_from_longitude_deg(lon_full)
                        zodiac_degree_full = round(_degree_within_sign(lon_full), 4)
                        _LOGGER.debug(
                            "Next full moon t=%s lon_moon_geo=%.6f° lat_moon_geo=%.6f° sign=%s deg_in_sign=%s",
                            _safe_time_iso(t_full, ZoneInfo("UTC")),
                            lon_full,
                            lat_full,
                            zodiac_full,
                            zodiac_degree_full,
                        )
                except (ValueError, ArithmeticError) as exc:
                    _LOGGER.debug("Failed zodiac_full computation: %r", exc)
                    zodiac_full = None
                    zodiac_degree_full = None

                zodiac_icon_new = _zodiac_icon(zodiac_new)
                zodiac_icon_full = _zodiac_icon(zodiac_full)

                above = el_deg > 0.0

                return {
                    KEY_PHASE: phase,
                    KEY_AZ: round(az_deg, 4),
                    KEY_EL: round(el_deg, 4),
                    KEY_ILLUM: round(illum_now, 3),
                    KEY_DISTANCE: round(dist_km, 3),
                    KEY_PARALLAX: round(parallax_deg, 4),
                    KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC: round(ecl_lon_topo, 6),
                    KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC: round(ecl_lat_topo, 6),
                    KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC: round(ecl_lon_geo, 6),
                    KEY_ECLIPTIC_LATITUDE_GEOCENTRIC: round(ecl_lat_geo, 6),
                    # The following four keys depend on optional event times; may be None/NaN
                    KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON: None
                    if lon_full is None
                    else round(lon_full, 6),
                    KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON: None
                    if lat_full is None
                    else round(lat_full, 6),
                    KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON: None
                    if lon_new is None
                    else round(lon_new, 6),
                    KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON: None
                    if lat_new is None
                    else round(lat_new, 6),
                    KEY_NEXT_RISE: _safe_time_iso(t_rise, self._tz),
                    KEY_NEXT_SET: _safe_time_iso(t_set, self._tz),
                    KEY_NEXT_APOGEE: _safe_time_iso(t_apogee, self._tz),
                    KEY_NEXT_PERIGEE: _safe_time_iso(t_perigee, self._tz),
                    KEY_NEXT_NEW_MOON: _safe_time_iso(t_new, self._tz),
                    KEY_NEXT_FIRST_QUARTER: _safe_time_iso(t_first, self._tz),
                    KEY_NEXT_FULL_MOON: _safe_time_iso(t_full, self._tz),
                    KEY_NEXT_LAST_QUARTER: _safe_time_iso(t_last, self._tz),
                    KEY_ZODIAC_SIGN_NEXT_NEW_MOON: zodiac_new,
                    KEY_ZODIAC_SIGN_NEXT_FULL_MOON: zodiac_full,
                    KEY_ZODIAC_DEGREE_NEXT_NEW_MOON: zodiac_degree_new,
                    KEY_ZODIAC_DEGREE_NEXT_FULL_MOON: zodiac_degree_full,
                    KEY_ZODIAC_ICON_NEXT_NEW_MOON: zodiac_icon_new,
                    KEY_ZODIAC_ICON_NEXT_FULL_MOON: zodiac_icon_full,
                    KEY_ABOVE_HORIZON: above,
                    KEY_WAXING: waxing,
                }

            return await self._hass.async_add_executor_job(_calc)
        except _RECOVERABLE_UPDATE_ERRORS as err:
            raise UpdateFailed(str(err)) from err
