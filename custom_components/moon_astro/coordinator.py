"""Data coordinator for Moon Astro.

This module implements the data coordinator responsible for computing high-precision
Moon position and related ephemerides for Home Assistant.

The code is structured to:
- keep deterministic results and minute-level timestamp stability
- limit redundant heavy computations
- keep full numerical precision while reducing CPU spikes on low-power devices
"""

from __future__ import annotations

import asyncio
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta, tzinfo
from itertools import pairwise
import logging
import math
from typing import Any

from skyfield import almanac
from skyfield.api import wgs84
from skyfield.jpllib import SpiceKernel
from skyfield.positionlib import Apparent
from skyfield.timelib import Time, Timescale
from skyfield.toposlib import GeographicPosition

from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant, callback
from homeassistant.helpers.event import async_track_point_in_time
from homeassistant.helpers.update_coordinator import DataUpdateCoordinator, UpdateFailed
from homeassistant.util import dt as dt_util

from .const import (
    CONF_ALT,
    CONF_HIGH_PRECISION,
    CONF_LAT,
    CONF_LON,
    DARK_MOON,
    DEFAULT_HIGH_PRECISION,
    DOMAIN,
    FIRST_QUARTER,
    FULL_MOON,
    FULL_MOON_NAMES,
    HIGH_PRECISION_BRACKET_EXPAND,
    HIGH_PRECISION_BRACKETS_TO_REFINE,
    HIGH_PRECISION_STEP_HOURS,
    HORIZON_ALTITUDE_DEG,
    KEY_ABOVE_HORIZON,
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
    KEY_WAXING,
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
    LAST_QUARTER,
    PHASE_CODES,
    PHASE_SEARCH_DAYS_AHEAD,
    PHASE_SEARCH_DAYS_BACK,
    PRINCIPAL_PHASE_WINDOW_DEG,
    RISE_SET_SEARCH_DAYS,
    STANDARD_PRECISION_BRACKET_EXPAND,
    STANDARD_PRECISION_STEP_HOURS,
    ZODIAC_SIGNS,
)

# Skyfield ships no type information; this alias names the loaded kernel.
type Ephemeris = SpiceKernel

# Errors raised by Skyfield searches and by numeric refinement; the affected value
# is reported as unavailable while the rest of the payload is kept.
_RECOVERABLE_SKYFIELD_ERRORS: tuple[type[Exception], ...] = (ValueError, RuntimeError)
_RECOVERABLE_NUMERIC_ERRORS: tuple[type[Exception], ...] = (
    ArithmeticError,
    ValueError,
    RuntimeError,
)
# Errors turning a whole coordinator update into a failure.
_RECOVERABLE_UPDATE_ERRORS: tuple[type[Exception], ...] = (
    OSError,
    ValueError,
    ArithmeticError,
    RuntimeError,
    KeyError,
    TypeError,
)

_LOGGER = logging.getLogger(__name__)

_EARTH_EQUATORIAL_RADIUS_KM = 6378.137

# Short phase names by Skyfield almanac.moon_phases value.
_PHASE_NAMES: dict[int, str] = {
    DARK_MOON: "new",
    FIRST_QUARTER: "first",
    FULL_MOON: "full",
    LAST_QUARTER: "last",
}

# Payload key of each phase event instant, by "<next|prev>_<phase>" source name.
_PHASE_TIMESTAMP_KEYS: dict[str, str] = {
    "next_new": KEY_NEXT_NEW_MOON,
    "next_first": KEY_NEXT_FIRST_QUARTER,
    "next_full": KEY_NEXT_FULL_MOON,
    "next_last": KEY_NEXT_LAST_QUARTER,
    "prev_new": KEY_PREVIOUS_NEW_MOON,
    "prev_first": KEY_PREVIOUS_FIRST_QUARTER,
    "prev_full": KEY_PREVIOUS_FULL_MOON,
    "prev_last": KEY_PREVIOUS_LAST_QUARTER,
}

_RISE_SET_KEYS: tuple[str, ...] = (
    KEY_NEXT_RISE,
    KEY_NEXT_SET,
    KEY_PREVIOUS_RISE,
    KEY_PREVIOUS_SET,
)

# Payload keys of the upcoming events used to schedule the next event-based refresh.
_NEXT_EVENT_KEYS: tuple[str, ...] = (
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_LAST_QUARTER,
    KEY_NEXT_APOGEE,
    KEY_NEXT_PERIGEE,
)

# Payload keys (longitude, latitude) of the ecliptic coordinates at lunations.
_LUNATION_ECLIPTIC_KEYS: dict[str, tuple[str, str]] = {
    "next_new": (
        KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON,
        KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON,
    ),
    "next_full": (
        KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON,
        KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON,
    ),
    "prev_new": (
        KEY_ECLIPTIC_LONGITUDE_PREVIOUS_NEW_MOON,
        KEY_ECLIPTIC_LATITUDE_PREVIOUS_NEW_MOON,
    ),
    "prev_full": (
        KEY_ECLIPTIC_LONGITUDE_PREVIOUS_FULL_MOON,
        KEY_ECLIPTIC_LATITUDE_PREVIOUS_FULL_MOON,
    ),
}

# Zodiac payload keys (sign, degree within sign) per longitude source.
_ZODIAC_KEYS: dict[str, tuple[str, str]] = {
    "current": (KEY_ZODIAC_SIGN_CURRENT_MOON, KEY_ZODIAC_DEGREE_CURRENT_MOON),
    "next_new": (KEY_ZODIAC_SIGN_NEXT_NEW_MOON, KEY_ZODIAC_DEGREE_NEXT_NEW_MOON),
    "next_full": (KEY_ZODIAC_SIGN_NEXT_FULL_MOON, KEY_ZODIAC_DEGREE_NEXT_FULL_MOON),
    "prev_new": (
        KEY_ZODIAC_SIGN_PREVIOUS_NEW_MOON,
        KEY_ZODIAC_DEGREE_PREVIOUS_NEW_MOON,
    ),
    "prev_full": (
        KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON,
        KEY_ZODIAC_DEGREE_PREVIOUS_FULL_MOON,
    ),
}

# -----------------------------------------------------------------------------
# Time conversion helpers
# -----------------------------------------------------------------------------


def _round_to_minute_utc(dt: datetime) -> datetime:
    """Round an aware datetime to the nearest minute and express it in UTC.

    Minute alignment keeps timestamp values stable from one refresh to the next.

    Args:
        dt: Aware datetime.

    Returns:
        Aware UTC datetime on a minute boundary.
    """
    utc = dt.astimezone(UTC)
    base = utc.replace(second=0, microsecond=0)
    return base + timedelta(minutes=1) if utc.second >= 30 else base


def _time_to_utc(t: Time | None) -> datetime | None:
    """Return a Skyfield Time as an aware UTC datetime rounded to the minute.

    Args:
        t: Skyfield Time, or None when the event was not found.

    Returns:
        Aware UTC datetime, or None.
    """
    return None if t is None else _round_to_minute_utc(t.utc_datetime())


def _time_to_local_datetime(t: Time, tz: tzinfo) -> datetime:
    """Return a Skyfield Time as an aware datetime in the given time zone.

    Args:
        t: Skyfield Time.
        tz: Target time zone.

    Returns:
        Aware datetime in tz.
    """
    return t.utc_datetime().astimezone(tz)


# -----------------------------------------------------------------------------
# Core astronomical helpers
# -----------------------------------------------------------------------------


def _topocentric_apparent(
    eph: Ephemeris, t: Time, observer: GeographicPosition
) -> tuple[Apparent, float, float]:
    """Return the topocentric apparent Moon position with its azimuth and altitude.

    Args:
        eph: Loaded ephemeris.
        t: Skyfield Time.
        observer: Observer position on the WGS84 ellipsoid.

    Returns:
        A tuple (apparent position, azimuth in degrees, altitude in degrees).
    """
    apparent = (eph["earth"] + observer).at(t).observe(eph["moon"]).apparent()
    alt, az, _distance = apparent.altaz()
    return apparent, float(az.degrees), float(alt.degrees)


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


def _julian_centuries_tt_from_tt(tt: float) -> float:
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
    T = _julian_centuries_tt_from_tt(tt)
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


def _moon_parallax_angle_deg(distance_km: float) -> float:
    """Return the equatorial horizontal parallax for a geocentric distance.

    Args:
        distance_km: Geocentric distance to the Moon in kilometers.

    Returns:
        Horizontal parallax in degrees.
    """
    return math.degrees(
        math.asin(
            _EARTH_EQUATORIAL_RADIUS_KM / max(distance_km, _EARTH_EQUATORIAL_RADIUS_KM)
        )
    )


# -----------------------------------------------------------------------------
# Phase naming
# -----------------------------------------------------------------------------


def _moon_phase_code(phase_deg: float) -> str:
    """Return the phase code for a Moon phase angle.

    A principal phase (new moon, first quarter, full moon, last quarter) is
    reported while the angle lies within PRINCIPAL_PHASE_WINDOW_DEG of 0, 90, 180
    or 270 degrees; the crescent and gibbous codes cover the rest of each quadrant.

    Args:
        phase_deg: Moon phase angle in degrees, 0 at new moon and 180 at full moon.

    Returns:
        A code listed in PHASE_CODES.
    """
    quadrant, offset = divmod(phase_deg % 360.0, 90.0)
    index = 2 * int(quadrant)
    if offset <= PRINCIPAL_PHASE_WINDOW_DEG:
        return PHASE_CODES[index % 8]
    if offset >= 90.0 - PRINCIPAL_PHASE_WINDOW_DEG:
        return PHASE_CODES[(index + 2) % 8]
    return PHASE_CODES[index + 1]


# -----------------------------------------------------------------------------
# Rise and set search
# -----------------------------------------------------------------------------


def _rise_set_around(
    eph: Ephemeris, t: Time, observer: GeographicPosition
) -> dict[str, Time | None]:
    """Return the moonrise and moonset instants closest to t on both sides.

    Args:
        eph: Loaded ephemeris.
        t: Reference time.
        observer: Observer position on the WGS84 ellipsoid.

    Returns:
        Mapping from the rise/set payload keys to event times, None when the event
        does not occur within RISE_SET_SEARCH_DAYS of t.
    """
    is_up = almanac.risings_and_settings(
        eph, eph["moon"], observer, horizon_degrees=HORIZON_ALTITUDE_DEG
    )
    times, rising = almanac.find_discrete(
        t - RISE_SET_SEARCH_DAYS, t + RISE_SET_SEARCH_DAYS, is_up
    )
    events: dict[str, Time | None] = dict.fromkeys(_RISE_SET_KEYS)
    for ti, is_rising in zip(times, rising, strict=True):
        if ti.tt <= t.tt:
            # Later events overwrite earlier ones: the last one before t is kept.
            events[KEY_PREVIOUS_RISE if is_rising else KEY_PREVIOUS_SET] = ti
        else:
            key = KEY_NEXT_RISE if is_rising else KEY_NEXT_SET
            if events[key] is None:
                events[key] = ti
    return events


# -----------------------------------------------------------------------------
# Numeric root/extremum helpers
# -----------------------------------------------------------------------------


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


def _brent_root(
    f: Callable[[float], float],
    a: float,
    b: float,
    *,
    tol: float = 1e-10,
    max_iter: int = 200,
) -> float | None:
    """Find a root of f(x)=0 on [a,b] using a Brent-style bracketing method.

    The function requires a sign change over the bracket. The implementation is
    designed to be deterministic and robust for numerical derivatives.

    Args:
        f: Scalar function.
        a: Left bracket bound.
        b: Right bracket bound.
        tol: Absolute tolerance on the abscissa.
        max_iter: Maximum iterations.

    Returns:
        The root location in the independent variable, or None if no sign change.
    """
    fa = float(f(a))
    fb = float(f(b))

    if math.isnan(fa) or math.isnan(fb):
        return None
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        return None

    c = a
    fc = fa
    d = e = b - a

    for _ in range(max_iter):
        if fb == 0.0:
            return b

        # Ensure |fb| <= |fc|
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb

        tol1 = 2.0 * 1e-12 + 0.5 * tol
        m = 0.5 * (c - b)

        if abs(m) <= tol1:
            return b

        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                # Secant
                p = 2.0 * m * s
                q = 1.0 - s
            else:
                # Inverse quadratic interpolation
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)

            if p > 0.0:
                q = -q
            p = abs(p)

            cond1 = 2.0 * p >= 3.0 * m * q - abs(tol1 * q)
            cond2 = p >= abs(0.5 * e * q)

            if cond1 or cond2 or q == 0.0:
                e = d
                d = m
            else:
                e = d
                d = p / q
        else:
            e = d
            d = m

        a = b
        fa = fb
        if abs(d) > tol1:
            b = b + d
        else:
            b = b + (tol1 if m > 0 else -tol1)

        fb = float(f(b))

        # Maintain the bracket on (b,c)
        if fb * fc > 0.0:
            c = a
            fc = fa
            d = e = b - a

    return None


def _derivative_distance_tt(
    ts: Timescale,
    f_tt: Callable[[float], float],
    tt: float,
    *,
    h_minutes: float,
) -> float:
    """Return a centered finite-difference derivative of distance wrt TT time.

    Args:
        ts: Skyfield Timescale (not used directly, kept for signature symmetry).
        f_tt: Distance function f(tt)->km.
        tt: Evaluation point (TT Julian date).
        h_minutes: Half-step in minutes for centered differences.

    Returns:
        Approximate derivative in km/day (TT days in denominator).
    """
    h_days = h_minutes / 1440.0
    da = float(f_tt(tt - h_days))
    db = float(f_tt(tt + h_days))
    return (db - da) / (2.0 * h_days)


def _find_derivative_sign_brackets(
    ts: Timescale,
    f_tt: Callable[[float], float],
    t_start: Time,
    *,
    search_backward: bool,
    days_window: float,
    step_minutes: float,
    h_minutes: float,
    max_candidates: int = 6,
) -> list[tuple[float, float]]:
    """Find brackets where the distance derivative changes sign.

    The function samples a centered finite-difference derivative on a regular TT grid
    and returns candidate intervals [tt_i, tt_{i+1}] where a sign change occurs.

    The returned brackets are ordered by proximity to t_start to reduce variability
    when multiple sign changes exist inside the search window.

    Args:
        ts: Skyfield Timescale.
        f_tt: Distance function f(tt)->km.
        t_start: Reference time.
        search_backward: Search backward if True, otherwise forward.
        days_window: Search window size in days.
        step_minutes: Sampling step in minutes.
        h_minutes: Half-step for centered derivative in minutes.
        max_candidates: Hard cap on returned brackets.

    Returns:
        A list of (tt_left, tt_right) candidate brackets, ordered by proximity to t_start.
    """
    step_days = step_minutes / 1440.0
    if step_days <= 0.0:
        return []

    tt_ref = float(t_start.tt)
    tt0 = tt_ref - float(days_window) if search_backward else tt_ref
    tt1 = tt_ref if search_backward else tt_ref + float(days_window)

    if not (tt1 > tt0):
        return []

    span = tt1 - tt0
    steps = int(math.ceil(span / step_days)) + 1
    if steps < 3:
        return []

    tts: list[float] = []
    vals: list[float] = []

    for i in range(steps):
        tt_i = min(tt0 + i * step_days, tt1)

        tts.append(tt_i)
        vals.append(_derivative_distance_tt(ts, f_tt, tt_i, h_minutes=h_minutes))

        if tt_i == tt1:
            break

    brackets: list[tuple[float, float]] = []
    for i in range(len(tts) - 1):
        a_tt = tts[i]
        b_tt = tts[i + 1]
        fa = vals[i]
        fb = vals[i + 1]

        if math.isnan(fa) or math.isnan(fb):
            continue

        if fa == 0.0:
            left = max(tt0, a_tt - step_days)
            right = min(tt1, a_tt + step_days)
            if right > left:
                brackets.append((left, right))
        elif fa * fb < 0.0:
            brackets.append((a_tt, b_tt))

    if not brackets:
        return []

    # Filter brackets by direction relative to t_start to avoid selecting "wrong side" candidates.
    # A bracket is kept if its center lies strictly on the expected side of t_start.
    filtered: list[tuple[float, float]] = []
    for a_tt, b_tt in brackets:
        center = 0.5 * (a_tt + b_tt)
        if search_backward:
            if center < tt_ref:
                filtered.append((a_tt, b_tt))
        elif center > tt_ref:
            filtered.append((a_tt, b_tt))

    if not filtered:
        return []

    # Sort by proximity to t_start (deterministic tie-breakers included).
    # Primary: |center - tt_ref|
    # Secondary: smaller bracket width first (tighter bracket is preferred)
    # Tertiary: chronological order for full determinism
    def _sort_key(item: tuple[float, float]) -> tuple[float, float, float]:
        """Return a deterministic sort key for a bracket."""
        a_tt, b_tt = item
        center = 0.5 * (a_tt + b_tt)
        width = max(0.0, b_tt - a_tt)
        return (abs(center - tt_ref), width, a_tt)

    filtered.sort(key=_sort_key)

    # Limit output size after sorting to keep CPU bounded.
    return filtered[: max(1, int(max_candidates))]


def _classify_extremum_kind_at_tt(
    ts: Timescale,
    f_tt: Callable[[float], float],
    tt_center: float,
) -> bool | None:
    """Classify an extremum as minimum or maximum around a candidate center.

    The classification uses a symmetric comparison around the center:
    - if f(center) <= f(center±delta): minimum
    - if f(center) >= f(center±delta): maximum

    Args:
        ts: Skyfield Timescale (not used directly, kept for signature symmetry).
        f_tt: Distance function f(tt)->km.
        tt_center: Candidate TT Julian date.

    Returns:
        True if the extremum is a minimum, False if maximum, or None if undecidable.
    """
    delta_days = 2.0 / 1440.0  # 2 minutes
    y0 = float(f_tt(tt_center))
    yl = float(f_tt(tt_center - delta_days))
    yr = float(f_tt(tt_center + delta_days))

    if math.isnan(y0) or math.isnan(yl) or math.isnan(yr):
        return None

    if y0 <= yl and y0 <= yr:
        return True
    if y0 >= yl and y0 >= yr:
        return False
    return None


def _select_extremum_candidate(
    ts: Timescale,
    f_tt: Callable[[float], float],
    candidates_tt: list[float],
    t_start: Time,
    *,
    is_min: bool,
    search_backward: bool,
) -> float | None:
    """Select the nearest candidate extremum time matching direction and kind.

    This selector is deterministic:
    - it filters candidates by direction relative to t_start
    - it filters by minimum/maximum classification
    - it selects the chronologically nearest valid candidate

    Args:
        ts: Skyfield Timescale.
        f_tt: Distance function f(tt)->km.
        candidates_tt: Candidate TT solutions.
        t_start: Reference time.
        is_min: True for perigee, False for apogee.
        search_backward: True for previous, False for next.

    Returns:
        Selected candidate TT or None.
    """
    ref_tt = float(t_start.tt)
    filtered: list[float] = []

    for tt in candidates_tt:
        if search_backward and not (tt < ref_tt):
            continue
        if (not search_backward) and not (tt > ref_tt):
            continue

        kind = _classify_extremum_kind_at_tt(ts, f_tt, tt)
        if kind is None:
            continue
        if bool(kind) != bool(is_min):
            continue

        filtered.append(tt)

    if not filtered:
        return None

    return max(filtered) if search_backward else min(filtered)


def _minute_validation_extremum(
    ts: Timescale,
    f_tt: Callable[[float], float],
    tt_center: float,
    *,
    is_min: bool,
) -> Time:
    """Snap the extremum to the most extreme minute among {t-1m, t, t+1m}.

    The input is a continuous-time solution (typically from Brent). The output is a
    Skyfield Time exactly on a minute boundary (TT), chosen by evaluating the function
    on the 3 neighboring minute instants.

    Args:
        ts: Skyfield Timescale.
        f_tt: Function evaluated on TT Julian dates.
        tt_center: Candidate solution (TT Julian date).
        is_min: True for perigee, False for apogee.

    Returns:
        A Skyfield Time located on the selected minute.
    """
    # Snap the continuous solution to the nearest UTC minute before probing neighbors.
    t_min = ts.from_datetime(_round_to_minute_utc(ts.tt_jd(tt_center).utc_datetime()))

    minute_days = 1.0 / 1440.0
    candidates = [
        t_min.tt - minute_days,
        t_min.tt,
        t_min.tt + minute_days,
    ]

    best_tt = candidates[1]
    best_val = float(f_tt(best_tt))

    for tt_i in candidates:
        val = float(f_tt(tt_i))
        if is_min:
            if val < best_val:
                best_val = val
                best_tt = tt_i
        elif val > best_val:
            best_val = val
            best_tt = tt_i

    return ts.tt_jd(best_tt)


# -----------------------------------------------------------------------------
# Extremum bracket detection for normal mode
# -----------------------------------------------------------------------------


def _refine_brackets(
    t_list: list[Time],
    y_list: list[float],
    kind: str = "max",
    *,
    expand: int = 1,
) -> list[tuple[float, float, int]]:
    """Return candidate TT brackets for local extrema in a sampled series.

    Args:
        t_list: Monotonic Time samples.
        y_list: Sample values.
        kind: "max" or "min".
        expand: Number of samples added on both sides of the bracket.

    Returns:
        List of (tt_a, tt_b, idx_center) ordered by idx_center.
    """
    if len(t_list) < 3 or len(y_list) < 3:
        return []

    eps = max(
        1e-12,
        1e-12 * max((abs(v) for v in y_list if math.isfinite(v)), default=0.0),
    )
    sgn = [
        1 if b - a > eps else -1 if b - a < -eps else 0
        for a, b in pairwise(y_list)
    ]

    want_left = 1 if kind == "max" else -1
    want_right = -1 if kind == "max" else 1

    n = len(y_list)
    candidates: list[int] = []

    for i in range(1, n - 1):
        k = i - 1
        while k >= 0 and sgn[k] == 0:
            k -= 1
        left_eff = sgn[k] if k >= 0 else 0

        j = i
        while j < len(sgn) and sgn[j] == 0:
            j += 1
        right_eff = sgn[j] if j < len(sgn) else 0

        if left_eff == want_left and right_eff == want_right:
            candidates.append(i)

    if not candidates:
        return []

    brackets: list[tuple[float, float, int]] = []
    for idx in candidates:
        left = idx
        while left > 0 and sgn[left - 1] == 0:
            left -= 1

        right = idx
        while right < len(sgn) and sgn[right] == 0:
            right += 1

        i0 = max(0, left - expand)
        i1 = min(n - 1, right + expand)
        if i1 <= i0:
            continue

        brackets.append((t_list[i0].tt, t_list[i1].tt, idx))

    brackets.sort(key=lambda item: item[2])
    return brackets


# -----------------------------------------------------------------------------
# Apsides (apogee/perigee) computation
# -----------------------------------------------------------------------------


def _make_geocentric_distance_tt_function(
    eph: Ephemeris,
    ts: Timescale,
) -> Callable[[float], float]:
    """Build a cached distance function f(tt)->km for geocentric Earth-Moon distance.

    The returned function is bounded-cached using quantized TT values to reduce
    repeated Skyfield evaluations during bracketing and iterative refinement.

    Args:
        eph: Loaded ephemeris.
        ts: Skyfield Timescale.

    Returns:
        A callable mapping TT Julian date to geocentric distance in kilometers.
    """
    earth, moon = eph["earth"], eph["moon"]

    dist_cache: dict[float, float] = {}
    dist_cache_order: list[float] = []
    dist_cache_max = 256

    def _quantize_tt(tt: float) -> float:
        """Quantize a TT Julian date for stable cache keys.

        Args:
            tt: TT Julian date.

        Returns:
            Quantized TT Julian date.
        """
        # 1e-10 day is ~8.64 microseconds; it avoids missing reuse due to tiny float jitter.
        return round(float(tt), 10)

    def _cache_put(key: float, value: float) -> None:
        """Insert a value into the bounded cache.

        Args:
            key: Quantized TT Julian date.
            value: Distance in km.

        Returns:
            None.
        """
        if key in dist_cache:
            return

        dist_cache[key] = value
        dist_cache_order.append(key)

        if len(dist_cache_order) > dist_cache_max:
            old = dist_cache_order.pop(0)
            dist_cache.pop(old, None)

    def f_tt(tt: float) -> float:
        """Return geocentric Earth-Moon distance in km as a function of TT Julian date.

        Args:
            tt: TT Julian date.

        Returns:
            Geocentric distance in kilometers.
        """
        key = _quantize_tt(tt)
        cached = dist_cache.get(key)
        if cached is not None:
            return cached

        value = float(earth.at(ts.tt_jd(tt)).observe(moon).distance().km)
        _cache_put(key, value)
        return value

    return f_tt


def _find_extremum_high_precision(
    ts: Timescale,
    f_tt: Callable[[float], float],
    t_start: Time,
    *,
    is_min: bool,
    search_backward: bool,
    days_window: float,
    tol: float,
    max_iter: int,
) -> Time | None:
    """Find an extremum using derivative sign-change brackets and root refinement.

    Args:
        ts: Skyfield Timescale.
        f_tt: Cached distance function f(tt)->km.
        t_start: Reference time.
        is_min: True to search for a minimum, False for a maximum.
        search_backward: True for previous, False for next.
        days_window: Search window in days.
        tol: Root solver tolerance.
        max_iter: Root solver max iterations.

    Returns:
        Skyfield Time for the extremum, or None if not found.
    """
    step_minutes = 10.0
    h_minutes = 20.0

    brackets = _find_derivative_sign_brackets(
        ts,
        f_tt,
        t_start,
        search_backward=search_backward,
        days_window=days_window,
        step_minutes=step_minutes,
        h_minutes=h_minutes,
        max_candidates=8,
    )

    _LOGGER.debug(
        "Distance extremum: mode=high_precision is_min=%s search_backward=%s brackets_found=%s",
        is_min,
        search_backward,
        len(brackets),
    )

    if not brackets:
        _LOGGER.debug("Distance extremum: no brackets found")
        return None

    def g_tt(tt: float) -> float:
        """Return derivative proxy d(distance)/d(TT day) at tt.

        Args:
            tt: TT Julian date.

        Returns:
            Derivative estimate.
        """
        return _derivative_distance_tt(ts, f_tt, tt, h_minutes=h_minutes)

    roots_tt: list[float] = []
    refined = 0

    for a_tt, b_tt in brackets:
        if refined >= HIGH_PRECISION_BRACKETS_TO_REFINE:
            break
        refined += 1

        root = _brent_root(g_tt, a_tt, b_tt, tol=tol, max_iter=max_iter)
        if root is None or not math.isfinite(root):
            continue

        roots_tt.append(float(root))

    _LOGGER.debug(
        "Distance extremum: refined_brackets=%s candidate_roots_found=%s",
        refined,
        len(roots_tt),
    )

    selected_tt = _select_extremum_candidate(
        ts,
        f_tt,
        roots_tt,
        t_start,
        is_min=is_min,
        search_backward=search_backward,
    )

    if selected_tt is None:
        _LOGGER.debug("Distance extremum: no candidate matched direction/kind filters")
        return None

    _LOGGER.debug(
        "Distance extremum: selected_candidate_tt=%s (search_backward=%s)",
        selected_tt,
        search_backward,
    )
    try:
        selected_distance_km = float(f_tt(selected_tt))
    except (ValueError, ArithmeticError, OverflowError) as exc:
        _LOGGER.debug(
            "Distance extremum: failed to evaluate distance at selected_tt (error=%r)",
            exc,
        )
    else:
        _LOGGER.debug(
            "Distance extremum: selected_distance_km=%.3f (expected_kind=%s)",
            selected_distance_km,
            "minimum" if is_min else "maximum",
        )

    t_ext = _minute_validation_extremum(ts, f_tt, selected_tt, is_min=is_min)

    _LOGGER.debug("Distance extremum: selected_time_utc=%s", t_ext.utc_iso())
    return t_ext


def _find_extremum_standard(
    ts: Timescale,
    f_tt: Callable[[float], float],
    eph: Ephemeris,
    t_start: Time,
    *,
    is_min: bool,
    search_backward: bool,
    days_window: float,
    step_hours: float,
    bracket_expand: int,
    tol: float,
    max_iter: int,
) -> Time | None:
    """Find an extremum by sampling distance and refining a local bracket.

    Args:
        ts: Skyfield Timescale.
        f_tt: Cached distance function f(tt)->km.
        eph: Loaded ephemeris (used for direct distance sampling on Time objects).
        t_start: Reference time.
        is_min: True to search for a minimum, False for a maximum.
        search_backward: True for previous, False for next.
        days_window: Search window in days.
        step_hours: Sampling step in hours.
        bracket_expand: Bracket expansion in samples.
        tol: Extremum solver tolerance.
        max_iter: Extremum solver max iterations.

    Returns:
        Skyfield Time for the extremum, or None if not found.
    """
    earth, moon = eph["earth"], eph["moon"]

    def geocentric_distance_km(t: Time) -> float:
        """Return geocentric Earth-Moon distance at time t in kilometers.

        Args:
            t: Skyfield Time.

        Returns:
            Geocentric distance in kilometers.
        """
        return earth.at(t).observe(moon).distance().km

    steps = int(days_window * 24.0 / step_hours) + 1
    t0 = (t_start - days_window) if search_backward else t_start

    t_list: list[Time] = []
    d_list: list[float] = []
    for i in range(steps):
        dt_days = i * (step_hours / 24.0)
        t_i = t0 + dt_days
        t_list.append(t_i)
        d_list.append(geocentric_distance_km(t_i))

    bracket_kind = "min" if is_min else "max"
    brackets = _refine_brackets(
        t_list, d_list, kind=bracket_kind, expand=bracket_expand
    )

    _LOGGER.debug(
        "Distance extremum: mode=standard is_min=%s search_backward=%s brackets_found=%s",
        is_min,
        search_backward,
        len(brackets),
    )
    if not brackets:
        return None

    t0_tt = t_start.tt
    chosen: tuple[float, float, int] | None = None

    if search_backward:
        for tt_a, tt_b, idx_center in brackets:
            if t_list[idx_center].tt < t0_tt:
                chosen = (tt_a, tt_b, idx_center)
            else:
                break
    else:
        for tt_a, tt_b, idx_center in brackets:
            if t_list[idx_center].tt > t0_tt:
                chosen = (tt_a, tt_b, idx_center)
                break

    if chosen is None:
        return None

    _LOGGER.debug(
        "Distance extremum: chosen_bracket_tt=(%s, %s) idx_center=%s",
        chosen[0],
        chosen[1],
        chosen[2],
    )

    tt_a, tt_b, _idx_center = chosen
    res = _brent_extremum(f_tt, tt_a, tt_b, is_min=is_min, tol=tol, max_iter=max_iter)
    t_ext = ts.tt_jd(res.tt)

    if search_backward:
        return t_ext if t_ext.tt < t_start.tt else None
    return t_ext if t_ext.tt > t_start.tt else None


def _find_geocentric_distance_extremum(
    eph: Ephemeris,
    ts: Timescale,
    t_start: Time,
    *,
    is_min: bool,
    search_backward: bool,
    days_window: float = 40.0,
    step_hours: float = STANDARD_PRECISION_STEP_HOURS,
    bracket_expand: int = STANDARD_PRECISION_BRACKET_EXPAND,
    tol: float = 1e-7,
    max_iter: int = 200,
) -> Time | None:
    """Find a local extremum of the geocentric Earth-Moon distance around a reference time.

    In high precision mode, a derivative-based strategy is used:
    - detect sign-change brackets for d(distance)/dt on a TT grid
    - refine root times using a bracketing root solver
    - classify each root as minimum/maximum
    - select the nearest candidate in the required direction
    - snap to a deterministic minute-aligned extremum

    In normal mode, the original distance-extremum strategy is preserved.

    Args:
        eph: Loaded ephemeris.
        ts: Skyfield Timescale.
        t_start: Reference time for the search.
        is_min: True to search for a minimum (perigee), False for a maximum (apogee).
        search_backward: True to search in the past window, False to search in the future window.
        days_window: Search window size in days.
        step_hours: Coarse sampling step in hours used to locate a bracketing interval (normal mode).
        bracket_expand: Extra samples included on each side of the detected bracket (normal mode).
        tol: Absolute tolerance on TT Julian date during refinement.
        max_iter: Maximum number of refinement iterations.

    Returns:
        A Skyfield Time for the extremum, or None if not found.
    """
    f_tt = _make_geocentric_distance_tt_function(eph, ts)

    use_high_precision_refine = (
        step_hours <= HIGH_PRECISION_STEP_HOURS
        and bracket_expand >= HIGH_PRECISION_BRACKET_EXPAND
    )

    if use_high_precision_refine:
        return _find_extremum_high_precision(
            ts,
            f_tt,
            t_start,
            is_min=is_min,
            search_backward=search_backward,
            days_window=days_window,
            tol=tol,
            max_iter=max_iter,
        )

    return _find_extremum_standard(
        ts,
        f_tt,
        eph,
        t_start,
        is_min=is_min,
        search_backward=search_backward,
        days_window=days_window,
        step_hours=step_hours,
        bracket_expand=bracket_expand,
        tol=tol,
        max_iter=max_iter,
    )


# -----------------------------------------------------------------------------
# Phase events and full moon naming
# -----------------------------------------------------------------------------


def _phase_events_around(
    t_ref: Time, times: Iterable[Time], phases: Iterable[int]
) -> dict[str, Time | None]:
    """Return the previous and next instant of each principal phase around t_ref.

    Args:
        t_ref: Reference time.
        times: Chronological event times found by almanac.find_discrete.
        phases: Skyfield phase values matching the times.

    Returns:
        Mapping from "<next|prev>_<phase>" source names to event times.
    """
    events: dict[str, Time | None] = dict.fromkeys(_PHASE_TIMESTAMP_KEYS)
    for ti, phase in zip(times, phases, strict=True):
        name = _PHASE_NAMES[int(phase)]
        if ti.tt <= t_ref.tt:
            events[f"prev_{name}"] = ti
        elif events[f"next_{name}"] is None:
            events[f"next_{name}"] = ti
    return events


def _full_moon_name_codes(
    full_moons: list[Time], t_ref: Time, tz: tzinfo
) -> tuple[str | None, str | None]:
    """Return the name codes of the next and previous full moons around t_ref.

    A full moon is named after its Gregorian month in the given time zone, or is a
    blue moon when it is the second one of that month.

    Args:
        full_moons: Chronological full moon instants surrounding t_ref.
        t_ref: Reference time.
        tz: Time zone defining the calendar month boundaries.

    Returns:
        A tuple (next full moon name, previous full moon name); an item is None when
        the corresponding full moon lies outside full_moons.
    """

    def name_at(index: int) -> str | None:
        """Return the name of the full moon at the given index, or None."""
        if not 0 <= index < len(full_moons):
            return None
        local = _time_to_local_datetime(full_moons[index], tz)
        if index > 0:
            before = _time_to_local_datetime(full_moons[index - 1], tz)
            if (before.year, before.month) == (local.year, local.month):
                return "blue_moon"
        return FULL_MOON_NAMES[local.month - 1]

    next_index = next(
        (i for i, ti in enumerate(full_moons) if ti.tt > t_ref.tt), len(full_moons)
    )
    return name_at(next_index), name_at(next_index - 1)


def _full_moon_alt_names_state_code(name_code: str | None) -> str | None:
    """Return the translation state code listing the alternative full moon names.

    A blue moon has no traditional alternative names and maps to an empty state.

    Args:
        name_code: Full moon name code.

    Returns:
        Translation state code, or None when no full moon is known.
    """
    if name_code is None:
        return None
    return "" if name_code == "blue_moon" else f"{name_code}_alt_names"


# -----------------------------------------------------------------------------
# Zodiac helpers
# -----------------------------------------------------------------------------


def _zodiac_sign_from_longitude_deg(lon_deg: float) -> str:
    """Map an ecliptic longitude to a zodiac sign code.

    Args:
        lon_deg: Ecliptic longitude in degrees.

    Returns:
        Lowercase zodiac sign code.
    """
    # The modulo guards against a normalized longitude rounding up to exactly 360.
    return ZODIAC_SIGNS[int((lon_deg % 360.0) // 30.0) % 12]


def _degree_within_sign(lon_deg: float) -> float:
    """Return the degree within the current zodiac sign.

    This function is intended to be fed with an unrounded ecliptic longitude to
    avoid boundary artifacts when the value is close to a sign cusp.

    Args:
        lon_deg: Ecliptic longitude in degrees.

    Returns:
        Degree within the sign in [0, 30).
    """
    lon_norm = (float(lon_deg) % 360.0 + 360.0) % 360.0
    return lon_norm - (math.floor(lon_norm / 30.0) * 30.0)


# -----------------------------------------------------------------------------
# Computation grouping helpers
# -----------------------------------------------------------------------------


class _Calc:
    """Private computation namespace.

    This class groups coordinator computation helpers to reduce the number of module-level
    symbols and to keep call sites compact.
    """

    @staticmethod
    def _round_or_none(value: float | None, ndigits: int) -> float | None:
        """Round a float value or return None.

        Args:
            value: Input float or None.
            ndigits: Decimal digits.

        Returns:
            Rounded float or None.
        """
        if value is None:
            return None
        return round(float(value), ndigits)

    @staticmethod
    def current(
        eph: Ephemeris, t: Time, observer: GeographicPosition
    ) -> tuple[dict[str, Any], float]:
        """Compute the current Moon position and the quantities derived from it.

        Args:
            eph: Loaded ephemeris.
            t: Current time.
            observer: Observer position on the WGS84 ellipsoid.

        Returns:
            The rounded payload fragment and the unrounded geocentric ecliptic
            longitude used by the zodiac computation.
        """
        topo_vec, az_deg, alt_deg = _topocentric_apparent(eph, t, observer)
        geo_vec = _geocentric_vector(eph, t)
        ecl_lon_topo, ecl_lat_topo = _ecliptic_lon_lat_deg_of_date(topo_vec)
        ecl_lon_geo, ecl_lat_geo = _ecliptic_lon_lat_deg_of_date(geo_vec)
        distance_km = float(geo_vec.distance().km)
        phase_deg = float(almanac.moon_phase(eph, t).degrees)
        payload = {
            KEY_PHASE: _moon_phase_code(phase_deg),
            KEY_AZIMUTH: round(az_deg, 4),
            KEY_ELEVATION: round(alt_deg, 4),
            KEY_ILLUMINATION: round(
                100.0 * float(geo_vec.fraction_illuminated(eph["sun"])), 3
            ),
            KEY_DISTANCE: round(distance_km, 3),
            KEY_PARALLAX: round(_moon_parallax_angle_deg(distance_km), 4),
            KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC: round(ecl_lon_topo, 6),
            KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC: round(ecl_lat_topo, 6),
            KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC: round(ecl_lon_geo, 6),
            KEY_ECLIPTIC_LATITUDE_GEOCENTRIC: round(ecl_lat_geo, 6),
            KEY_ABOVE_HORIZON: alt_deg > HORIZON_ALTITUDE_DEG,
            KEY_WAXING: phase_deg < 180.0,
        }
        return payload, ecl_lon_geo

    @staticmethod
    def rise_set(
        eph: Ephemeris, t: Time, observer: GeographicPosition
    ) -> dict[str, datetime | None]:
        """Compute the previous and next moonrise and moonset.

        Args:
            eph: Loaded ephemeris.
            t: Reference time.
            observer: Observer position on the WGS84 ellipsoid.

        Returns:
            A payload fragment with the rise and set instants.
        """
        try:
            events = _rise_set_around(eph, t, observer)
        except _RECOVERABLE_SKYFIELD_ERRORS:
            events = dict.fromkeys(_RISE_SET_KEYS)
        return {key: _time_to_utc(ti) for key, ti in events.items()}

    @staticmethod
    def apsis(
        eph: Ephemeris, ts: Timescale, t: Time, *, high_precision: bool
    ) -> dict[str, datetime | None]:
        """Compute the previous and next apogee and perigee.

        Args:
            eph: Loaded ephemeris.
            ts: Skyfield Timescale.
            t: Reference time.
            high_precision: Use the finer, more CPU intensive search.

        Returns:
            A payload fragment with the apsis instants.
        """
        step_hours = (
            HIGH_PRECISION_STEP_HOURS if high_precision else STANDARD_PRECISION_STEP_HOURS
        )
        bracket_expand = (
            HIGH_PRECISION_BRACKET_EXPAND
            if high_precision
            else STANDARD_PRECISION_BRACKET_EXPAND
        )
        payload: dict[str, datetime | None] = {}
        for key, is_min, search_backward in (
            (KEY_NEXT_APOGEE, False, False),
            (KEY_NEXT_PERIGEE, True, False),
            (KEY_PREVIOUS_APOGEE, False, True),
            (KEY_PREVIOUS_PERIGEE, True, True),
        ):
            try:
                t_event = _find_geocentric_distance_extremum(
                    eph,
                    ts,
                    t,
                    is_min=is_min,
                    search_backward=search_backward,
                    step_hours=step_hours,
                    bracket_expand=bracket_expand,
                )
            except _RECOVERABLE_NUMERIC_ERRORS:
                t_event = None
            payload[key] = _time_to_utc(t_event)
        _LOGGER.debug(
            "Apsis computation (high_precision=%s): %s", high_precision, payload
        )
        return payload

    @staticmethod
    def phases_and_names(
        eph: Ephemeris, t: Time, tz: tzinfo
    ) -> tuple[dict[str, Any], dict[str, Time | None]]:
        """Compute the phase event instants and the full moon name codes.

        A single discrete search covers the previous two lunations and the next one.

        Args:
            eph: Loaded ephemeris.
            t: Reference time.
            tz: Time zone defining the calendar months used for full moon names.

        Returns:
            The payload fragment and the raw event times by source name, reused by
            the lunation ecliptic computation.
        """
        events: dict[str, Time | None]
        try:
            times, phases = almanac.find_discrete(
                t - PHASE_SEARCH_DAYS_BACK,
                t + PHASE_SEARCH_DAYS_AHEAD,
                almanac.moon_phases(eph),
            )
        except _RECOVERABLE_SKYFIELD_ERRORS:
            events = dict.fromkeys(_PHASE_TIMESTAMP_KEYS)
            next_name = prev_name = None
        else:
            events = _phase_events_around(t, times, phases)
            next_name, prev_name = _full_moon_name_codes(
                [ti for ti, phase in zip(times, phases, strict=True) if phase == FULL_MOON],
                t,
                tz,
            )
        payload: dict[str, Any] = {
            key: _time_to_utc(events[source])
            for source, key in _PHASE_TIMESTAMP_KEYS.items()
        }
        payload |= {
            KEY_NEXT_FULL_MOON_NAME: next_name,
            KEY_NEXT_FULL_MOON_ALT_NAMES: _full_moon_alt_names_state_code(next_name),
            KEY_PREVIOUS_FULL_MOON_NAME: prev_name,
            KEY_PREVIOUS_FULL_MOON_ALT_NAMES: _full_moon_alt_names_state_code(
                prev_name
            ),
        }
        return payload, events

    @staticmethod
    def lunation_ecliptics(
        eph: Ephemeris, events: Mapping[str, Time | None]
    ) -> tuple[dict[str, Any], dict[str, float | None]]:
        """Compute the geocentric ecliptic coordinates at the surrounding lunations.

        Args:
            eph: Loaded ephemeris.
            events: Raw event times by source name.

        Returns:
            The rounded payload fragment and the unrounded longitudes by source name,
            reused by the zodiac computation.
        """
        payload: dict[str, Any] = {}
        raw_lons: dict[str, float | None] = {}
        for source, (lon_key, lat_key) in _LUNATION_ECLIPTIC_KEYS.items():
            lon: float | None = None
            lat: float | None = None
            if (t_event := events[source]) is not None:
                try:
                    lon, lat = _ecliptic_lon_lat_deg_of_date(
                        _geocentric_vector(eph, t_event)
                    )
                except _RECOVERABLE_NUMERIC_ERRORS as exc:
                    _LOGGER.debug(
                        "Ecliptic coordinates at %s unavailable: %r", source, exc
                    )
            payload[lon_key] = _Calc._round_or_none(lon, 6)
            payload[lat_key] = _Calc._round_or_none(lat, 6)
            raw_lons[source] = lon
        return payload, raw_lons

    @staticmethod
    def zodiac(longitudes: Mapping[str, float | None]) -> dict[str, Any]:
        """Compute zodiac sign and degree-in-sign keys for a set of longitudes.

        Args:
            longitudes: Unrounded geocentric ecliptic longitudes in degrees (or None
                when unavailable), keyed by a source name listed in _ZODIAC_KEYS.

        Returns:
            A payload dictionary fragment with the sign and degree keys of each source.
        """
        payload: dict[str, Any] = {}
        for source, lon_deg in longitudes.items():
            sign_key, degree_key = _ZODIAC_KEYS[source]
            payload[sign_key] = (
                None if lon_deg is None else _zodiac_sign_from_longitude_deg(lon_deg)
            )
            payload[degree_key] = (
                None if lon_deg is None else round(_degree_within_sign(lon_deg), 4)
            )
        return payload


# -----------------------------------------------------------------------------
# Coordinator
# -----------------------------------------------------------------------------


class MoonAstroCoordinator(DataUpdateCoordinator[dict[str, Any]]):
    """Coordinator computing lunar ephemerides and derived quantities for HA sensors."""

    def __init__(
        self,
        hass: HomeAssistant,
        entry: MoonAstroConfigEntry,
        *,
        eph: Ephemeris,
        ts: Timescale,
        interval: timedelta,
    ) -> None:
        """Initialize the coordinator with the observer and scheduling settings.

        Args:
            hass: Home Assistant instance.
            entry: Config entry providing the observer coordinates.
            eph: Loaded ephemeris shared with the other coordinators.
            ts: Loaded timescale shared with the other coordinators.
            interval: Update interval for the coordinator.
        """
        super().__init__(
            hass,
            logger=_LOGGER,
            config_entry=entry,
            name="Moon Astro",
            update_interval=interval,
        )
        self._eph = eph
        self._ts = ts
        self._observer = wgs84.latlon(
            latitude_degrees=float(entry.data.get(CONF_LAT, hass.config.latitude)),
            longitude_degrees=float(entry.data.get(CONF_LON, hass.config.longitude)),
            elevation_m=float(entry.data.get(CONF_ALT, hass.config.elevation)),
        )

    async def _async_compute_payload(self) -> dict[str, Any]:
        """Compute the payload from two concurrent executor jobs.

        The position and rise/set searches are independent and run in parallel; the
        zodiac derivation is plain arithmetic and runs inline.

        Returns:
            A dictionary with all computed keys ready to be exposed by entities.
        """
        t = self._ts.from_datetime(_round_to_minute_utc(datetime.now(UTC)))
        (current_payload, current_lon_geo), rise_set_payload = await asyncio.gather(
            self.hass.async_add_executor_job(
                _Calc.current, self._eph, t, self._observer
            ),
            self.hass.async_add_executor_job(
                _Calc.rise_set, self._eph, t, self._observer
            ),
        )
        return {
            **current_payload,
            **rise_set_payload,
            **_Calc.zodiac({"current": current_lon_geo}),
        }

    async def _async_update_data(self) -> dict[str, Any]:
        """Compute current lunar data for sensors.

        Returns:
            A dictionary with all computed keys ready to be exposed by entities.

        Raises:
            UpdateFailed: If a recoverable error occurs during calculations.
        """
        try:
            return await self._async_compute_payload()
        except _RECOVERABLE_UPDATE_ERRORS as err:
            raise UpdateFailed(f"Moon position computation failed: {err}") from err


# -----------------------------------------------------------------------------
# Events Coordinator
# -----------------------------------------------------------------------------


class MoonAstroEventsCoordinator(DataUpdateCoordinator[dict[str, Any]]):
    """Coordinator computing rare event-based lunar values.

    This coordinator focuses on values that only change when crossing an astronomical
    event boundary. A lightweight periodic fallback is kept to ensure eventual
    resynchronization.
    """

    def __init__(
        self,
        hass: HomeAssistant,
        entry: MoonAstroConfigEntry,
        *,
        eph: Ephemeris,
        ts: Timescale,
        interval: timedelta,
        tz: tzinfo,
    ) -> None:
        """Initialize the coordinator with computation and scheduling settings.

        Args:
            hass: Home Assistant instance.
            entry: Config entry providing the precision options.
            eph: Loaded ephemeris shared with the other coordinators.
            ts: Loaded timescale shared with the other coordinators.
            interval: Fallback update interval.
            tz: Time zone used to localize event timestamps.
        """
        super().__init__(
            hass,
            logger=_LOGGER,
            config_entry=entry,
            name="Moon Astro Events",
            update_interval=interval,
        )
        self._eph = eph
        self._ts = ts
        self._tz = tz
        self._high_precision = bool(
            entry.options.get(CONF_HIGH_PRECISION, DEFAULT_HIGH_PRECISION)
        )

        self._unsub_next_event: Callable[[], None] | None = None
        self._next_refresh_utc: datetime | None = None

        # Dedicated executor to avoid monopolizing Home Assistant's shared thread pool.
        # A single worker enforces determinism and prevents concurrent heavy computations.
        self._executor = ThreadPoolExecutor(
            max_workers=1,
            thread_name_prefix="moon_astro_events",
        )

        # Track an in-flight refresh task to avoid spawning multiple long-running
        # computations when many refresh requests happen close together.
        self._inflight_task: asyncio.Task[None] | None = None

    @property
    def next_refresh_utc(self) -> datetime | None:
        """Return the next scheduled refresh time in UTC.

        This timestamp represents when the coordinator expects to refresh event-based
        values next (either after the next event boundary or via the periodic fallback).

        Returns:
            A timezone-aware UTC datetime, or None if not scheduled.
        """
        if self._next_refresh_utc is not None:
            return self._next_refresh_utc

        interval = self.update_interval
        if interval is None:
            return None

        return datetime.now(UTC) + interval

    def _cancel_next_event_timer(self) -> None:
        """Cancel the scheduled next event refresh timer.

        Returns:
            None.
        """
        if self._unsub_next_event is not None:
            self._unsub_next_event()
            self._unsub_next_event = None
        self._next_refresh_utc = None

    def _schedule_next_event_refresh(self, data: dict[str, Any]) -> None:
        """Schedule a refresh shortly after the earliest upcoming event.

        Args:
            data: Coordinator data containing the next event timestamps.

        Returns:
            None.
        """
        self._cancel_next_event_timer()

        candidates = [
            dt for key in _NEXT_EVENT_KEYS if (dt := data.get(key)) is not None
        ]
        if not candidates:
            self._next_refresh_utc = (
                None
                if (interval := self.update_interval) is None
                else datetime.now(UTC) + interval
            )
            return

        # Refresh slightly after the boundary to avoid edge instability. An event
        # already in the past (clock jump, delayed startup) triggers a prompt refresh.
        next_event = min(candidates)
        now = dt_util.utcnow()
        when = next_event + timedelta(minutes=2)
        if when <= now:
            when = now + timedelta(seconds=30)
        self._next_refresh_utc = when

        @callback
        def _on_event_boundary(_: datetime) -> None:
            """Request a refresh once the event boundary has passed."""
            self._unsub_next_event = None
            self.hass.async_create_background_task(
                self.async_request_refresh(), name=f"{DOMAIN}-events-boundary-refresh"
            )

        _LOGGER.debug(
            "Events scheduler: refresh scheduled at %s for the event at %s",
            when.isoformat(),
            next_event.isoformat(),
        )
        self._unsub_next_event = async_track_point_in_time(
            self.hass, _on_event_boundary, when
        )

    async def _async_compute_events_payload(self) -> dict[str, Any]:
        """Compute the event-based payload in the dedicated executor.

        Returns:
            A dictionary containing only event-based keys.
        """
        eph, ts, tz = self._eph, self._ts, self._tz
        high_precision = self._high_precision
        t = ts.from_datetime(_round_to_minute_utc(datetime.now(UTC)))

        def _calc() -> dict[str, Any]:
            """Run the heavy event computations in the executor thread."""
            phase_payload, events = _Calc.phases_and_names(eph, t, tz)
            ecliptic_payload, raw_lons = _Calc.lunation_ecliptics(eph, events)
            return {
                **phase_payload,
                **_Calc.apsis(eph, ts, t, high_precision=high_precision),
                **ecliptic_payload,
                **_Calc.zodiac(raw_lons),
            }

        return await asyncio.get_running_loop().run_in_executor(self._executor, _calc)

    async def _async_run_refresh_job(self) -> None:
        """Compute the event-based payload and publish it.

        The previous data is kept when the computation fails; listeners are then
        notified of the unavailability through the coordinator error path.

        Returns:
            None.
        """
        try:
            data = await self._async_compute_events_payload()
        except _RECOVERABLE_UPDATE_ERRORS as err:
            self.async_set_update_error(err)
            return
        self.async_set_updated_data(data)
        self._schedule_next_event_refresh(data)

    async def _async_update_data(self) -> dict[str, Any]:
        """Start a background computation and return the latest available data.

        The heavy computation never runs inside the coordinator update path so that
        scheduled and requested refreshes return immediately; the background job
        publishes the new payload through async_set_updated_data.

        Returns:
            The latest event-based data, possibly stale while a job is running.
        """
        if self._inflight_task is None or self._inflight_task.done():
            self._inflight_task = self.hass.async_create_background_task(
                self._async_run_refresh_job(), name=f"{DOMAIN}-events-refresh"
            )
            if self._next_refresh_utc is None and (
                interval := self.update_interval
            ) is not None:
                self._next_refresh_utc = datetime.now(UTC) + interval
        else:
            _LOGGER.debug("Events coordinator: computation already running")
        return self.data or {}

    async def async_shutdown(self) -> None:
        """Stop scheduled refreshes and release the dedicated executor.

        Home Assistant calls this method through the unload callback registered on
        the config entry; it may run more than once and stays idempotent.

        Returns:
            None.
        """
        await super().async_shutdown()
        self._cancel_next_event_timer()

        if (task := self._inflight_task) is not None and not task.done():
            task.cancel()

        # Non-blocking: pending jobs are dropped and the worker exits after the
        # currently running job, if any.
        self._executor.shutdown(wait=False, cancel_futures=True)


# -----------------------------------------------------------------------------
# Config entry runtime data
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MoonAstroRuntimeData:
    """Runtime objects attached to a loaded config entry."""

    coordinator: MoonAstroCoordinator
    events_coordinator: MoonAstroEventsCoordinator


type MoonAstroConfigEntry = ConfigEntry[MoonAstroRuntimeData]