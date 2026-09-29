"""Skyfield computations for Moon Astro.

The functions of this module run in the executor. They take the loaded DE440
kernel and a reference time and return the payload fragments published by the
coordinators: angles in degrees, instants as aware UTC datetimes rounded to the
minute, state codes listed in const.
"""

from __future__ import annotations

from collections.abc import Iterable
from datetime import UTC, datetime, timedelta, tzinfo
import math
from typing import Any

from skyfield import almanac
from skyfield.framelib import ecliptic_frame
from skyfield.jpllib import SpiceKernel
from skyfield.positionlib import ICRF
from skyfield.searchlib import find_maxima, find_minima
from skyfield.timelib import Time
from skyfield.toposlib import GeographicPosition

from .const import (
    APSIS_EPSILON_DAYS,
    APSIS_SEARCH_DAYS,
    APSIS_STEP_DAYS,
    DARK_MOON,
    FIRST_QUARTER,
    FULL_MOON,
    FULL_MOON_NAMES,
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
    ZODIAC_SIGNS,
)

_EARTH_EQUATORIAL_RADIUS_KM = 6378.137

# Payload keys of the (previous, next) instant of each kind of event.
_RISE_KEYS = (KEY_PREVIOUS_RISE, KEY_NEXT_RISE)
_SET_KEYS = (KEY_PREVIOUS_SET, KEY_NEXT_SET)
_PERIGEE_KEYS = (KEY_PREVIOUS_PERIGEE, KEY_NEXT_PERIGEE)
_APOGEE_KEYS = (KEY_PREVIOUS_APOGEE, KEY_NEXT_APOGEE)
_PHASE_KEYS: dict[int, tuple[str, str]] = {
    DARK_MOON: (KEY_PREVIOUS_NEW_MOON, KEY_NEXT_NEW_MOON),
    FIRST_QUARTER: (KEY_PREVIOUS_FIRST_QUARTER, KEY_NEXT_FIRST_QUARTER),
    FULL_MOON: (KEY_PREVIOUS_FULL_MOON, KEY_NEXT_FULL_MOON),
    LAST_QUARTER: (KEY_PREVIOUS_LAST_QUARTER, KEY_NEXT_LAST_QUARTER),
}

# Full moon instant key with its (name, alternative names) keys.
_FULL_MOON_NAME_KEYS: tuple[tuple[str, str, str], ...] = (
    (
        KEY_PREVIOUS_FULL_MOON,
        KEY_PREVIOUS_FULL_MOON_NAME,
        KEY_PREVIOUS_FULL_MOON_ALT_NAMES,
    ),
    (KEY_NEXT_FULL_MOON, KEY_NEXT_FULL_MOON_NAME, KEY_NEXT_FULL_MOON_ALT_NAMES),
)

# Ecliptic (longitude, latitude) and zodiac (sign, degree) keys of each lunation key.
_LUNATION_KEYS: dict[str, tuple[str, str, str, str]] = {
    KEY_PREVIOUS_NEW_MOON: (
        KEY_ECLIPTIC_LONGITUDE_PREVIOUS_NEW_MOON,
        KEY_ECLIPTIC_LATITUDE_PREVIOUS_NEW_MOON,
        KEY_ZODIAC_SIGN_PREVIOUS_NEW_MOON,
        KEY_ZODIAC_DEGREE_PREVIOUS_NEW_MOON,
    ),
    KEY_PREVIOUS_FULL_MOON: (
        KEY_ECLIPTIC_LONGITUDE_PREVIOUS_FULL_MOON,
        KEY_ECLIPTIC_LATITUDE_PREVIOUS_FULL_MOON,
        KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON,
        KEY_ZODIAC_DEGREE_PREVIOUS_FULL_MOON,
    ),
    KEY_NEXT_NEW_MOON: (
        KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON,
        KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON,
        KEY_ZODIAC_SIGN_NEXT_NEW_MOON,
        KEY_ZODIAC_DEGREE_NEXT_NEW_MOON,
    ),
    KEY_NEXT_FULL_MOON: (
        KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON,
        KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON,
        KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
        KEY_ZODIAC_DEGREE_NEXT_FULL_MOON,
    ),
}

# -----------------------------------------------------------------------------
# Time helpers
# -----------------------------------------------------------------------------


def round_to_minute_utc(dt: datetime) -> datetime:
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
    return None if t is None else round_to_minute_utc(t.utc_datetime())


def _surrounding(
    t_ref: Time, times: Iterable[Time]
) -> tuple[Time | None, Time | None]:
    """Return the last event at or before t_ref and the first event after it.

    Args:
        t_ref: Reference time.
        times: Chronological event times.

    Returns:
        A tuple (previous, next); an item is None when no such event exists.
    """
    previous = following = None
    for ti in times:
        if ti.tt > t_ref.tt:
            following = ti
            break
        previous = ti
    return previous, following


# -----------------------------------------------------------------------------
# Positions and coordinates
# -----------------------------------------------------------------------------


def _topocentric_apparent(
    eph: SpiceKernel, t: Time, observer: GeographicPosition
) -> tuple[ICRF, float, float]:
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


def _geocentric_apparent(eph: SpiceKernel, t: Time) -> ICRF:
    """Return the geocentric apparent Moon position.

    Args:
        eph: Loaded ephemeris.
        t: Skyfield Time.

    Returns:
        Apparent position seen from the geocenter.
    """
    return eph["earth"].at(t).observe(eph["moon"]).apparent()


def _ecliptic_lon_lat_deg(position: ICRF) -> tuple[float, float]:
    """Return the coordinates of a position in the true ecliptic of date.

    Skyfield rotates the ICRS vector with its IAU 2006/2000A precession-nutation
    matrix, then by the true obliquity of date.

    Args:
        position: Position carrying its ICRS vector and its time.

    Returns:
        A tuple (longitude in [0, 360), latitude) in degrees.
    """
    lat, lon, _distance = position.frame_latlon(ecliptic_frame)
    return float(lon.degrees) % 360.0, float(lat.degrees)


class _GeocentricDistance:
    """Geometric Earth-Moon distance in kilometers as a function of time.

    Instances are the callables expected by the Skyfield extremum searches, which
    read the step_days attribute to size their coarse sampling grid.
    """

    step_days = APSIS_STEP_DAYS

    def __init__(self, eph: SpiceKernel) -> None:
        """Bind the Earth and Moon segments of the ephemeris."""
        self._earth_to_moon = eph["moon"] - eph["earth"]

    def __call__(self, t: Time) -> Any:
        """Return the distance at t; an array when t holds several times."""
        return self._earth_to_moon.at(t).distance().km


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


def _round_or_none(value: float | None, ndigits: int) -> float | None:
    """Round a value, passing None through.

    Args:
        value: Value to round, or None.
        ndigits: Number of decimals.

    Returns:
        The rounded value, or None.
    """
    return None if value is None else round(value, ndigits)


# -----------------------------------------------------------------------------
# Phase, zodiac and full moon naming
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


def _zodiac(lon_deg: float | None, sign_key: str, degree_key: str) -> dict[str, Any]:
    """Return the zodiac sign and the degree within the sign of an ecliptic longitude.

    Args:
        lon_deg: Unrounded ecliptic longitude in [0, 360), or None when unavailable.
        sign_key: Payload key of the sign code.
        degree_key: Payload key of the degree within the sign.

    Returns:
        A payload fragment with both keys, None-valued when the longitude is None.
    """
    if lon_deg is None:
        return {sign_key: None, degree_key: None}
    return {
        sign_key: ZODIAC_SIGNS[int(lon_deg // 30.0) % 12],
        degree_key: round(lon_deg % 30.0, 4),
    }


def _full_moon_name_codes(full_moons: Iterable[Time], tz: tzinfo) -> dict[datetime, str]:
    """Return the name code of each full moon of a chronological sequence.

    A full moon is named after its Gregorian month in the given time zone, or is a
    blue moon when the preceding full moon fell in the same month; the first full
    moon of the sequence can therefore never be recognized as a blue moon.

    Args:
        full_moons: Chronological full moon instants.
        tz: Time zone defining the calendar month boundaries.

    Returns:
        Name codes keyed by the minute-rounded UTC instant published in the payload.
    """
    names: dict[datetime, str] = {}
    previous_month: tuple[int, int] | None = None
    for ti in full_moons:
        utc = ti.utc_datetime()
        local = utc.astimezone(tz)
        month = (local.year, local.month)
        names[round_to_minute_utc(utc)] = (
            "blue_moon" if month == previous_month else FULL_MOON_NAMES[local.month - 1]
        )
        previous_month = month
    return names


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


def _snap_to_minute(
    distance: _GeocentricDistance, t: Time, *, is_min: bool
) -> datetime:
    """Return the minute boundary closest to a distance extremum located near t.

    The three minute boundaries around t are compared on the distance itself, so
    the published minute is that of the true extremum whatever the sampling grid of
    the search that produced t.

    Args:
        distance: Geometric Earth-Moon distance function.
        t: Extremum instant located within a few seconds.
        is_min: True for a perigee, False for an apogee.

    Returns:
        Aware UTC datetime on a minute boundary.
    """
    nearest = round_to_minute_utc(t.utc_datetime())
    candidates = [nearest + timedelta(minutes=offset) for offset in (-1, 0, 1)]
    values = distance(t.ts.from_datetimes(candidates))
    return candidates[int(values.argmin() if is_min else values.argmax())]


# -----------------------------------------------------------------------------
# Payloads
# -----------------------------------------------------------------------------


def compute_current_payload(
    eph: SpiceKernel, t: Time, observer: GeographicPosition
) -> dict[str, Any]:
    """Compute the current Moon position and the surrounding moonrises and moonsets.

    Args:
        eph: Loaded ephemeris.
        t: Reference time.
        observer: Observer position on the WGS84 ellipsoid.

    Returns:
        The payload of the main coordinator.
    """
    topocentric, az_deg, alt_deg = _topocentric_apparent(eph, t, observer)
    geocentric = _geocentric_apparent(eph, t)
    ecl_lon_topo, ecl_lat_topo = _ecliptic_lon_lat_deg(topocentric)
    ecl_lon_geo, ecl_lat_geo = _ecliptic_lon_lat_deg(geocentric)
    distance_km = float(_GeocentricDistance(eph)(t))
    payload: dict[str, Any] = {
        KEY_PHASE: _moon_phase_code(float(almanac.moon_phase(eph, t).degrees)),
        KEY_AZIMUTH: round(az_deg, 4),
        KEY_ELEVATION: round(alt_deg, 4),
        KEY_ILLUMINATION: round(
            100.0 * float(geocentric.fraction_illuminated(eph["sun"])), 3
        ),
        KEY_DISTANCE: round(distance_km, 3),
        KEY_PARALLAX: round(_moon_parallax_angle_deg(distance_km), 4),
        KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC: round(ecl_lon_topo, 6),
        KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC: round(ecl_lat_topo, 6),
        KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC: round(ecl_lon_geo, 6),
        KEY_ECLIPTIC_LATITUDE_GEOCENTRIC: round(ecl_lat_geo, 6),
        KEY_ABOVE_HORIZON: alt_deg > HORIZON_ALTITUDE_DEG,
        **_zodiac(
            ecl_lon_geo, KEY_ZODIAC_SIGN_CURRENT_MOON, KEY_ZODIAC_DEGREE_CURRENT_MOON
        ),
    }
    times, rising = almanac.find_discrete(
        t - RISE_SET_SEARCH_DAYS,
        t + RISE_SET_SEARCH_DAYS,
        almanac.risings_and_settings(
            eph, eph["moon"], observer, horizon_degrees=HORIZON_ALTITUDE_DEG
        ),
    )
    events = list(zip(times, rising, strict=True))
    for keys, rises in ((_RISE_KEYS, True), (_SET_KEYS, False)):
        surrounding = _surrounding(t, (ti for ti, up in events if bool(up) == rises))
        payload.update(zip(keys, map(_time_to_utc, surrounding), strict=True))
    return payload


def compute_events_payload(eph: SpiceKernel, t: Time, tz: tzinfo) -> dict[str, Any]:
    """Compute the values that only change at astronomical events.

    Args:
        eph: Loaded ephemeris.
        t: Reference time.
        tz: Time zone defining the calendar months that name full moons.

    Returns:
        The previous and next principal phases, full moon names and apsides, and the
        ecliptic and zodiac position of the Moon at the surrounding new and full moons.
    """
    times, phases = almanac.find_discrete(
        t - PHASE_SEARCH_DAYS_BACK,
        t + PHASE_SEARCH_DAYS_AHEAD,
        almanac.moon_phases(eph),
    )
    phase_events = list(zip(times, phases, strict=True))
    lunations: dict[str, Time | None] = {}
    for phase, keys in _PHASE_KEYS.items():
        surrounding = _surrounding(t, (ti for ti, p in phase_events if p == phase))
        lunations.update(zip(keys, surrounding, strict=True))
    payload: dict[str, Any] = {key: _time_to_utc(ti) for key, ti in lunations.items()}

    names = _full_moon_name_codes(
        (ti for ti, p in phase_events if p == FULL_MOON), tz
    )
    for instant_key, name_key, alt_names_key in _FULL_MOON_NAME_KEYS:
        name = names.get(payload[instant_key])
        payload[name_key] = name
        payload[alt_names_key] = _full_moon_alt_names_state_code(name)

    distance = _GeocentricDistance(eph)
    for keys, search, is_min in (
        (_PERIGEE_KEYS, find_minima, True),
        (_APOGEE_KEYS, find_maxima, False),
    ):
        apsides, _distances = search(
            t - APSIS_SEARCH_DAYS,
            t + APSIS_SEARCH_DAYS,
            distance,
            epsilon=APSIS_EPSILON_DAYS,
        )
        payload.update(
            zip(
                keys,
                (
                    None if ti is None else _snap_to_minute(distance, ti, is_min=is_min)
                    for ti in _surrounding(t, apsides)
                ),
                strict=True,
            )
        )

    for instant_key, (lon_key, lat_key, sign_key, degree_key) in _LUNATION_KEYS.items():
        lon, lat = (
            (None, None)
            if (ti := lunations[instant_key]) is None
            else _ecliptic_lon_lat_deg(_geocentric_apparent(eph, ti))
        )
        payload[lon_key] = _round_or_none(lon, 6)
        payload[lat_key] = _round_or_none(lat, 6)
        payload.update(_zodiac(lon, sign_key, degree_key))
    return payload
