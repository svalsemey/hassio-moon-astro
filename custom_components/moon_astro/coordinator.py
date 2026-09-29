"""Coordinators for Moon Astro.

Two coordinators share the loaded DE440 kernel: one recomputes the current Moon
position on the configured scan interval, the other recomputes the values that only
change at astronomical events and schedules its own refresh after the next event.
All Skyfield work runs in the executor; instants are published as aware UTC datetimes
rounded to the minute.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta, tzinfo
import logging
import math
from typing import Any

from skyfield import almanac
from skyfield.api import wgs84
from skyfield.jpllib import SpiceKernel
from skyfield.positionlib import Apparent
from skyfield.searchlib import find_maxima, find_minima
from skyfield.timelib import Time, Timescale
from skyfield.toposlib import GeographicPosition

from homeassistant.config_entries import ConfigEntry
from homeassistant.core import CALLBACK_TYPE, HomeAssistant
from homeassistant.helpers.event import async_track_point_in_time
from homeassistant.helpers.update_coordinator import DataUpdateCoordinator, UpdateFailed
from homeassistant.util import dt as dt_util

from .const import (
    APSIS_EPSILON_DAYS,
    APSIS_SEARCH_DAYS,
    APSIS_STEP_DAYS,
    CONF_ALT,
    CONF_LAT,
    CONF_LON,
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

_LOGGER = logging.getLogger(__name__)

# Errors turning a refresh into a failure: a kernel read error, or a Skyfield search or
# numeric failure such as an instant outside the ephemeris range.
_RECOVERABLE_UPDATE_ERRORS: tuple[type[Exception], ...] = (
    OSError,
    ValueError,
    ArithmeticError,
    RuntimeError,
)

_EARTH_EQUATORIAL_RADIUS_KM = 6378.137

# Margin after an event boundary before the event-based values are recomputed.
_BOUNDARY_REFRESH_DELAY = timedelta(minutes=2)

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

# Upcoming event keys driving the boundary refresh of the events coordinator.
_NEXT_EVENT_KEYS: tuple[str, ...] = (
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_LAST_QUARTER,
    KEY_NEXT_APOGEE,
    KEY_NEXT_PERIGEE,
)

# -----------------------------------------------------------------------------
# Time helpers
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
# Positions
# -----------------------------------------------------------------------------


def _topocentric_apparent(
    eph: SpiceKernel, t: Time, observer: GeographicPosition
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


def _geocentric_apparent(eph: SpiceKernel, t: Time) -> Apparent:
    """Return the geocentric apparent Moon position.

    Args:
        eph: Loaded ephemeris.
        t: Skyfield Time.

    Returns:
        Apparent position seen from the geocenter.
    """
    return eph["earth"].at(t).observe(eph["moon"]).apparent()


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


# ---------- High-accuracy nutation (IAU 1980: full 106-term) and ecliptic-of-date ----------

# >>> KEEP THE EXISTING BLOCK UNCHANGED HERE: _julian_centuries_tt_from_tt, _deg_to_rad,
# >>> _arcsec_to_rad, _mean_obliquity_arcsec, _fundamental_arguments_deg, _IAU1980_TERMS,
# >>> _nutation_iau1980, _true_obliquity_rad and _ecliptic_lon_lat_deg_of_date.


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
        names[_round_to_minute_utc(utc)] = (
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


# -----------------------------------------------------------------------------
# Payload computations (executor)
# -----------------------------------------------------------------------------


def _snap_to_minute(
    distance: _GeocentricDistance, ts: Timescale, t: Time, *, is_min: bool
) -> datetime:
    """Return the minute boundary closest to a distance extremum located near t.

    The three minute boundaries around t are compared on the distance itself, so
    the published minute is that of the true extremum whatever the sampling grid of
    the search that produced t.

    Args:
        distance: Geometric Earth-Moon distance function.
        ts: Skyfield timescale.
        t: Extremum instant located within a few seconds.
        is_min: True for a perigee, False for an apogee.

    Returns:
        Aware UTC datetime on a minute boundary.
    """
    nearest = _round_to_minute_utc(t.utc_datetime())
    candidates = [nearest + timedelta(minutes=offset) for offset in (-1, 0, 1)]
    values = distance(ts.from_datetimes(candidates))
    return candidates[int(values.argmin() if is_min else values.argmax())]


def _compute_current_payload(
    eph: SpiceKernel,
    t: Time,
    observer: GeographicPosition,
    distance: _GeocentricDistance,
) -> dict[str, Any]:
    """Compute the current Moon position, its derived quantities and the surrounding
    moonrises and moonsets.

    Args:
        eph: Loaded ephemeris.
        t: Reference time.
        observer: Observer position on the WGS84 ellipsoid.
        distance: Geometric Earth-Moon distance function.

    Returns:
        The payload of the main coordinator.
    """
    topocentric, az_deg, alt_deg = _topocentric_apparent(eph, t, observer)
    geocentric = _geocentric_apparent(eph, t)
    ecl_lon_topo, ecl_lat_topo = _ecliptic_lon_lat_deg_of_date(topocentric)
    ecl_lon_geo, ecl_lat_geo = _ecliptic_lon_lat_deg_of_date(geocentric)
    distance_km = float(distance(t))
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
    for keys, wanted in ((_RISE_KEYS, True), (_SET_KEYS, False)):
        events = (ti for ti, up in zip(times, rising, strict=True) if bool(up) is wanted)
        payload.update(zip(keys, map(_time_to_utc, _surrounding(t, events)), strict=True))
    return payload


def _compute_events_payload(
    eph: SpiceKernel,
    ts: Timescale,
    t: Time,
    tz: tzinfo,
    distance: _GeocentricDistance,
) -> dict[str, Any]:
    """Compute the values that only change at astronomical events.

    Args:
        eph: Loaded ephemeris.
        ts: Skyfield timescale.
        t: Reference time.
        tz: Time zone defining the calendar months that name full moons.
        distance: Geometric Earth-Moon distance function.

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
        lunations.update(
            zip(
                keys,
                _surrounding(t, (ti for ti, p in phase_events if p == phase)),
                strict=True,
            )
        )
    payload: dict[str, Any] = {key: _time_to_utc(ti) for key, ti in lunations.items()}

    names = _full_moon_name_codes(
        (ti for ti, p in phase_events if p == FULL_MOON), tz
    )
    for instant_key, name_key, alt_names_key in _FULL_MOON_NAME_KEYS:
        name = names.get(payload[instant_key])
        payload[name_key] = name
        payload[alt_names_key] = _full_moon_alt_names_state_code(name)

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
                    None if ti is None else _snap_to_minute(distance, ts, ti, is_min=is_min)
                    for ti in _surrounding(t, apsides)
                ),
                strict=True,
            )
        )

    for instant_key, (lon_key, lat_key, sign_key, degree_key) in _LUNATION_KEYS.items():
        lon, lat = (
            (None, None)
            if (ti := lunations[instant_key]) is None
            else _ecliptic_lon_lat_deg_of_date(_geocentric_apparent(eph, ti))
        )
        payload[lon_key] = _round_or_none(lon, 6)
        payload[lat_key] = _round_or_none(lat, 6)
        payload.update(_zodiac(lon, sign_key, degree_key))
    return payload


# -----------------------------------------------------------------------------
# Coordinators
# -----------------------------------------------------------------------------


class _MoonAstroBaseCoordinator(DataUpdateCoordinator[dict[str, Any]]):
    """Coordinator running a Skyfield computation in the executor at each refresh."""

    def __init__(
        self,
        hass: HomeAssistant,
        entry: MoonAstroConfigEntry,
        *,
        name: str,
        eph: SpiceKernel,
        ts: Timescale,
        interval: timedelta,
    ) -> None:
        """Initialize the coordinator.

        Args:
            hass: Home Assistant instance.
            entry: Config entry owning the coordinator.
            name: Coordinator name used in logs.
            eph: Loaded ephemeris shared by the coordinators.
            ts: Loaded timescale shared by the coordinators.
            interval: Update interval.
        """
        super().__init__(
            hass,
            logger=_LOGGER,
            config_entry=entry,
            name=name,
            update_interval=interval,
        )
        self._eph = eph
        self._ts = ts
        self._distance = _GeocentricDistance(eph)

    def _compute_payload(self, t: Time) -> dict[str, Any]:
        """Compute the payload for the reference time t; runs in the executor."""
        raise NotImplementedError

    async def _async_update_data(self) -> dict[str, Any]:
        """Compute the payload for the current minute.

        Returns:
            The payload exposed by the entities.

        Raises:
            UpdateFailed: If the computation fails.
        """
        t = self._ts.from_datetime(_round_to_minute_utc(dt_util.utcnow()))
        try:
            return await self.hass.async_add_executor_job(self._compute_payload, t)
        except _RECOVERABLE_UPDATE_ERRORS as err:
            raise UpdateFailed(f"{self.name} computation failed: {err}") from err


class MoonAstroCoordinator(_MoonAstroBaseCoordinator):
    """Coordinator computing the current Moon position and the surrounding rises and sets."""

    def __init__(
        self,
        hass: HomeAssistant,
        entry: MoonAstroConfigEntry,
        *,
        eph: SpiceKernel,
        ts: Timescale,
        interval: timedelta,
    ) -> None:
        """Initialize the coordinator with the observer location of the entry.

        Args:
            hass: Home Assistant instance.
            entry: Config entry providing the observer coordinates.
            eph: Loaded ephemeris shared by the coordinators.
            ts: Loaded timescale shared by the coordinators.
            interval: Update interval.
        """
        super().__init__(
            hass, entry, name="Moon Astro", eph=eph, ts=ts, interval=interval
        )
        self._observer = wgs84.latlon(
            latitude_degrees=float(entry.data.get(CONF_LAT, hass.config.latitude)),
            longitude_degrees=float(entry.data.get(CONF_LON, hass.config.longitude)),
            elevation_m=float(entry.data.get(CONF_ALT, hass.config.elevation)),
        )

    def _compute_payload(self, t: Time) -> dict[str, Any]:
        """Compute the current position payload; runs in the executor."""
        return _compute_current_payload(self._eph, t, self._observer, self._distance)


class MoonAstroEventsCoordinator(_MoonAstroBaseCoordinator):
    """Coordinator computing the values that only change at astronomical events.

    Besides the periodic fallback refresh, a refresh is scheduled shortly after the
    earliest upcoming event so that the previous and next values switch on time.
    """

    def __init__(
        self,
        hass: HomeAssistant,
        entry: MoonAstroConfigEntry,
        *,
        eph: SpiceKernel,
        ts: Timescale,
        interval: timedelta,
        tz: tzinfo,
    ) -> None:
        """Initialize the coordinator.

        Args:
            hass: Home Assistant instance.
            entry: Config entry owning the coordinator.
            eph: Loaded ephemeris shared by the coordinators.
            ts: Loaded timescale shared by the coordinators.
            interval: Fallback update interval.
            tz: Time zone defining the calendar months that name full moons.
        """
        super().__init__(
            hass, entry, name="Moon Astro Events", eph=eph, ts=ts, interval=interval
        )
        self._tz = tz
        self._unsub_boundary: CALLBACK_TYPE | None = None
        self.next_refresh_utc: datetime | None = None

    def _compute_payload(self, t: Time) -> dict[str, Any]:
        """Compute the event-based payload; runs in the executor."""
        return _compute_events_payload(self._eph, self._ts, t, self._tz, self._distance)

    async def _async_update_data(self) -> dict[str, Any]:
        """Compute the payload and schedule the refresh following the next event."""
        data = await super()._async_update_data()
        self._schedule_boundary_refresh(data)
        return data

    def _schedule_boundary_refresh(self, data: dict[str, Any]) -> None:
        """Schedule a refresh shortly after the earliest upcoming event.

        Also records the instant of the next refresh, whether it comes from the event
        boundary or from the periodic fallback interval.

        Args:
            data: Freshly computed payload.
        """
        self._cancel_boundary_refresh()
        fallback = (
            None
            if (interval := self.update_interval) is None
            else dt_util.utcnow() + interval
        )
        upcoming = [dt for key in _NEXT_EVENT_KEYS if (dt := data[key]) is not None]
        boundary = None if not upcoming else min(upcoming) + _BOUNDARY_REFRESH_DELAY
        self.next_refresh_utc = min(
            (dt for dt in (fallback, boundary) if dt is not None), default=None
        )
        if boundary is not None:
            self._unsub_boundary = async_track_point_in_time(
                self.hass, self._async_refresh_at_boundary, boundary
            )

    async def _async_refresh_at_boundary(self, _now: datetime) -> None:
        """Request a refresh once the event boundary has passed."""
        self._unsub_boundary = None
        await self.async_request_refresh()

    def _cancel_boundary_refresh(self) -> None:
        """Cancel the pending boundary refresh, if any."""
        if self._unsub_boundary is not None:
            self._unsub_boundary()
            self._unsub_boundary = None

    async def async_shutdown(self) -> None:
        """Cancel the boundary refresh and stop the coordinator."""
        self._cancel_boundary_refresh()
        await super().async_shutdown()


# -----------------------------------------------------------------------------
# Config entry runtime data
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MoonAstroRuntimeData:
    """Runtime objects attached to a loaded config entry."""

    coordinator: MoonAstroCoordinator
    events_coordinator: MoonAstroEventsCoordinator


type MoonAstroConfigEntry = ConfigEntry[MoonAstroRuntimeData]


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
