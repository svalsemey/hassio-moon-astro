"""Data coordinator for Moon Astro."""

from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Any

import numpy as np
from zoneinfo import ZoneInfo
from timezonefinder import TimezoneFinder

from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant
from homeassistant.helpers.update_coordinator import DataUpdateCoordinator, UpdateFailed

from skyfield.api import Loader, wgs84
from skyfield import almanac

from .const import (
    DE440_FILE,
    DARK_MOON,
    FIRST_QUARTER,
    FULL_MOON,
    LAST_QUARTER,
    CONF_LAT,
    CONF_LON,
    CONF_ALT,
    CONF_USE_HA_TZ,
    KEY_PHASE,
    KEY_AZ,
    KEY_EL,
    KEY_ILLUM,
    KEY_DISTANCE,
    KEY_PARALLAX,
    KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
    KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
    KEY_ECLIPTIC_LONGITUDE_GEOCENTRIC,
    KEY_ECLIPTIC_LATITUDE_GEOCENTRIC,
    KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON,
    KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON,
    KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON,
    KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON,
    KEY_NEXT_RISE,
    KEY_NEXT_SET,
    KEY_NEXT_APOGEE,
    KEY_NEXT_PERIGEE,
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_LAST_QUARTER,
    KEY_ABOVE_HORIZON,
    KEY_WAXING,
    KEY_ZODIAC_SIGN_NEXT_NEW_MOON,
    KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
    KEY_ZODIAC_DEGREE_NEXT_NEW_MOON,
    KEY_ZODIAC_DEGREE_NEXT_FULL_MOON,
    KEY_ZODIAC_ICON_NEXT_NEW_MOON,
    KEY_ZODIAC_ICON_NEXT_FULL_MOON,
    NEW_MOON_STRICT_PCT,
    FULL_MOON_STRICT_PCT,
    QUARTER_TOL_PCT,
)

_LOGGER = logging.getLogger(__name__)


def _detect_timezone(lat: float, lon: float):
    tf = TimezoneFinder()
    tzname = tf.timezone_at(lat=lat, lng=lon)
    try:
        return ZoneInfo(tzname) if tzname else ZoneInfo("UTC")
    except Exception:
        return ZoneInfo("UTC")


def _tz_for_hass(hass: HomeAssistant):
    tzname = hass.config.time_zone
    try:
        return ZoneInfo(tzname) if tzname else ZoneInfo("UTC")
    except Exception:
        return ZoneInfo("UTC")


def _to_local_iso(dt_utc, tz) -> str | None:
    if dt_utc is None:
        return None
    if dt_utc.tzinfo is None:
        dt_utc = dt_utc.replace(tzinfo=timezone.utc)
    return dt_utc.astimezone(tz).isoformat()


def _safe_time_iso(t_obj, tz):
    if t_obj is None:
        return None
    return _to_local_iso(t_obj.utc_datetime().replace(tzinfo=timezone.utc), tz)


def _moon_illumination_percentage(eph, t):
    frac = almanac.fraction_illuminated(eph, "moon", t)
    return float(frac * 100.0)


def _moon_phase_name(eph, t, ts=None):
    illum_now = float(almanac.fraction_illuminated(eph, "moon", t) * 100.0)

    if ts is not None:
        t_future = ts.tt_jd(t.tt + 6.0 / 24.0)
        illum_future = float(almanac.fraction_illuminated(eph, "moon", t_future) * 100.0)
        waxing = illum_future > illum_now + 1e-6
    else:
        waxing = True

    f = almanac.moon_phases(eph)
    t0 = t - 30.0
    t1 = t + 30.0
    times, phases = almanac.find_discrete(t0, t1, f)

    last_ev = None
    next_ev = None
    for ti, pv in zip(times, phases):
        if ti.tt <= t.tt:
            last_ev = (ti, int(pv))
        elif next_ev is None and ti.tt > t.tt:
            next_ev = (ti, int(pv))

    def close_to(pct, target, tol=QUARTER_TOL_PCT):
        return abs(pct - target) <= tol

    if last_ev is not None:
        _, phase_value = last_ev

        if phase_value == DARK_MOON:
            if waxing and illum_now <= NEW_MOON_STRICT_PCT:
                return "new_moon"

        elif phase_value == FIRST_QUARTER:
            if close_to(illum_now, 50.0) and waxing:
                return "first_quarter"

        elif phase_value == FULL_MOON:
            if (not waxing) or illum_now >= FULL_MOON_STRICT_PCT:
                return "full_moon"

        elif phase_value == LAST_QUARTER:
            if close_to(illum_now, 50.0) and (not waxing):
                return "last_quarter"

    if illum_now <= NEW_MOON_STRICT_PCT:
        return "new_moon" if waxing else "waning_crescent"
    if illum_now <= 45.0:
        return "waxing_crescent" if waxing else "waning_crescent"
    if 45.0 < illum_now < 55.0:
        return "first_quarter" if waxing else "last_quarter"
    if illum_now < FULL_MOON_STRICT_PCT:
        return "waxing_gibbous" if waxing else "waning_gibbous"
    return "full_moon"


def _topocentric_vectors(eph, t, lat, lon, alt_m):
    earth = eph["earth"]
    observer = wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt_m)
    ap = earth + observer
    apparent = ap.at(t).observe(eph["moon"]).apparent()
    alt, az, distance = apparent.altaz()
    return apparent, float(az.degrees), float(alt.degrees), float(distance.km)


def _geocentric_vector(eph, t):
    earth = eph["earth"]
    apparent = earth.at(t).observe(eph["moon"]).apparent()
    return apparent


# ---------- High-accuracy nutation (IAU 1980: full 106-term) and ecliptic-of-date ----------

def _julian_centuries_TT_from_tt(tt: float) -> float:
    # tt = Julian Date (TT). Skyfield Time.tt is JD(TT)
    return (tt - 2451545.0) / 36525.0

def _deg_to_rad(x: float) -> float:
    return x * math.pi / 180.0

def _arcsec_to_rad(x: float) -> float:
    return _deg_to_rad(x / 3600.0)

def _mean_obliquity_arcsec(T: float) -> float:
    # IAU 2006 polynomial for mean obliquity; accurate for centuries near J2000
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

# Fundamental Delaunay arguments (degrees) — consistent with IAU 1980 nutation
def _fundamental_arguments_deg(T: float):
    # Mean anomaly of the Moon (M'), of the Sun (M), Moon's argument of latitude (F),
    # Moon's elongation from the Sun (D), and longitude of the ascending node (Ω).
    Lm = (134.96298139 + (1325.0 * 360.0 + 198.8673981) * T
          + 0.0086972 * T**2 + T**3 / 56250.0)
    Ls = (357.52772333 + (99.0 * 360.0 + 359.0503400) * T
          - 0.0001603 * T**2 - T**3 / 300000.0)
    F  = (93.27191028 + (1342.0 * 360.0 + 82.0175381) * T
          - 0.0036825 * T**2 + T**3 / 327270.0)
    D  = (297.85036306 + (1236.0 * 360.0 + 307.1114800) * T
          - 0.0019142 * T**2 + T**3 / 189474.0)
    Om = (125.04452222 - (5.0 * 360.0 + 134.1362608) * T
          + 0.0020708 * T**2 + T**3 / 450000.0)

    def norm(x): return (x % 360.0 + 360.0) % 360.0
    return norm(Lm), norm(Ls), norm(F), norm(D), norm(Om)

# Full IAU 1980 106-term nutation series
# Format: (D, M, M', F, Ω, Δψ_arcsec, Δψ_t_arcsec_per_century, Δε_arcsec, Δε_t_arcsec_per_century)
_IAU1980_TERMS = [
    ( 0,  0,  0,  0,  1, -171996.0, -174.2,  92025.0,    8.9),
    ( 0,  0,  2, -2,  2,  -13187.0,   -1.6,   5736.0,   -3.1),
    ( 0,  0,  2,  0,  2,   -2274.0,   -0.2,    977.0,   -0.5),
    ( 0,  0,  0,  0,  2,    2062.0,    0.2,   -895.0,    0.5),
    ( 0,  1,  0,  0,  0,    1426.0,   -3.4,     54.0,   -0.1),
    ( 1,  0,  0,  0,  0,     712.0,    0.1,     -7.0,    0.0),
    ( 0,  1,  2, -2,  2,    -517.0,    1.2,    224.0,   -0.6),
    ( 0,  0,  2,  0,  1,    -386.0,   -0.4,    200.0,    0.0),
    ( 1,  0,  2,  0,  2,    -301.0,    0.0,    129.0,   -0.1),
    ( 0, -1,  2, -2,  2,     217.0,   -0.5,    -95.0,    0.3),
    ( 1,  0,  0, -2,  0,    -158.0,    0.0,      0.0,    0.0),
    ( 0,  0,  2, -2,  1,     129.0,    0.1,    -70.0,    0.0),
    (-1,  0,  2,  0,  2,     123.0,    0.0,    -53.0,    0.0),
    ( 0,  0,  0,  2,  0,      63.0,    0.0,      0.0,    0.0),
    ( 1,  0,  2, -2,  2,      63.0,    0.1,    -33.0,    0.0),
    (-1,  0,  0,  2,  0,     -58.0,   -0.1,      0.0,    0.0),
    (-1,  0,  2,  2,  2,     -51.0,    0.0,     27.0,    0.0),
    ( 1,  0,  2,  0,  1,      48.0,    0.0,    -24.0,    0.0),
    ( 0,  0,  2,  2,  2,     -38.0,    0.0,     16.0,    0.0),
    ( 2,  0,  2,  0,  2,     -31.0,    0.0,     13.0,    0.0),
    ( 2,  0,  0,  0,  0,      29.0,    0.0,      0.0,    0.0),
    ( 0,  0,  2,  0,  0,      29.0,    0.0,      0.0,    0.0),
    ( 0,  0,  2, -2,  0,      26.0,    0.0,      0.0,    0.0),
    (-1,  0,  2,  0,  1,      21.0,    0.0,    -10.0,    0.0),
    ( 0,  2,  0,  0,  0,     -16.0,    0.0,      0.0,    0.0),
    ( 1,  0,  0,  0,  1,      16.0,    0.0,     -8.0,    0.0),
    ( 0,  0,  0,  0,  3,     -15.0,    0.0,      9.0,    0.0),
    ( 1,  0,  2, -2,  1,     -13.0,    0.0,      7.0,    0.0),
    ( 0,  1,  0,  0,  1,     -12.0,    0.0,      6.0,    0.0),
    (-1,  0,  0,  0,  1,      11.0,    0.0,     -5.0,    0.0),
    ( 0,  1,  2,  0,  2,     -10.0,    0.0,      5.0,    0.0),
    ( 0, -1,  2,  0,  2,      -8.0,    0.0,      3.0,    0.0),
    ( 2,  0,  2, -2,  2,      -7.0,    0.0,      3.0,    0.0),
    ( 1,  1,  0,  0,  0,      -7.0,    0.0,      0.0,    0.0),
    (-1,  1,  0,  0,  0,      -7.0,    0.0,      0.0,    0.0),
    ( 0,  1,  2, -2,  2,      -7.0,    0.0,      3.0,    0.0),
    ( 0,  0,  0,  2,  1,       6.0,    0.0,     -3.0,    0.0),
    ( 1,  0,  2,  2,  2,      -6.0,    0.0,      3.0,    0.0),
    ( 1,  0,  0,  2,  0,       6.0,    0.0,      0.0,    0.0),
    ( 2,  0,  2,  0,  1,      -6.0,    0.0,      3.0,    0.0),
    ( 0,  0,  0,  2,  0,      -5.0,    0.0,      0.0,    0.0),
    ( 0, -1,  2, -2,  2,       5.0,    0.0,     -3.0,    0.0),
    ( 2,  0,  2, -2,  1,       5.0,    0.0,     -3.0,    0.0),
    ( 0,  1,  0,  0,  2,      -5.0,    0.0,      3.0,    0.0),
    ( 1,  0,  2, -2,  0,      -4.0,    0.0,      0.0,    0.0),
    ( 0,  0,  0,  1,  0,      -4.0,    0.0,      0.0,    0.0),
    ( 1,  1,  0,  0,  1,      -4.0,    0.0,      2.0,    0.0),
    ( 1, -1,  0,  0,  1,      -4.0,    0.0,      2.0,    0.0),
    ( 1,  0,  0, -1,  0,      -4.0,    0.0,      0.0,    0.0),
    ( 0,  0,  2,  1,  2,      -4.0,    0.0,      2.0,    0.0),
    ( 1,  0,  0,  1,  0,       3.0,    0.0,      0.0,    0.0),
    ( 1, -1,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0),
    ( 0, -1,  2,  0,  1,      -3.0,    0.0,      1.0,    0.0),
    ( 1,  1,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0),
    (-1,  1,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0),
    ( 3,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0),
    ( 0,  0,  0,  0,  1,       3.0,    0.0,     -1.0,    0.0),
    (-1,  0,  2,  2,  1,      -3.0,    0.0,      1.0,    0.0),
    ( 0,  0,  2,  2,  1,      -3.0,    0.0,      1.0,    0.0),
    ( 1,  0,  2,  2,  1,      -3.0,    0.0,      1.0,    0.0),
    (-1,  0,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0),
    ( 2,  0,  0,  0,  1,       2.0,    0.0,     -1.0,    0.0),
    ( 1,  0,  0,  0,  2,      -2.0,    0.0,      1.0,    0.0),
    ( 2,  0,  2, -2,  2,      -2.0,    0.0,      1.0,    0.0),
    ( 0, -1,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0),
    ( 0,  1,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0),
    ( 0, -2,  0,  2,  0,      -2.0,    0.0,      0.0,    0.0),
    ( 2,  0,  0, -2,  1,       2.0,    0.0,     -1.0,    0.0),
    (-2,  0,  2,  0,  1,       2.0,    0.0,     -1.0,    0.0),
    ( 0,  0,  2,  0,  1,       2.0,    0.0,     -1.0,    0.0),
    ( 2,  0,  2,  0,  1,       2.0,    0.0,     -1.0,    0.0),
    ( 0,  0,  0,  2,  1,       2.0,    0.0,     -1.0,    0.0),
    ( 1,  0,  2, -2,  2,      -1.0,    0.0,      0.0,    0.0),
    ( 1,  0,  0,  0,  0,      -1.0,    0.0,      0.0,    0.0),
    (-1,  0,  0,  0,  2,       1.0,    0.0,      0.0,    0.0),
    ( 1,  0,  0, -2,  1,       1.0,    0.0,      0.0,    0.0),
    ( 0,  0,  0,  2,  2,      -1.0,    0.0,      0.0,    0.0),
    ( 0,  0,  2,  2,  2,      -1.0,    0.0,      0.0,    0.0),
    ( 1,  0,  2,  0,  2,      -1.0,    0.0,      0.0,    0.0),
    ( 0,  0,  2,  0,  2,      -1.0,    0.0,      0.0,    0.0),
    ( 1,  0,  0,  2,  0,       1.0,    0.0,      0.0,    0.0),
    ( 0,  0,  0,  2,  1,      -1.0,    0.0,      0.0,    0.0),
    ( 1,  0,  2, -2,  1,      -1.0,    0.0,      0.0,    0.0),
    ( 1,  1,  0, -2,  0,      -1.0,    0.0,      0.0,    0.0),
    ( 1, -1,  0, -2,  0,      -1.0,    0.0,      0.0,    0.0),
    ( 2,  0,  0,  0,  0,      -1.0,    0.0,      0.0,    0.0),
    ( 0,  1,  2,  0,  1,      -1.0,    0.0,      0.0,    0.0),
    (-1,  0,  2,  2,  2,      -1.0,    0.0,      0.0,    0.0),
    ( 0, -1,  2,  2,  2,      -1.0,    0.0,      0.0,    0.0),
    ( 1, -1,  2,  0,  1,      -1.0,    0.0,      0.0,    0.0),
    ( 0,  0,  2, -1,  2,      -1.0,    0.0,      0.0,    0.0),
    ( 1,  0,  0,  0,  1,      -1.0,    0.0,      0.0,    0.0),
    ( 1,  0,  0, -1,  1,      -1.0,    0.0,      0.0,    0.0),
    ( 0,  1,  0,  1,  0,      -1.0,    0.0,      0.0,    0.0),
    ( 0, -1,  0,  1,  0,      -1.0,    0.0,      0.0,    0.0),

    # Additional very small terms with time rates (kept for canonical completeness)
    ( 0,  0,  1,  0,  1,       0.0,    0.1,      0.0,    0.0),
    ( 0,  0,  1,  0,  1,       0.0,   -0.1,      0.0,    0.0),
    ( 0,  2,  2, -2,  2,       0.0,    0.1,      0.0,    0.0),
    ( 0, -2,  2, -2,  2,       0.0,    0.1,      0.0,    0.0),
    ( 2,  0,  0, -2,  0,       0.0,    0.1,      0.0,    0.0),
    ( 2,  0,  2, -2,  1,       0.0,    0.1,      0.0,    0.0),
    ( 2,  0,  2, -2,  1,       0.0,   -0.1,      0.0,    0.0),
]

def _nutation_iau1980(T: float, Lm_deg: float, Ls_deg: float, F_deg: float, D_deg: float, Om_deg: float):
    # Convert arguments to radians
    Lm = _deg_to_rad(Lm_deg)
    Ls = _deg_to_rad(Ls_deg)
    F  = _deg_to_rad(F_deg)
    D  = _deg_to_rad(D_deg)
    Om = _deg_to_rad(Om_deg)

    dpsi_as = 0.0
    deps_as = 0.0
    for (cD, cM, cMp, cF, cOm, ps, ps_t, pe, pe_t) in _IAU1980_TERMS:
        arg = cD * D + cM * Ls + cMp * Lm + cF * F + cOm * Om
        s = math.sin(arg)
        c = math.cos(arg)
        # Linear time dependence per century (IAU 1980 convention)
        dpsi_as += ps * s + ps_t * s * T
        deps_as += pe * c + pe_t * c * T
    return dpsi_as, deps_as  # arcseconds

def _true_obliquity_rad(tt: float) -> float:
    T = _julian_centuries_TT_from_tt(tt)
    Lm, Ls, F, D, Om = _fundamental_arguments_deg(T)
    _dpsi_as, deps_as = _nutation_iau1980(T, Lm, Ls, F, D, Om)
    eps0_as = _mean_obliquity_arcsec(T)
    eps_true = _arcsec_to_rad(eps0_as + deps_as)
    return eps_true

def _ecliptic_lon_lat_deg_of_date(apparent_vector):
    """Apparent true ecliptic-of-date longitude/latitude in degrees (tropical)."""
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


def _moon_parallax_angle_deg(distance_km):
    Re_km = 6378.137
    x = Re_km / max(distance_km, 1e-6)
    x = min(1.0, max(0.0, x))
    return math.degrees(math.asin(x))


def _next_rise_set(eph, lat, lon, alt_m, t_start):
    observer = wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt_m)
    f = almanac.risings_and_settings(eph, eph["moon"], observer)
    t0 = t_start
    t1 = t_start + 7.0
    times, events = almanac.find_discrete(t0, t1, f)
    next_rise = None
    next_set = None
    for ti, ei in zip(times, events):
        if ei and next_rise is None and ti.tt > t_start.tt:
            next_rise = ti
        if (not ei) and next_set is None and ti.tt > t_start.tt:
            next_set = ti
        if next_rise is not None and next_set is not None:
            break
    return next_rise, next_set


@dataclass
class _BrentResult:
    tt: float
    fval: float
    iterations: int


def _brent_extremum(f, a, b, is_min=True, tol=1e-6, max_iter=100):
    g = (lambda x: -f(x)) if not is_min else f

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
            parabolic_ok = (abs(p) < abs(0.5 * q * e)) and (p > q * (a - x)) and (p < q * (b - x))
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
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu

    return _BrentResult(tt=x, fval=(fx if is_min else -fx), iterations=max_iter)


def _refine_bracket(t_list, y_list, kind="max"):
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


def _find_next_apogee(eph, ts, t_start):
    earth, moon = eph["earth"], eph["moon"]

    def geocentric_distance_km(t):
        return earth.at(t).observe(moon).distance().km

    days_window = 40.0
    step_hours = 2.0
    steps = int(days_window * 24 / step_hours) + 1
    t_list = []
    d_list = []
    for i in range(steps):
        dt_days = i * (step_hours / 24.0)
        t_i = t_start + dt_days
        t_list.append(t_i)
        d_list.append(geocentric_distance_km(t_i))

    bracket = _refine_bracket(t_list, d_list, kind="max")
    if bracket is None:
        return None
    tt_a, tt_b = bracket

    def f_tt(tt):
        return geocentric_distance_km(ts.tt_jd(tt))

    res = _brent_extremum(f_tt, tt_a, tt_b, is_min=False, tol=1e-7, max_iter=200)
    t_ext = ts.tt_jd(res.tt)
    return t_ext if t_ext.tt > t_start.tt else None


def _find_next_perigee(eph, ts, t_start):
    earth, moon = eph["earth"], eph["moon"]

    def geocentric_distance_km(t):
        return earth.at(t).observe(moon).distance().km

    days_window = 40.0
    step_hours = 2.0
    steps = int(days_window * 24 / step_hours) + 1
    t_list = []
    d_list = []
    for i in range(steps):
        dt_days = i * (step_hours / 24.0)
        t_i = t_start + dt_days
        t_list.append(t_i)
        d_list.append(geocentric_distance_km(t_i))

    bracket = _refine_bracket(t_list, d_list, kind="min")
    if bracket is None:
        return None
    tt_a, tt_b = bracket

    def f_tt(tt):
        return geocentric_distance_km(ts.tt_jd(tt))

    res = _brent_extremum(f_tt, tt_a, tt_b, is_min=True, tol=1e-7, max_iter=200)
    t_ext = ts.tt_jd(res.tt)
    return t_ext if t_ext.tt > t_start.tt else None


def _find_phase_next(ts, f, t_start, phase_value):
    t0 = t_start
    t1 = t_start + 40.0
    times, phases = almanac.find_discrete(t0, t1, f)
    for ti, pv in zip(times, phases):
        if int(pv) == int(phase_value) and ti.tt > t_start.tt:
            return ti
    return None


def _zodiac_sign_from_longitude_deg(lon_deg: float) -> str | None:
    if lon_deg != lon_deg:  # NaN
        return None
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
    return float((lon_deg % 360.0) % 30.0)


def _zodiac_icon(sign: str | None) -> str | None:
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
    def __init__(self, hass: HomeAssistant, lat: float, lon: float, elev: float, interval: timedelta) -> None:
        super().__init__(hass, logger=logging.getLogger(__name__), name="Moon Astro", update_interval=interval)
        self._lat = float(lat)
        self._lon = float(lon)
        self._elev = float(elev)
        self._eph = None
        self._ts = None
        self._tz = None
        self._use_ha_tz = False
        self._hass = hass

    @classmethod
    def from_config_entry(cls, hass: HomeAssistant, entry: ConfigEntry, interval: timedelta):
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

    async def _async_load_ephemeris(self):
        def _load():
            cache_dir = self._hass.config.path(".skyfield")
            os.makedirs(cache_dir, exist_ok=True)
            load = Loader(cache_dir)
            eph = load(DE440_FILE)
            ts = load.timescale()
            return eph, ts

        return await self._hass.async_add_executor_job(_load)

    async def _async_update_data(self) -> dict[str, Any]:
        try:
            if self._eph is None or self._ts is None:
                self._eph, self._ts = await self._async_load_ephemeris()

            now_utc = datetime.now(timezone.utc)
            t = self._ts.from_datetime(now_utc)
            t_future = self._ts.from_datetime(now_utc + timedelta(hours=6))

            def _calc():
                topo_vec, az_deg, el_deg, dist_km = _topocentric_vectors(self._eph, t, self._lat, self._lon, self._elev)
                geo_vec = _geocentric_vector(self._eph, t)

                ecl_lon_topo, ecl_lat_topo = _ecliptic_lon_lat_deg_of_date(topo_vec)
                ecl_lon_geo, ecl_lat_geo = _ecliptic_lon_lat_deg_of_date(geo_vec)

                illum_now = _moon_illumination_percentage(self._eph, t)
                parallax_deg = _moon_parallax_angle_deg(dist_km)

                phase = _moon_phase_name(self._eph, t, self._ts)

                illum_future = _moon_illumination_percentage(self._eph, t_future)
                waxing = illum_future > illum_now + 1e-6

                try:
                    t_rise, t_set = _next_rise_set(self._eph, self._lat, self._lon, self._elev, t)
                except Exception:
                    t_rise, t_set = None, None

                try:
                    f = almanac.moon_phases(self._eph)
                    t_new = _find_phase_next(self._ts, f, t, DARK_MOON)
                    t_first = _find_phase_next(self._ts, f, t, FIRST_QUARTER)
                    t_full = _find_phase_next(self._ts, f, t, FULL_MOON)
                    t_last = _find_phase_next(self._ts, f, t, LAST_QUARTER)
                except Exception:
                    t_new = t_first = t_full = t_last = None

                try:
                    t_apogee = _find_next_apogee(self._eph, self._ts, t)
                except Exception:
                    t_apogee = None
                try:
                    t_perigee = _find_next_perigee(self._eph, self._ts, t)
                except Exception:
                    t_perigee = None

                # Compute ecliptic lon/lat at the lunation instants (geocentric, true-of-date)
                lon_new = lat_new = lon_full = lat_full = None
                try:
                    if t_new is not None:
                        v_new = _geocentric_vector(self._eph, t_new)
                        lon_new, lat_new = _ecliptic_lon_lat_deg_of_date(v_new)
                except Exception as exc:
                    _LOGGER.debug("Failed ecliptic lon/lat at new moon: %r", exc)
                    lon_new = lat_new = None
                try:
                    if t_full is not None:
                        v_full = _geocentric_vector(self._eph, t_full)
                        lon_full, lat_full = _ecliptic_lon_lat_deg_of_date(v_full)
                except Exception as exc:
                    _LOGGER.debug("Failed ecliptic lon/lat at full moon: %r", exc)
                    lon_full = lat_full = None

                zodiac_new = None
                zodiac_full = None
                zodiac_degree_new = None
                zodiac_degree_full = None
                try:
                    if lon_new is not None and lon_new == lon_new:
                        zodiac_new = _zodiac_sign_from_longitude_deg(lon_new)
                        zodiac_degree_new = round(_degree_within_sign(lon_new), 4)
                        _LOGGER.debug(
                            "Next new moon t=%s lon_moon_geo=%.6f° lat_moon_geo=%.6f° sign=%s deg_in_sign=%s",
                            _safe_time_iso(t_new, timezone.utc),
                            lon_new,
                            lat_new,
                            zodiac_new,
                            zodiac_degree_new,
                        )
                except Exception as exc:
                    _LOGGER.debug("Failed zodiac_new computation: %r", exc)
                    zodiac_new = None
                    zodiac_degree_new = None
                try:
                    if lon_full is not None and lon_full == lon_full:
                        zodiac_full = _zodiac_sign_from_longitude_deg(lon_full)
                        zodiac_degree_full = round(_degree_within_sign(lon_full), 4)
                        _LOGGER.debug(
                            "Next full moon t=%s lon_moon_geo=%.6f° lat_moon_geo=%.6f° sign=%s deg_in_sign=%s",
                            _safe_time_iso(t_full, timezone.utc),
                            lon_full,
                            lat_full,
                            zodiac_full,
                            zodiac_degree_full,
                        )
                except Exception as exc:
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
                    KEY_ECLIPTIC_LONGITUDE_NEXT_FULL_MOON: None if lon_full is None or lon_full != lon_full else round(lon_full, 6),
                    KEY_ECLIPTIC_LATITUDE_NEXT_FULL_MOON: None if lat_full is None or lat_full != lat_full else round(lat_full, 6),
                    KEY_ECLIPTIC_LONGITUDE_NEXT_NEW_MOON: None if lon_new is None or lon_new != lon_new else round(lon_new, 6),
                    KEY_ECLIPTIC_LATITUDE_NEXT_NEW_MOON: None if lat_new is None or lat_new != lat_new else round(lat_new, 6),
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
        except Exception as err:
            raise UpdateFailed(str(err)) from err