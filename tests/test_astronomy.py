"""Tests for the Moon Astro computations, run on the DE440 excerpt."""

from datetime import UTC, datetime, timedelta, tzinfo
from typing import Any
from unittest.mock import patch
from zoneinfo import ZoneInfo

from custom_components.moon_astro.astronomy import (
    _full_moon_name_codes,
    _moon_phase_code,
    _zodiac,
    compute_events_payload,
    round_to_minute_utc,
)
from custom_components.moon_astro.const import (
    HORIZON_ALTITUDE_DEG,
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
    KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
    KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON,
)
from custom_components.moon_astro.utils import EphemerisKernel
import numpy as np
import pytest
from skyfield import almanac
from skyfield.api import wgs84
from skyfield.toposlib import GeographicPosition

# Paris Observatory, the location of the USNO reference moonrises and moonsets.
PARIS = wgs84.latlon(48.8363, 2.3366, elevation_m=67.0)

# Locations from the equator to the poles, where the Moon can stay above or below the
# horizon for about two weeks.
LOCATIONS = [
    pytest.param(PARIS, id="paris"),
    pytest.param(wgs84.latlon(78.2232, 15.6267), id="svalbard"),
    pytest.param(wgs84.latlon(-33.8688, 151.2093), id="sydney"),
    pytest.param(wgs84.latlon(-0.1807, -78.4678), id="quito"),
    pytest.param(wgs84.latlon(89.5, 0.0), id="arctic"),
    pytest.param(wgs84.latlon(-89.5, 0.0), id="antarctic"),
]

PREVIOUS_AND_NEXT_KEYS = [
    (KEY_PREVIOUS_RISE, KEY_NEXT_RISE),
    (KEY_PREVIOUS_SET, KEY_NEXT_SET),
    (KEY_PREVIOUS_NEW_MOON, KEY_NEXT_NEW_MOON),
    (KEY_PREVIOUS_FIRST_QUARTER, KEY_NEXT_FIRST_QUARTER),
    (KEY_PREVIOUS_FULL_MOON, KEY_NEXT_FULL_MOON),
    (KEY_PREVIOUS_LAST_QUARTER, KEY_NEXT_LAST_QUARTER),
    (KEY_PREVIOUS_PERIGEE, KEY_NEXT_PERIGEE),
    (KEY_PREVIOUS_APOGEE, KEY_NEXT_APOGEE),
]


def _utc(text: str) -> datetime:
    """Return the aware UTC datetime of an ISO 8601 date and time without offset."""
    return datetime.fromisoformat(text).replace(tzinfo=UTC)


def _events(
    ephemeris: EphemerisKernel,
    reference: datetime,
    observer: GeographicPosition = PARIS,
    tz: tzinfo = UTC,
) -> dict[str, Any]:
    """Compute the event-based payload for a reference instant."""
    eph, ts = ephemeris
    return compute_events_payload(eph, ts.from_datetime(reference), observer, tz)


@pytest.mark.parametrize(
    ("dt", "expected"),
    [
        pytest.param(
            datetime(2026, 11, 2, 22, 52, 29, 999999, tzinfo=UTC),
            _utc("2026-11-02T22:52"),
            id="down",
        ),
        pytest.param(
            datetime(2026, 11, 2, 22, 52, 30, tzinfo=UTC),
            _utc("2026-11-02T22:53"),
            id="up",
        ),
        pytest.param(
            datetime(2026, 12, 31, 23, 59, 45, tzinfo=ZoneInfo("Europe/Paris")),
            _utc("2026-12-31T23:00"),
            id="to_utc",
        ),
    ],
)
def test_round_to_minute_utc(dt: datetime, expected: datetime) -> None:
    """Test the rounding of instants to the nearest UTC minute."""
    assert round_to_minute_utc(dt) == expected


@pytest.mark.parametrize(
    ("phase_angle", "expected"),
    [
        (0.0, "new_moon"),
        (3.0, "new_moon"),
        (3.1, "waxing_crescent"),
        (86.9, "waxing_crescent"),
        (87.0, "first_quarter"),
        (93.0, "first_quarter"),
        (93.1, "waxing_gibbous"),
        (180.0, "full_moon"),
        (183.1, "waning_gibbous"),
        (270.0, "last_quarter"),
        (356.9, "waning_crescent"),
        (357.0, "new_moon"),
        (360.0, "new_moon"),
    ],
)
def test_moon_phase_code(phase_angle: float, expected: str) -> None:
    """Test that principal phases are named within 3 degrees of phase angle."""
    assert _moon_phase_code(phase_angle) == expected


@pytest.mark.parametrize(
    ("longitude", "expected"),
    [
        (0.0, {"sign": "aries", "degree": 0.0}),
        (30.0, {"sign": "taurus", "degree": 0.0}),
        (45.5, {"sign": "taurus", "degree": 15.5}),
        (359.5, {"sign": "pisces", "degree": 29.5}),
        (None, {"sign": None, "degree": None}),
    ],
)
def test_zodiac(longitude: float | None, expected: dict[str, Any]) -> None:
    """Test the zodiac sign and degree of an ecliptic longitude."""
    assert _zodiac(longitude, "sign", "degree") == expected


@pytest.mark.parametrize(
    ("reference", "expected"),
    [
        pytest.param(
            "2026-09-01",
            {
                KEY_NEXT_LAST_QUARTER: "2026-09-04T07:51",
                KEY_NEXT_NEW_MOON: "2026-09-11T03:27",
                KEY_NEXT_FIRST_QUARTER: "2026-09-18T20:44",
                KEY_NEXT_FULL_MOON: "2026-09-26T16:49",
            },
            id="2026-09",
        ),
        pytest.param(
            "2026-10-01",
            {
                KEY_NEXT_LAST_QUARTER: "2026-10-03T13:25",
                KEY_NEXT_NEW_MOON: "2026-10-10T15:50",
                KEY_NEXT_FIRST_QUARTER: "2026-10-18T16:12",
                KEY_NEXT_FULL_MOON: "2026-10-26T04:12",
            },
            id="2026-10",
        ),
        pytest.param(
            "2026-10-31",
            {
                KEY_NEXT_LAST_QUARTER: "2026-11-01T20:28",
                KEY_NEXT_NEW_MOON: "2026-11-09T07:02",
                KEY_NEXT_FIRST_QUARTER: "2026-11-17T11:48",
                KEY_NEXT_FULL_MOON: "2026-11-24T14:53",
            },
            id="2026-11",
        ),
        pytest.param(
            "2026-11-30",
            {
                KEY_NEXT_LAST_QUARTER: "2026-12-01T06:08",
                KEY_NEXT_NEW_MOON: "2026-12-09T00:52",
                KEY_NEXT_FIRST_QUARTER: "2026-12-17T05:42",
                KEY_NEXT_FULL_MOON: "2026-12-24T01:28",
            },
            id="2026-12",
        ),
        pytest.param(
            "2026-12-29",
            {
                KEY_NEXT_LAST_QUARTER: "2026-12-30T18:59",
                KEY_NEXT_NEW_MOON: "2027-01-07T20:24",
                KEY_NEXT_FIRST_QUARTER: "2027-01-15T20:34",
                KEY_NEXT_FULL_MOON: "2027-01-22T12:17",
            },
            id="2027-01",
        ),
    ],
)
def test_phases_match_usno(
    ephemeris: EphemerisKernel, reference: str, expected: dict[str, str]
) -> None:
    """Test the phase instants against those published by the USNO.

    The USNO computes the instants with its own model and rounds them to the minute:
    https://aa.usno.navy.mil/api/moon/phases/date?date=2026-09-01&nump=20
    """
    payload = _events(ephemeris, _utc(reference))

    for key, instant in expected.items():
        assert abs(payload[key] - _utc(instant)) <= timedelta(minutes=1), key


@pytest.mark.parametrize(
    ("reference", "expected"),
    [
        pytest.param(
            "2026-11-02T00:00",
            {KEY_NEXT_SET: "2026-11-02T13:45", KEY_NEXT_RISE: "2026-11-02T23:50"},
            id="2026-11-02",
        ),
        pytest.param(
            "2026-11-10T00:00",
            {KEY_NEXT_RISE: "2026-11-10T08:29", KEY_NEXT_SET: "2026-11-10T16:20"},
            id="2026-11-10",
        ),
        pytest.param(
            "2026-11-20T00:00",
            {KEY_NEXT_SET: "2026-11-20T01:48", KEY_NEXT_RISE: "2026-11-20T13:43"},
            id="2026-11-20",
        ),
    ],
)
def test_rises_and_sets_match_usno(
    ephemeris: EphemerisKernel, reference: str, expected: dict[str, str]
) -> None:
    """Test the moonrises and moonsets in Paris against those published by the USNO.

    The USNO times refer to the upper limb of the Moon, those of Moon Astro to its
    center: in Paris, the two differ by one to two minutes.
    https://aa.usno.navy.mil/api/rstt/oneday?date=2026-11-02&coords=48.8363,2.3366&tz=0
    """
    payload = _events(ephemeris, _utc(reference))

    for key, instant in expected.items():
        assert abs(payload[key] - _utc(instant)) <= timedelta(minutes=3), key


@pytest.mark.parametrize("observer", LOCATIONS)
def test_events_surround_reference_time(
    ephemeris: EphemerisKernel, now: datetime, observer: GeographicPosition
) -> None:
    """Test that every kind of event is found before and after the reference time."""
    payload = _events(ephemeris, now, observer)

    for previous_key, next_key in PREVIOUS_AND_NEXT_KEYS:
        assert payload[previous_key] <= now <= payload[next_key], next_key


@pytest.mark.parametrize("observer", LOCATIONS)
def test_rises_and_sets_cross_the_horizon(
    ephemeris: EphemerisKernel, now: datetime, observer: GeographicPosition
) -> None:
    """Test that the Moon center crosses the horizon within the published minute."""
    eph, ts = ephemeris
    payload = _events(ephemeris, now, observer)
    site = eph["earth"] + observer

    for key, rising in (
        (KEY_PREVIOUS_RISE, True),
        (KEY_NEXT_RISE, True),
        (KEY_PREVIOUS_SET, False),
        (KEY_NEXT_SET, False),
    ):
        around = ts.from_datetimes(
            [payload[key] - timedelta(minutes=1), payload[key] + timedelta(minutes=1)]
        )
        altitude = site.at(around).observe(eph["moon"]).apparent().altaz()[0].degrees
        assert list(altitude > HORIZON_ALTITUDE_DEG) == [not rising, rising], key


def test_rises_and_sets_beyond_search_span(
    ephemeris: EphemerisKernel, now: datetime
) -> None:
    """Test that moonrises and moonsets beyond the widest search span are None."""
    with patch("custom_components.moon_astro.astronomy.RISE_SET_SEARCH_DAYS", (0.01,)):
        payload = _events(ephemeris, now)

    assert [
        payload[key]
        for key in (KEY_PREVIOUS_RISE, KEY_NEXT_RISE, KEY_PREVIOUS_SET, KEY_NEXT_SET)
    ] == [None] * 4


def test_full_moons_beyond_search_window(
    ephemeris: EphemerisKernel, now: datetime
) -> None:
    """Test that full moons beyond the phase search window yield no values."""
    with (
        patch("custom_components.moon_astro.astronomy.PHASE_SEARCH_DAYS_BACK", 1.0),
        patch("custom_components.moon_astro.astronomy.PHASE_SEARCH_DAYS_AHEAD", 1.0),
    ):
        payload = _events(ephemeris, now)

    assert [
        payload[key]
        for key in (
            KEY_PREVIOUS_FULL_MOON,
            KEY_PREVIOUS_FULL_MOON_NAME,
            KEY_PREVIOUS_FULL_MOON_ALT_NAMES,
            KEY_ZODIAC_SIGN_PREVIOUS_FULL_MOON,
            KEY_NEXT_FULL_MOON,
            KEY_NEXT_FULL_MOON_NAME,
            KEY_NEXT_FULL_MOON_ALT_NAMES,
            KEY_ZODIAC_SIGN_NEXT_FULL_MOON,
        )
    ] == [None] * 8


def test_phase_angles_at_phase_instants(
    ephemeris: EphemerisKernel, now: datetime
) -> None:
    """Test that the phase angle reaches its target within the published minute."""
    eph, ts = ephemeris
    payload = _events(ephemeris, now)

    for key, target in (
        (KEY_PREVIOUS_NEW_MOON, 0.0),
        (KEY_NEXT_NEW_MOON, 0.0),
        (KEY_PREVIOUS_FIRST_QUARTER, 90.0),
        (KEY_NEXT_FIRST_QUARTER, 90.0),
        (KEY_PREVIOUS_FULL_MOON, 180.0),
        (KEY_NEXT_FULL_MOON, 180.0),
        (KEY_PREVIOUS_LAST_QUARTER, 270.0),
        (KEY_NEXT_LAST_QUARTER, 270.0),
    ):
        angle = almanac.moon_phase(eph, ts.from_datetime(payload[key])).degrees
        # The phase angle grows by 0.0085 degree per minute.
        assert abs((angle - target + 180.0) % 360.0 - 180.0) < 0.005, key


def test_apsides_are_distance_extremes(
    ephemeris: EphemerisKernel, now: datetime
) -> None:
    """Test that the apsides are the distance extremes within three days."""
    eph, ts = ephemeris
    payload = _events(ephemeris, now)
    earth_to_moon = eph["moon"] - eph["earth"]
    offsets = [timedelta(minutes=minutes) for minutes in range(-4320, 4321)]

    for key, extreme in (
        (KEY_PREVIOUS_PERIGEE, np.argmin),
        (KEY_NEXT_PERIGEE, np.argmin),
        (KEY_PREVIOUS_APOGEE, np.argmax),
        (KEY_NEXT_APOGEE, np.argmax),
    ):
        around = ts.from_datetimes([payload[key] + offset for offset in offsets])
        distances = earth_to_moon.at(around).distance().km
        assert offsets[int(extreme(distances))] == timedelta(0), key


def test_blue_moon(ephemeris: EphemerisKernel) -> None:
    """Test that the second full moon of May 2026 is named a blue moon."""
    payload = _events(ephemeris, _utc("2026-06-10T00:00"))

    assert payload[KEY_PREVIOUS_FULL_MOON] == _utc("2026-05-31T08:45")
    assert payload[KEY_PREVIOUS_FULL_MOON_NAME] == "blue_moon"
    assert payload[KEY_PREVIOUS_FULL_MOON_ALT_NAMES] == ""
    assert payload[KEY_NEXT_FULL_MOON_NAME] == "strawberry_moon"
    assert payload[KEY_NEXT_FULL_MOON_ALT_NAMES] == "strawberry_moon_alt_names"


@pytest.mark.parametrize(
    ("tz", "expected"),
    [
        pytest.param(UTC, ["worm_moon", "blue_moon"], id="utc"),
        pytest.param(ZoneInfo("Europe/Paris"), ["worm_moon", "pink_moon"], id="paris"),
    ],
)
def test_full_moon_names_follow_time_zone(
    ephemeris: EphemerisKernel, tz: tzinfo, expected: list[str]
) -> None:
    """Test that the calendar month naming a full moon is that of the time zone.

    A full moon at 23:30 UTC on March 31 falls on April 1 in Paris.
    """
    _eph, ts = ephemeris
    full_moons = [ts.utc(2026, 3, 1, 12), ts.utc(2026, 3, 31, 23, 30)]

    assert list(_full_moon_name_codes(full_moons, tz).values()) == expected
