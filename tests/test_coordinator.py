"""Tests for the refresh scheduling of the Moon Astro coordinators."""

from datetime import datetime, timedelta
from typing import Any
from unittest.mock import patch

from custom_components.moon_astro.astronomy import compute_events_payload
from custom_components.moon_astro.const import (
    CONF_EVENTS_REFRESH_FALLBACK,
    DEFAULT_EVENTS_REFRESH_FALLBACK,
    KEY_NEXT_APOGEE,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_LAST_QUARTER,
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_PERIGEE,
    KEY_NEXT_RISE,
    KEY_NEXT_SET,
)
from custom_components.moon_astro.coordinator import ComputationStats
from freezegun.api import FrozenDateTimeFactory
import pytest
from pytest_homeassistant_custom_component.common import (
    MockConfigEntry,
    async_fire_time_changed,
)

from homeassistant.core import HomeAssistant
from homeassistant.util import dt as dt_util

NEXT_RISE = "sensor.moon_astro_next_rise"
PREVIOUS_RISE = "sensor.moon_astro_previous_rise"
NEXT_EVENT_ENTITIES = (
    NEXT_RISE,
    "sensor.moon_astro_next_set",
    "sensor.moon_astro_next_new_moon",
    "sensor.moon_astro_next_first_quarter",
    "sensor.moon_astro_next_full_moon",
    "sensor.moon_astro_next_last_quarter",
    "sensor.moon_astro_next_apogee",
    "sensor.moon_astro_next_perigee",
)
NEXT_EVENT_KEYS = (
    KEY_NEXT_RISE,
    KEY_NEXT_SET,
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_LAST_QUARTER,
    KEY_NEXT_APOGEE,
    KEY_NEXT_PERIGEE,
)


def _instant(hass: HomeAssistant, entity_id: str) -> datetime:
    """Return the instant held by a timestamp sensor."""
    return datetime.fromisoformat(hass.states.get(entity_id).state)


def _next_update(hass: HomeAssistant) -> datetime:
    """Return the next recomputation instant of the event-based sensors."""
    return hass.states.get(NEXT_RISE).attributes["next_update"]


@pytest.mark.usefixtures("init_integration")
async def test_refresh_after_next_event(
    hass: HomeAssistant, freezer: FrozenDateTimeFactory
) -> None:
    """Test that the event-based sensors are recomputed after the next event."""
    next_rise = _instant(hass, NEXT_RISE)
    earliest = min(_instant(hass, entity_id) for entity_id in NEXT_EVENT_ENTITIES)
    # The moonrise is the earliest upcoming event at the test location and time.
    assert earliest == next_rise
    assert _next_update(hass) == next_rise + timedelta(minutes=2)

    freezer.move_to(_next_update(hass))
    async_fire_time_changed(hass)
    await hass.async_block_till_done(wait_background_tasks=True)

    assert _instant(hass, PREVIOUS_RISE) == next_rise
    assert _instant(hass, NEXT_RISE) > next_rise
    assert _next_update(hass) > dt_util.utcnow()


@pytest.mark.usefixtures("mock_get_ephemeris")
async def test_fallback_refresh_before_next_event(
    hass: HomeAssistant,
    freezer: FrozenDateTimeFactory,
    mock_config_entry: MockConfigEntry,
    now: datetime,
) -> None:
    """Test that the fallback interval forces a recomputation before the next event."""
    mock_config_entry.add_to_hass(hass)
    hass.config_entries.async_update_entry(
        mock_config_entry, options={CONF_EVENTS_REFRESH_FALLBACK: 3600}
    )
    with patch(
        "custom_components.moon_astro.coordinator.compute_events_payload",
        wraps=compute_events_payload,
    ) as mock_compute:
        assert await hass.config_entries.async_setup(mock_config_entry.entry_id)
        await hass.async_block_till_done()
        assert _next_update(hass) == now + timedelta(hours=1)

        freezer.tick(timedelta(hours=1))
        async_fire_time_changed(hass)
        await hass.async_block_till_done(wait_background_tasks=True)

    assert mock_compute.call_count == 2


@pytest.mark.usefixtures("mock_get_ephemeris")
async def test_fallback_refresh_without_upcoming_event(
    hass: HomeAssistant, mock_config_entry: MockConfigEntry, now: datetime
) -> None:
    """Test that the fallback interval drives the refresh when no event is found."""

    def _without_next_events(*args: Any) -> dict[str, Any]:
        return compute_events_payload(*args) | dict.fromkeys(NEXT_EVENT_KEYS)

    mock_config_entry.add_to_hass(hass)
    with patch(
        "custom_components.moon_astro.coordinator.compute_events_payload",
        side_effect=_without_next_events,
    ):
        assert await hass.config_entries.async_setup(mock_config_entry.entry_id)
        await hass.async_block_till_done()

    assert _next_update(hass) == now + timedelta(
        seconds=DEFAULT_EVENTS_REFRESH_FALLBACK
    )


@pytest.mark.usefixtures("mock_get_ephemeris")
async def test_no_refresh_after_unload(
    hass: HomeAssistant,
    freezer: FrozenDateTimeFactory,
    mock_config_entry: MockConfigEntry,
) -> None:
    """Test that unloading the entry cancels the scheduled recomputation."""
    mock_config_entry.add_to_hass(hass)
    with patch(
        "custom_components.moon_astro.coordinator.compute_events_payload",
        wraps=compute_events_payload,
    ) as mock_compute:
        assert await hass.config_entries.async_setup(mock_config_entry.entry_id)
        await hass.async_block_till_done()
        next_update = _next_update(hass)
        assert await hass.config_entries.async_unload(mock_config_entry.entry_id)

        freezer.move_to(next_update)
        async_fire_time_changed(hass)
        await hass.async_block_till_done(wait_background_tasks=True)

    assert mock_compute.call_count == 1


def test_computation_stats(now: datetime) -> None:
    """Test the aggregation of the computation durations."""
    stats = ComputationStats()

    for minutes, duration in ((0, 0.004), (5, 0.002), (10, 0.006)):
        stats.record(now + timedelta(minutes=minutes), duration)

    assert stats.count == 3
    assert stats.last == 0.006
    assert stats.minimum == 0.002
    assert stats.maximum == 0.006
    assert stats.total == pytest.approx(0.012)
    assert stats.reference_time == now + timedelta(minutes=10)
