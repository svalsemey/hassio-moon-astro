"""Tests for the Moon Astro sensors."""

from datetime import timedelta
from unittest.mock import patch

from custom_components.moon_astro.const import DEFAULT_EVENTS_REFRESH_FALLBACK
from freezegun.api import FrozenDateTimeFactory
import pytest
from pytest_homeassistant_custom_component.common import (
    MockConfigEntry,
    async_fire_time_changed,
    snapshot_platform,
)
from syrupy.assertion import SnapshotAssertion

from homeassistant.const import STATE_UNAVAILABLE, Platform
from homeassistant.core import HomeAssistant
from homeassistant.helpers import entity_registry as er

# Beyond the periodic interval of both coordinators, including the random fraction of
# a second that Home Assistant adds when it schedules a refresh.
REFRESH_DELAY = timedelta(seconds=DEFAULT_EVENTS_REFRESH_FALLBACK, minutes=1)


@pytest.fixture
def platforms() -> list[Platform]:
    """Return the platforms to set up."""
    return [Platform.SENSOR]


@pytest.mark.usefixtures("entity_registry_enabled_by_default")
async def test_sensors(
    hass: HomeAssistant,
    entity_registry: er.EntityRegistry,
    snapshot: SnapshotAssertion,
    init_integration: MockConfigEntry,
) -> None:
    """Test the sensor entities and their states."""
    await snapshot_platform(hass, entity_registry, snapshot, init_integration.entry_id)


@pytest.mark.usefixtures("init_integration")
@pytest.mark.parametrize(
    "entity_id",
    [
        "sensor.moon_astro_ecliptic_longitude_topocentric",
        "sensor.moon_astro_ecliptic_latitude_topocentric",
        "sensor.moon_astro_ecliptic_longitude_at_next_new_moon",
        "sensor.moon_astro_ecliptic_latitude_at_next_new_moon",
        "sensor.moon_astro_ecliptic_longitude_at_next_full_moon",
        "sensor.moon_astro_ecliptic_latitude_at_next_full_moon",
        "sensor.moon_astro_ecliptic_longitude_at_previous_new_moon",
        "sensor.moon_astro_ecliptic_latitude_at_previous_new_moon",
        "sensor.moon_astro_ecliptic_longitude_at_previous_full_moon",
        "sensor.moon_astro_ecliptic_latitude_at_previous_full_moon",
        "sensor.moon_astro_zodiac_degree_at_next_new_moon",
        "sensor.moon_astro_zodiac_degree_at_next_full_moon",
        "sensor.moon_astro_zodiac_degree_at_previous_new_moon",
        "sensor.moon_astro_zodiac_degree_at_previous_full_moon",
    ],
)
async def test_detailed_sensors_disabled_by_default(
    hass: HomeAssistant, entity_registry: er.EntityRegistry, entity_id: str
) -> None:
    """Test that the detailed coordinates are registered disabled."""
    assert (entry := entity_registry.async_get(entity_id))
    assert entry.disabled_by is er.RegistryEntryDisabler.INTEGRATION
    assert hass.states.get(entity_id) is None


@pytest.mark.usefixtures("init_integration")
@pytest.mark.parametrize(
    ("computation", "entity_id"),
    [
        pytest.param(
            "compute_current_payload", "sensor.moon_astro_azimuth", id="current"
        ),
        pytest.param(
            "compute_events_payload", "sensor.moon_astro_next_new_moon", id="events"
        ),
    ],
)
async def test_unavailable_while_computation_fails(
    hass: HomeAssistant,
    freezer: FrozenDateTimeFactory,
    computation: str,
    entity_id: str,
) -> None:
    """Test that sensors are unavailable while their computation fails."""
    state = hass.states.get(entity_id)
    assert state.state != STATE_UNAVAILABLE

    with patch(
        f"custom_components.moon_astro.coordinator.{computation}",
        side_effect=ArithmeticError("diverged"),
    ):
        freezer.tick(REFRESH_DELAY)
        async_fire_time_changed(hass)
        await hass.async_block_till_done(wait_background_tasks=True)

    assert hass.states.get(entity_id).state == STATE_UNAVAILABLE

    freezer.tick(REFRESH_DELAY)
    async_fire_time_changed(hass)
    await hass.async_block_till_done(wait_background_tasks=True)

    assert hass.states.get(entity_id).state != STATE_UNAVAILABLE
