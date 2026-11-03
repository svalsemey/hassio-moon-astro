"""Tests for the Moon Astro binary sensor."""

from custom_components.moon_astro.const import CONF_ALT, CONF_LAT, CONF_LON, DOMAIN
import pytest
from pytest_homeassistant_custom_component.common import (
    MockConfigEntry,
    snapshot_platform,
)
from syrupy.assertion import SnapshotAssertion

from homeassistant.const import STATE_OFF, STATE_ON, Platform
from homeassistant.core import HomeAssistant
from homeassistant.helpers import entity_registry as er

ENTITY_ID = "binary_sensor.moon_astro_above_horizon"


@pytest.fixture
def platforms() -> list[Platform]:
    """Return the platforms to set up."""
    return [Platform.BINARY_SENSOR]


async def test_binary_sensor(
    hass: HomeAssistant,
    entity_registry: er.EntityRegistry,
    snapshot: SnapshotAssertion,
    init_integration: MockConfigEntry,
) -> None:
    """Test the binary sensor entity and its state."""
    await snapshot_platform(hass, entity_registry, snapshot, init_integration.entry_id)


@pytest.mark.usefixtures("mock_get_ephemeris")
@pytest.mark.parametrize(
    ("latitude", "longitude", "expected"),
    [
        pytest.param(-33.8688, 151.2093, STATE_ON, id="sydney"),
        pytest.param(40.7128, -74.006, STATE_OFF, id="new_york"),
    ],
)
async def test_above_horizon_depends_on_location(
    hass: HomeAssistant, latitude: float, longitude: float, expected: str
) -> None:
    """Test that the state follows the Moon altitude at the observer location."""
    entry = MockConfigEntry(
        domain=DOMAIN,
        unique_id=DOMAIN,
        data={CONF_LAT: latitude, CONF_LON: longitude, CONF_ALT: 0.0},
    )
    entry.add_to_hass(hass)

    assert await hass.config_entries.async_setup(entry.entry_id)
    await hass.async_block_till_done()

    assert hass.states.get(ENTITY_ID).state == expected
