"""Tests for the Moon Astro diagnostics."""

from collections.abc import Generator
from datetime import timedelta
from unittest.mock import MagicMock, patch

from custom_components.moon_astro.const import (
    CONF_ALT,
    CONF_LAT,
    CONF_LON,
    DEFAULT_SCAN_INTERVAL,
    KEY_AZIMUTH,
    KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
    KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
    KEY_ELEVATION,
)
from freezegun.api import FrozenDateTimeFactory
import pytest
from pytest_homeassistant_custom_component.common import (
    MockConfigEntry,
    async_fire_time_changed,
)
from pytest_homeassistant_custom_component.components.diagnostics import (
    get_diagnostics_for_config_entry,
)
from pytest_homeassistant_custom_component.typing import ClientSessionGenerator
from syrupy.assertion import SnapshotAssertion

from homeassistant.core import HomeAssistant


@pytest.fixture
def mock_perf_counter() -> Generator[MagicMock]:
    """Make the computation durations deterministic.

    The setup computes the current position (4 ms) then the events (50 ms); the next
    periodic refresh computes the current position again (2 ms).
    """
    with patch("custom_components.moon_astro.coordinator.time") as mock_time:
        mock_time.perf_counter.side_effect = [1.0, 1.004, 2.0, 2.05, 3.0, 3.002]
        yield mock_time.perf_counter


@pytest.mark.usefixtures("mock_perf_counter")
async def test_diagnostics(
    hass: HomeAssistant,
    hass_client: ClientSessionGenerator,
    freezer: FrozenDateTimeFactory,
    snapshot: SnapshotAssertion,
    init_integration: MockConfigEntry,
) -> None:
    """Test the config entry diagnostics."""
    freezer.tick(timedelta(seconds=DEFAULT_SCAN_INTERVAL))
    async_fire_time_changed(hass)
    await hass.async_block_till_done(wait_background_tasks=True)

    diagnostics = await get_diagnostics_for_config_entry(
        hass, hass_client, init_integration
    )

    assert diagnostics == snapshot
    # The location must not identify the home: about 10 km and 100 m only.
    assert diagnostics["entry"]["data"] == {
        CONF_LAT: 48.8,
        CONF_LON: 2.3,
        CONF_ALT: 100.0,
    }
    current = diagnostics["coordinators"]["current"]["data"]
    for key in (
        KEY_AZIMUTH,
        KEY_ELEVATION,
        KEY_ECLIPTIC_LONGITUDE_TOPOCENTRIC,
        KEY_ECLIPTIC_LATITUDE_TOPOCENTRIC,
    ):
        assert current[key] == round(current[key], 1)
