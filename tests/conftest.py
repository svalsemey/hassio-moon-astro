"""Fixtures for the Moon Astro tests."""

from collections.abc import Generator
from datetime import UTC, datetime
from pathlib import Path
from unittest.mock import AsyncMock, patch

from custom_components.moon_astro.const import (
    CONF_ALT,
    CONF_LAT,
    CONF_LON,
    DOMAIN,
    NAME,
)
from custom_components.moon_astro.utils import EphemerisKernel
from freezegun.api import FrozenDateTimeFactory
import pytest
from pytest_homeassistant_custom_component.common import MockConfigEntry
from pytest_homeassistant_custom_component.syrupy import HomeAssistantSnapshotExtension
from skyfield.api import load
from skyfield.jpllib import SpiceKernel
from syrupy.assertion import SnapshotAssertion

from homeassistant.const import Platform
from homeassistant.core import HomeAssistant


@pytest.fixture(autouse=True)
def auto_enable_custom_integrations(enable_custom_integrations: None) -> None:
    """Enable the loading of custom integrations."""


@pytest.fixture(autouse=True)
def now(freezer: FrozenDateTimeFactory) -> datetime:
    """Freeze the time at an instant whose search windows lie within the excerpt."""
    instant = datetime(2026, 11, 2, 22, 52, tzinfo=UTC)
    freezer.move_to(instant)
    return instant


@pytest.fixture
def snapshot(snapshot: SnapshotAssertion) -> SnapshotAssertion:
    """Serialize Home Assistant objects as the core tests do."""
    return snapshot.use_extension(HomeAssistantSnapshotExtension)


@pytest.fixture
def entity_registry_enabled_by_default() -> Generator[None]:
    """Register every entity enabled, including those disabled by default."""
    with patch(
        "homeassistant.helpers.entity.Entity.entity_registry_enabled_default",
        return_value=True,
    ):
        yield


@pytest.fixture
def ephemeris_path() -> Path:
    """Return the excerpt of DE440 holding its coefficients for 2026 and 2027.

    Results match those of the full kernel. The file is produced from it with:
        python -m jplephem excerpt 2026/1/1 2028/1/1 de440.bsp de440_2026_2027.bsp
    """
    return Path(__file__).parent / "fixtures" / "de440_2026_2027.bsp"


@pytest.fixture
def ephemeris(ephemeris_path: Path) -> Generator[EphemerisKernel]:
    """Return the test ephemeris kernel and a timescale."""
    kernel = SpiceKernel(str(ephemeris_path))
    yield kernel, load.timescale()
    kernel.close()


@pytest.fixture
def mock_get_ephemeris(ephemeris: EphemerisKernel) -> Generator[AsyncMock]:
    """Serve the test ephemeris to the setup and the config flow."""
    with (
        patch(
            "custom_components.moon_astro.async_get_ephemeris", return_value=ephemeris
        ) as mock,
        patch("custom_components.moon_astro.config_flow.async_get_ephemeris", mock),
    ):
        yield mock


@pytest.fixture
def location() -> dict[str, float]:
    """Return the location of the Paris Observatory."""
    return {CONF_LAT: 48.8363, CONF_LON: 2.3366, CONF_ALT: 67.0}


@pytest.fixture
def mock_config_entry(location: dict[str, float]) -> MockConfigEntry:
    """Return the Moon Astro config entry."""
    return MockConfigEntry(
        domain=DOMAIN,
        title=NAME,
        unique_id=DOMAIN,
        data=location,
        # Fixed, as the unique ids of the entities derive from it.
        entry_id="01JZ6V8T4S3AN2W8MBY5HK6X7Q",
    )


@pytest.fixture
def platforms() -> list[Platform]:
    """Return the platforms to set up."""
    return [Platform.BINARY_SENSOR, Platform.SENSOR]


@pytest.fixture
async def init_integration(
    hass: HomeAssistant,
    mock_config_entry: MockConfigEntry,
    mock_get_ephemeris: AsyncMock,
    platforms: list[Platform],
) -> MockConfigEntry:
    """Set up the integration with the test ephemeris."""
    mock_config_entry.add_to_hass(hass)
    with patch("custom_components.moon_astro.PLATFORMS", platforms):
        assert await hass.config_entries.async_setup(mock_config_entry.entry_id)
        await hass.async_block_till_done()
    return mock_config_entry
