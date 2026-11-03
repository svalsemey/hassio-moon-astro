"""Tests for the Moon Astro setup."""

from datetime import timedelta
from unittest.mock import AsyncMock, patch

from custom_components.moon_astro.const import (
    CONF_EVENTS_REFRESH_FALLBACK,
    CONF_SCAN_INTERVAL,
    CONF_USE_HA_TZ,
)
from custom_components.moon_astro.utils import EphemerisError
import pytest
from pytest_homeassistant_custom_component.common import MockConfigEntry

from homeassistant.config_entries import ConfigEntryState
from homeassistant.core import HomeAssistant


async def test_setup_and_unload(
    hass: HomeAssistant, init_integration: MockConfigEntry
) -> None:
    """Test that the entry loads and unloads."""
    assert init_integration.state is ConfigEntryState.LOADED

    assert await hass.config_entries.async_unload(init_integration.entry_id)
    assert init_integration.state is ConfigEntryState.NOT_LOADED


async def test_setup_retry_without_ephemeris(
    hass: HomeAssistant,
    mock_config_entry: MockConfigEntry,
    mock_get_ephemeris: AsyncMock,
) -> None:
    """Test that the setup is retried when the ephemeris is not available."""
    mock_get_ephemeris.side_effect = EphemerisError("download failed")
    mock_config_entry.add_to_hass(hass)

    await hass.config_entries.async_setup(mock_config_entry.entry_id)

    assert mock_config_entry.state is ConfigEntryState.SETUP_RETRY
    assert mock_config_entry.error_reason_translation_key == "ephemeris_unavailable"


@pytest.mark.usefixtures("mock_get_ephemeris")
@pytest.mark.parametrize(
    "computation", ["compute_current_payload", "compute_events_payload"]
)
async def test_setup_retry_on_computation_failure(
    hass: HomeAssistant, mock_config_entry: MockConfigEntry, computation: str
) -> None:
    """Test that the setup is retried when a first computation fails."""
    mock_config_entry.add_to_hass(hass)

    with patch(
        f"custom_components.moon_astro.coordinator.{computation}",
        side_effect=ValueError("outside the ephemeris range"),
    ):
        await hass.config_entries.async_setup(mock_config_entry.entry_id)

    assert mock_config_entry.state is ConfigEntryState.SETUP_RETRY


async def test_options_update_reloads_entry(
    hass: HomeAssistant, init_integration: MockConfigEntry
) -> None:
    """Test that new options are applied by a reload."""
    hass.config_entries.async_update_entry(
        init_integration,
        options={
            CONF_SCAN_INTERVAL: 120,
            CONF_USE_HA_TZ: True,
            CONF_EVENTS_REFRESH_FALLBACK: 7200,
        },
    )
    await hass.async_block_till_done()

    assert init_integration.state is ConfigEntryState.LOADED
    runtime = init_integration.runtime_data
    assert runtime.coordinator.update_interval == timedelta(seconds=120)
    assert runtime.events_coordinator.update_interval == timedelta(seconds=7200)


async def test_remove_entry_discards_ephemeris(
    hass: HomeAssistant, init_integration: MockConfigEntry
) -> None:
    """Test that removing the entry deletes the cached ephemeris."""
    with patch("custom_components.moon_astro.async_discard_ephemeris") as mock_discard:
        assert await hass.config_entries.async_remove(init_integration.entry_id)

    mock_discard.assert_awaited_once_with(hass)
