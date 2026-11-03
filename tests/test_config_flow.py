"""Tests for the Moon Astro config and options flows."""

from typing import Any
from unittest.mock import AsyncMock

from custom_components.moon_astro.const import (
    CONF_ALT,
    CONF_EVENTS_REFRESH_FALLBACK,
    CONF_LAT,
    CONF_LON,
    CONF_SCAN_INTERVAL,
    CONF_TIME_ZONE,
    CONF_USE_HA_TZ,
    DOMAIN,
    NAME,
)
from custom_components.moon_astro.utils import EphemerisError, EphemerisKernel
import pytest
from pytest_homeassistant_custom_component.common import (
    MockConfigEntry,
    get_schema_suggested_value,
)

from homeassistant.config_entries import SOURCE_USER, ConfigEntryState, ConfigFlowResult
from homeassistant.core import HomeAssistant
from homeassistant.data_entry_flow import FlowResultType, InvalidData

OPTIONS = {
    CONF_SCAN_INTERVAL: 120,
    CONF_USE_HA_TZ: False,
    CONF_TIME_ZONE: "Asia/Tokyo",
    CONF_EVENTS_REFRESH_FALLBACK: 7200,
}


async def _async_submit_location(
    hass: HomeAssistant, flow_id: str, location: dict[str, Any]
) -> ConfigFlowResult:
    """Submit the location form and let the ephemeris preparation complete."""
    result = await hass.config_entries.flow.async_configure(flow_id, location)
    assert result["type"] is FlowResultType.SHOW_PROGRESS
    assert result["progress_action"] == "download_ephemeris"
    await hass.async_block_till_done()
    return await hass.config_entries.flow.async_configure(flow_id)


@pytest.mark.usefixtures("mock_get_ephemeris")
async def test_user_flow(hass: HomeAssistant, location: dict[str, float]) -> None:
    """Test the creation of the entry from the Home Assistant location."""
    result = await hass.config_entries.flow.async_init(
        DOMAIN, context={"source": SOURCE_USER}
    )
    assert result["type"] is FlowResultType.FORM
    assert result["step_id"] == "user"
    schema = result["data_schema"].schema
    assert get_schema_suggested_value(schema, CONF_LAT) == hass.config.latitude
    assert get_schema_suggested_value(schema, CONF_LON) == hass.config.longitude
    assert get_schema_suggested_value(schema, CONF_ALT) == hass.config.elevation

    result = await _async_submit_location(hass, result["flow_id"], location)

    assert result["type"] is FlowResultType.CREATE_ENTRY
    assert result["title"] == NAME
    assert result["data"] == location
    assert result["result"].unique_id == DOMAIN
    assert result["result"].state is ConfigEntryState.LOADED


async def test_user_flow_recovers_from_download_failure(
    hass: HomeAssistant,
    mock_get_ephemeris: AsyncMock,
    ephemeris: EphemerisKernel,
    location: dict[str, float],
) -> None:
    """Test that a failed download shows the location form again with an error."""
    mock_get_ephemeris.side_effect = [EphemerisError("unreachable"), ephemeris]
    result = await hass.config_entries.flow.async_init(
        DOMAIN, context={"source": SOURCE_USER}
    )

    result = await _async_submit_location(hass, result["flow_id"], location)

    assert result["type"] is FlowResultType.FORM
    assert result["step_id"] == "user"
    assert result["errors"] == {"base": "download_failed"}
    assert (
        get_schema_suggested_value(result["data_schema"].schema, CONF_LAT)
        == (location[CONF_LAT])
    )

    result = await _async_submit_location(hass, result["flow_id"], location)

    assert result["type"] is FlowResultType.CREATE_ENTRY
    assert result["data"] == location


@pytest.mark.parametrize(
    "invalid",
    [
        pytest.param({CONF_LAT: 90.5}, id="latitude"),
        pytest.param({CONF_LON: -180.5}, id="longitude"),
        pytest.param({CONF_ALT: 10001.0}, id="elevation_high"),
        pytest.param({CONF_ALT: -501.0}, id="elevation_low"),
    ],
)
async def test_user_flow_rejects_invalid_location(
    hass: HomeAssistant, location: dict[str, float], invalid: dict[str, float]
) -> None:
    """Test that out of range coordinates are rejected."""
    result = await hass.config_entries.flow.async_init(
        DOMAIN, context={"source": SOURCE_USER}
    )

    with pytest.raises(InvalidData):
        await hass.config_entries.flow.async_configure(
            result["flow_id"], location | invalid
        )


async def test_user_flow_single_instance(
    hass: HomeAssistant, mock_config_entry: MockConfigEntry
) -> None:
    """Test that a second entry cannot be created."""
    mock_config_entry.add_to_hass(hass)

    result = await hass.config_entries.flow.async_init(
        DOMAIN, context={"source": SOURCE_USER}
    )

    assert result["type"] is FlowResultType.ABORT
    assert result["reason"] == "single_instance_allowed"


async def test_reconfigure_flow(
    hass: HomeAssistant,
    init_integration: MockConfigEntry,
    location: dict[str, float],
) -> None:
    """Test that the reconfiguration updates the location and reloads the entry."""
    result = await init_integration.start_reconfigure_flow(hass)
    assert result["type"] is FlowResultType.FORM
    assert result["step_id"] == "reconfigure"
    assert (
        get_schema_suggested_value(result["data_schema"].schema, CONF_LAT)
        == (location[CONF_LAT])
    )
    svalbard = {CONF_LAT: 78.2232, CONF_LON: 15.6267, CONF_ALT: 10.0}

    result = await hass.config_entries.flow.async_configure(result["flow_id"], svalbard)
    await hass.async_block_till_done()

    assert result["type"] is FlowResultType.ABORT
    assert result["reason"] == "reconfigure_successful"
    assert init_integration.data == svalbard
    assert init_integration.state is ConfigEntryState.LOADED


async def test_options_flow(
    hass: HomeAssistant, init_integration: MockConfigEntry
) -> None:
    """Test that the options are stored."""
    result = await hass.config_entries.options.async_init(init_integration.entry_id)
    assert result["type"] is FlowResultType.FORM
    assert result["step_id"] == "init"

    result = await hass.config_entries.options.async_configure(
        result["flow_id"], OPTIONS
    )
    await hass.async_block_till_done()

    assert result["type"] is FlowResultType.CREATE_ENTRY
    assert init_integration.options == OPTIONS


async def test_options_flow_suggests_stored_options(
    hass: HomeAssistant,
    mock_config_entry: MockConfigEntry,
    mock_get_ephemeris: AsyncMock,
) -> None:
    """Test that the stored options pre-fill the options form."""
    mock_config_entry.add_to_hass(hass)
    hass.config_entries.async_update_entry(mock_config_entry, options=OPTIONS)
    assert await hass.config_entries.async_setup(mock_config_entry.entry_id)
    await hass.async_block_till_done()

    result = await hass.config_entries.options.async_init(mock_config_entry.entry_id)

    schema = result["data_schema"].schema
    assert {key: get_schema_suggested_value(schema, key) for key in OPTIONS} == OPTIONS


@pytest.mark.parametrize(
    "options",
    [
        pytest.param({**OPTIONS, CONF_SCAN_INTERVAL: 59}, id="scan_interval_low"),
        pytest.param({**OPTIONS, CONF_SCAN_INTERVAL: 21601}, id="scan_interval_high"),
        pytest.param({**OPTIONS, CONF_TIME_ZONE: "Mars/Olympus_Mons"}, id="time_zone"),
        pytest.param(
            {**OPTIONS, CONF_EVENTS_REFRESH_FALLBACK: 3599}, id="fallback_low"
        ),
        pytest.param(
            {**OPTIONS, CONF_EVENTS_REFRESH_FALLBACK: 604801}, id="fallback_high"
        ),
    ],
)
async def test_options_flow_rejects_invalid_options(
    hass: HomeAssistant, init_integration: MockConfigEntry, options: dict[str, Any]
) -> None:
    """Test that out of range options are rejected."""
    result = await hass.config_entries.options.async_init(init_integration.entry_id)

    with pytest.raises(InvalidData):
        await hass.config_entries.options.async_configure(result["flow_id"], options)
