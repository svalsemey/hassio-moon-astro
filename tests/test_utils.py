"""Tests for the Moon Astro ephemeris lifecycle and configuration helpers."""

import asyncio
from collections.abc import Generator
from pathlib import Path
import shutil
from typing import Any
from unittest.mock import patch

import aiohttp
from custom_components.moon_astro.const import (
    CACHE_DIR_NAME,
    CONF_TIME_ZONE,
    CONF_USE_HA_TZ,
    DE440_FILE,
    DE440_URL,
    DOMAIN,
)
from custom_components.moon_astro.utils import (
    EphemerisError,
    async_discard_ephemeris,
    async_get_ephemeris,
    async_resolve_time_zone,
    available_time_zones,
)
import pytest
from pytest_homeassistant_custom_component.common import MockConfigEntry
from pytest_homeassistant_custom_component.test_util.aiohttp import AiohttpClientMocker

from homeassistant.core import HomeAssistant

# Below the size of the excerpt, which then passes for a complete kernel.
MIN_SIZE = 100_000

DAMAGED_FILES = [
    pytest.param(b"\0" * 10, id="truncated"),
    pytest.param(b"\0" * MIN_SIZE, id="unreadable"),
]


@pytest.fixture
def hass_config_dir(hass_tmp_config_dir: str) -> str:
    """Use a writable configuration directory."""
    return hass_tmp_config_dir


@pytest.fixture(autouse=True)
def accept_excerpt() -> Generator[None]:
    """Accept the excerpt as a complete kernel."""
    with patch("custom_components.moon_astro.utils.MIN_EPHEMERIS_SIZE_BYTES", MIN_SIZE):
        yield


@pytest.fixture
def cache_dir(hass_config_dir: str) -> Path:
    """Return the directory holding the cached kernel."""
    return Path(hass_config_dir) / CACHE_DIR_NAME


@pytest.fixture
def kernel_path(cache_dir: Path) -> Path:
    """Return the path of the cached kernel."""
    return cache_dir / DE440_FILE


@pytest.fixture
def cached_kernel(kernel_path: Path, ephemeris_path: Path) -> Path:
    """Store the excerpt as the cached kernel."""
    kernel_path.parent.mkdir()
    shutil.copyfile(ephemeris_path, kernel_path)
    return kernel_path


@pytest.fixture
def jpl_server(
    aioclient_mock: AiohttpClientMocker, ephemeris_path: Path
) -> AiohttpClientMocker:
    """Serve the excerpt from the JPL URL."""
    aioclient_mock.get(DE440_URL, content=ephemeris_path.read_bytes())
    return aioclient_mock


def _files(directory: Path) -> list[str]:
    """Return the sorted names of the files of a directory."""
    return sorted(path.name for path in directory.iterdir())


@pytest.mark.usefixtures("cached_kernel")
async def test_cached_kernel_loaded_once(
    hass: HomeAssistant, aioclient_mock: AiohttpClientMocker
) -> None:
    """Test that a valid cached kernel is loaded once and not downloaded."""
    first = await async_get_ephemeris(hass)

    assert await async_get_ephemeris(hass) is first
    assert aioclient_mock.call_count == 0


async def test_missing_kernel_downloaded(
    hass: HomeAssistant, jpl_server: AiohttpClientMocker, cache_dir: Path
) -> None:
    """Test that a missing kernel is downloaded once, even for concurrent callers."""
    first, second = await asyncio.gather(
        async_get_ephemeris(hass), async_get_ephemeris(hass)
    )

    assert first is second
    assert jpl_server.call_count == 1
    assert await hass.async_add_executor_job(_files, cache_dir) == [DE440_FILE]


@pytest.mark.parametrize("content", DAMAGED_FILES)
async def test_damaged_kernel_downloaded_again(
    hass: HomeAssistant,
    jpl_server: AiohttpClientMocker,
    kernel_path: Path,
    ephemeris_path: Path,
    caplog: pytest.LogCaptureFixture,
    content: bytes,
) -> None:
    """Test that a damaged cached kernel is replaced by a fresh copy."""
    kernel_path.parent.mkdir()
    kernel_path.write_bytes(content)

    await async_get_ephemeris(hass)

    assert "Discarding the ephemeris file" in caplog.text
    assert jpl_server.call_count == 1
    assert kernel_path.stat().st_size == ephemeris_path.stat().st_size


@pytest.mark.parametrize(
    "response",
    [
        pytest.param({"status": 503}, id="http_error"),
        pytest.param({"exc": aiohttp.ClientError()}, id="client_error"),
        pytest.param({"exc": TimeoutError()}, id="timeout"),
    ],
)
async def test_download_failure(
    hass: HomeAssistant,
    aioclient_mock: AiohttpClientMocker,
    cache_dir: Path,
    response: dict[str, Any],
) -> None:
    """Test that a failed download raises and leaves no file behind."""
    aioclient_mock.get(DE440_URL, **response)

    with pytest.raises(EphemerisError, match="Download of"):
        await async_get_ephemeris(hass)

    assert await hass.async_add_executor_job(_files, cache_dir) == []


@pytest.mark.parametrize("content", DAMAGED_FILES)
async def test_downloaded_kernel_rejected(
    hass: HomeAssistant,
    aioclient_mock: AiohttpClientMocker,
    cache_dir: Path,
    content: bytes,
) -> None:
    """Test that a damaged download raises and is deleted."""
    aioclient_mock.get(DE440_URL, content=content)

    with pytest.raises(EphemerisError, match="failed validation"):
        await async_get_ephemeris(hass)

    assert await hass.async_add_executor_job(_files, cache_dir) == []


async def test_cache_dir_inaccessible(hass: HomeAssistant, cache_dir: Path) -> None:
    """Test that an unusable cache directory raises."""
    cache_dir.write_bytes(b"")

    with pytest.raises(EphemerisError, match="Cannot access"):
        await async_get_ephemeris(hass)


@pytest.mark.usefixtures("jpl_server")
async def test_discard_ephemeris(hass: HomeAssistant, cache_dir: Path) -> None:
    """Test that discarding deletes the cached kernel and its directory."""
    first = await async_get_ephemeris(hass)

    await async_discard_ephemeris(hass)

    assert not cache_dir.exists()
    assert await async_get_ephemeris(hass) is not first


@pytest.mark.usefixtures("cached_kernel")
async def test_discard_ephemeris_keeps_other_files(
    hass: HomeAssistant, cache_dir: Path
) -> None:
    """Test that discarding keeps the directory when it holds other files."""
    (cache_dir / "finals2000A.all").write_bytes(b"")

    await async_discard_ephemeris(hass)

    assert await hass.async_add_executor_job(_files, cache_dir) == ["finals2000A.all"]


async def test_discard_ephemeris_cancels_download(
    hass: HomeAssistant, cache_dir: Path
) -> None:
    """Test that discarding cancels a download in progress."""
    started = asyncio.Event()

    async def _endless_download(*_args: Any) -> None:
        started.set()
        await asyncio.Event().wait()

    with patch(
        "custom_components.moon_astro.utils._async_download_kernel",
        side_effect=_endless_download,
    ):
        caller = hass.async_create_task(async_get_ephemeris(hass))
        await started.wait()

        await async_discard_ephemeris(hass)

    with pytest.raises(asyncio.CancelledError):
        await caller
    assert not cache_dir.exists()


@pytest.mark.parametrize(
    ("options", "expected"),
    [
        pytest.param({}, "Europe/Paris", id="default"),
        pytest.param(
            {CONF_USE_HA_TZ: True, CONF_TIME_ZONE: "Asia/Tokyo"},
            "Europe/Paris",
            id="home_assistant",
        ),
        pytest.param(
            {CONF_USE_HA_TZ: False, CONF_TIME_ZONE: "Asia/Tokyo"},
            "Asia/Tokyo",
            id="custom",
        ),
    ],
)
async def test_resolve_time_zone(
    hass: HomeAssistant, options: dict[str, Any], expected: str
) -> None:
    """Test the time zone naming full moons."""
    await hass.config.async_set_time_zone("Europe/Paris")
    entry = MockConfigEntry(domain=DOMAIN, options=options)

    assert str(await async_resolve_time_zone(entry)) == expected


@pytest.mark.parametrize(
    "options",
    [
        pytest.param({CONF_USE_HA_TZ: False}, id="missing"),
        pytest.param(
            {CONF_USE_HA_TZ: False, CONF_TIME_ZONE: "Mars/Olympus_Mons"}, id="unknown"
        ),
    ],
)
async def test_resolve_time_zone_fallback(
    hass: HomeAssistant, caplog: pytest.LogCaptureFixture, options: dict[str, Any]
) -> None:
    """Test that an unusable time zone falls back to that of Home Assistant."""
    await hass.config.async_set_time_zone("Europe/Paris")
    entry = MockConfigEntry(domain=DOMAIN, options=options)

    assert str(await async_resolve_time_zone(entry)) == "Europe/Paris"
    assert "is not available on this system" in caplog.text


async def test_available_time_zones(hass: HomeAssistant) -> None:
    """Test the list of time zones offered by the options flow."""
    zones = await hass.async_add_executor_job(available_time_zones)

    assert zones == sorted(zones)
    assert "Europe/Paris" in zones
