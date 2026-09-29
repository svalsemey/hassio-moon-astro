"""Ephemeris lifecycle and configuration helpers for Moon Astro.

The DE440 kernel is downloaded once into the Home Assistant configuration
directory, validated, loaded a single time per Home Assistant instance and then
shared by every coordinator. Filesystem and Skyfield calls run in the executor.
"""

from __future__ import annotations

import asyncio
from contextlib import suppress
from datetime import tzinfo
from functools import cache, partial
import logging
from pathlib import Path
import zoneinfo

import aiohttp
from skyfield.api import Loader
from skyfield.jpllib import SpiceKernel
from skyfield.timelib import Timescale

from homeassistant.config_entries import ConfigEntry
from homeassistant.core import HomeAssistant
from homeassistant.exceptions import HomeAssistantError
from homeassistant.helpers.aiohttp_client import async_get_clientsession
from homeassistant.util import dt as dt_util
from homeassistant.util.hass_dict import HassKey

from .const import (
    CACHE_DIR_NAME,
    CONF_TIME_ZONE,
    CONF_USE_HA_TZ,
    DE440_FILE,
    DE440_URL,
    DEFAULT_USE_HA_TZ,
    DOMAIN,
    MIN_EPHEMERIS_SIZE_BYTES,
)

_LOGGER = logging.getLogger(__name__)

type EphemerisKernel = tuple[SpiceKernel, Timescale]

_SHARED_KERNEL_KEY: HassKey[EphemerisKernel] = HassKey(f"{DOMAIN}_ephemeris_kernel")
_PREPARE_TASK_KEY: HassKey[asyncio.Task[EphemerisKernel]] = HassKey(
    f"{DOMAIN}_ephemeris_task"
)

_DOWNLOAD_CHUNK_BYTES = 1024 * 1024
# No overall limit: the 115 MB transfer may legitimately take many minutes.
_DOWNLOAD_TIMEOUT = aiohttp.ClientTimeout(total=None, connect=30, sock_read=60)
# Errors raised by Skyfield and jplephem when a kernel file is damaged.
_KERNEL_LOAD_ERRORS: tuple[type[Exception], ...] = (
    OSError,
    ValueError,
    KeyError,
    IndexError,
    TypeError,
)


class EphemerisError(HomeAssistantError):
    """Raised when the DE440 ephemeris cannot be made available."""


def get_cache_dir(hass: HomeAssistant) -> Path:
    """Return the directory holding the Skyfield resources.

    Args:
        hass: Home Assistant instance.

    Returns:
        Path to the cache directory.
    """
    return Path(hass.config.path(CACHE_DIR_NAME))


def _unlink_quietly(path: Path) -> None:
    """Delete a file, ignoring a missing file or a filesystem error."""
    with suppress(OSError):
        path.unlink()


def _open_kernel(cache_dir: Path) -> EphemerisKernel:
    """Open the cached kernel and check that it is usable.

    Runs in the executor.

    Args:
        cache_dir: Directory holding the kernel file.

    Returns:
        The loaded kernel and timescale.

    Raises:
        FileNotFoundError: If the kernel file is missing.
        EphemerisError: If the file is truncated or unreadable.
    """
    path = cache_dir / DE440_FILE
    if (size := path.stat().st_size) < MIN_EPHEMERIS_SIZE_BYTES:
        raise EphemerisError(f"{path} is truncated ({size} bytes)")
    try:
        timescale = Loader(str(cache_dir), verbose=False).timescale()
        kernel = SpiceKernel(str(path))
        # Exercise every segment used at runtime so that a damaged file is rejected.
        kernel["earth"].at(timescale.now()).observe(
            kernel["moon"]
        ).apparent().fraction_illuminated(kernel["sun"])
    except _KERNEL_LOAD_ERRORS as err:
        raise EphemerisError(f"{path} is unreadable: {err}") from err
    return kernel, timescale


def _load_kernel(cache_dir: Path) -> EphemerisKernel | None:
    """Return the cached kernel, or None when it is missing or was found damaged.

    A damaged file is deleted so that a fresh copy gets downloaded. Runs in the
    executor.

    Args:
        cache_dir: Directory holding the kernel file.

    Returns:
        The loaded kernel and timescale, or None.

    Raises:
        EphemerisError: If the cache directory cannot be accessed.
    """
    try:
        kernel = _open_kernel(cache_dir)
    except FileNotFoundError:
        return None
    except OSError as err:
        raise EphemerisError(f"Cannot access {cache_dir / DE440_FILE}: {err}") from err
    except EphemerisError as err:
        _LOGGER.warning("Discarding the ephemeris file: %s", err)
        _unlink_quietly(cache_dir / DE440_FILE)
        return None
    return kernel


async def _async_download_kernel(hass: HomeAssistant, cache_dir: Path) -> None:
    """Stream the DE440 kernel from JPL into the cache directory.

    The data goes through a temporary file renamed once complete, so that an
    interrupted transfer never leaves a partial kernel in place.

    Args:
        hass: Home Assistant instance.
        cache_dir: Destination directory, created when missing.

    Raises:
        EphemerisError: If the download fails.
    """
    path = cache_dir / DE440_FILE
    temp_path = path.with_name(f"{DE440_FILE}.download")
    try:
        await hass.async_add_executor_job(
            partial(cache_dir.mkdir, parents=True, exist_ok=True)
        )
        async with async_get_clientsession(hass).get(
            DE440_URL, timeout=_DOWNLOAD_TIMEOUT
        ) as response:
            response.raise_for_status()
            handle = await hass.async_add_executor_job(temp_path.open, "wb")
            try:
                async for chunk in response.content.iter_chunked(_DOWNLOAD_CHUNK_BYTES):
                    await hass.async_add_executor_job(handle.write, chunk)
            finally:
                await hass.async_add_executor_job(handle.close)
        await hass.async_add_executor_job(temp_path.replace, path)
    except (aiohttp.ClientError, TimeoutError, OSError) as err:
        await hass.async_add_executor_job(_unlink_quietly, temp_path)
        raise EphemerisError(f"Download of {DE440_URL} failed: {err}") from err


async def _async_prepare_kernel(hass: HomeAssistant) -> EphemerisKernel:
    """Load the cached kernel, downloading a fresh copy first when needed.

    Args:
        hass: Home Assistant instance.

    Returns:
        The loaded kernel and timescale.

    Raises:
        EphemerisError: If no usable kernel can be obtained.
    """
    cache_dir = get_cache_dir(hass)
    if (kernel := await hass.async_add_executor_job(_load_kernel, cache_dir)) is None:
        _LOGGER.info(
            "Downloading the DE440 ephemeris (about 115 MB) from %s", DE440_URL
        )
        await _async_download_kernel(hass, cache_dir)
        if (
            kernel := await hass.async_add_executor_job(_load_kernel, cache_dir)
        ) is None:
            raise EphemerisError("The downloaded ephemeris file failed validation")
        _LOGGER.info("DE440 ephemeris stored in %s", cache_dir)
    hass.data[_SHARED_KERNEL_KEY] = kernel
    return kernel


async def async_get_ephemeris(hass: HomeAssistant) -> EphemerisKernel:
    """Return the shared ephemeris kernel and timescale, preparing them on first use.

    Concurrent callers (config flow, entry setup) share a single preparation task.
    That task is shielded so that a cancelled caller, typically an abandoned config
    flow, does not abort a download another caller may be waiting for.

    Args:
        hass: Home Assistant instance.

    Returns:
        The loaded kernel and timescale.

    Raises:
        EphemerisError: If the kernel cannot be made available.
    """
    if (kernel := hass.data.get(_SHARED_KERNEL_KEY)) is not None:
        return kernel
    if (task := hass.data.get(_PREPARE_TASK_KEY)) is None or task.done():
        task = hass.data[_PREPARE_TASK_KEY] = hass.async_create_background_task(
            _async_prepare_kernel(hass), name=f"{DOMAIN}-ephemeris"
        )
    return await asyncio.shield(task)


async def async_discard_ephemeris(hass: HomeAssistant) -> None:
    """Drop the shared kernel and delete the cached files.

    A preparation still in flight is cancelled and awaited first so that no file
    reappears after the cleanup.

    Args:
        hass: Home Assistant instance.
    """
    if (task := hass.data.pop(_PREPARE_TASK_KEY, None)) is not None and not task.done():
        task.cancel()
        await asyncio.wait([task])
    hass.data.pop(_SHARED_KERNEL_KEY, None)
    cache_dir = get_cache_dir(hass)

    def _remove() -> None:
        """Delete the kernel and any partial download, then the empty directory."""
        for name in (DE440_FILE, f"{DE440_FILE}.download"):
            _unlink_quietly(cache_dir / name)
        # Only succeeds when nothing else is stored in the directory.
        with suppress(OSError):
            cache_dir.rmdir()

    await hass.async_add_executor_job(_remove)


@cache
def available_time_zones() -> tuple[str, ...]:
    """Return the sorted IANA time zone names available on this system.

    The first call walks the tzdata directories and must run in the executor.

    Returns:
        Sorted time zone names.
    """
    return tuple(sorted(zoneinfo.available_timezones()))


async def async_resolve_time_zone(entry: ConfigEntry) -> tzinfo:
    """Return the time zone used for calendar-based values of a config entry.

    The Home Assistant time zone is used unless the entry explicitly opts out and
    provides a valid IANA time zone name. An unresolvable name falls back to the
    Home Assistant time zone so that the coordinators always have a zone.

    Args:
        entry: Config entry providing the time zone options.

    Returns:
        A tzinfo instance.
    """
    if entry.options.get(CONF_USE_HA_TZ, DEFAULT_USE_HA_TZ):
        return dt_util.get_default_time_zone()

    name = entry.options.get(CONF_TIME_ZONE)
    if name and (tz := await dt_util.async_get_time_zone(name)) is not None:
        return tz

    _LOGGER.warning(
        "Time zone %r is not available on this system; using the Home Assistant time zone",
        name,
    )
    return dt_util.get_default_time_zone()
