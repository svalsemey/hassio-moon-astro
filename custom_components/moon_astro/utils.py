"""Utility functions for Moon Astro integration.

This module centralizes helpers for ephemeris lifecycle management.
All filesystem and Skyfield operations are executed in an executor to avoid
blocking the event loop.
"""

from __future__ import annotations

from contextlib import suppress
import logging
from pathlib import Path

from homeassistant.core import HomeAssistant

from .const import CACHE_DIR_NAME, DE440_FILE

_LOGGER = logging.getLogger(__name__)

try:
    from skyfield.api import Loader
    from skyfield.jpllib import SpiceKernel

    SKYFIELD_AVAILABLE = True
except ImportError:  # pragma: no cover
    SKYFIELD_AVAILABLE = False


def get_cache_dir(hass: HomeAssistant) -> Path:
    """Return the cache directory used to store Skyfield resources.

    Args:
        hass: Home Assistant instance.

    Returns:
        Path to the cache directory.
    """
    return Path(hass.config.path(CACHE_DIR_NAME))


def get_ephemeris_path(hass: HomeAssistant) -> Path:
    """Return the expected full path of the ephemeris file.

    Args:
        hass: Home Assistant instance.

    Returns:
        Path to the ephemeris file.
    """
    return get_cache_dir(hass) / DE440_FILE


async def cleanup_cache_dir(
    hass: HomeAssistant,
    *,
    remove_empty_dir: bool = False,
    remove_ephemeris: bool = False,
) -> None:
    """Clean up Skyfield cache content.

    Args:
        hass: Home Assistant instance.
        remove_empty_dir: If True, remove the cache directory if empty after cleanup.
        remove_ephemeris: If True, also remove the main ephemeris file.

    Returns:
        None.
    """

    def _blocking_cleanup() -> None:
        """Execute cleanup in executor to avoid blocking the event loop."""
        cache_dir = get_cache_dir(hass)
        if not cache_dir.exists():
            return

        with suppress(OSError):
            for temp_file in cache_dir.glob("*.download*"):
                with suppress(OSError):
                    temp_file.unlink()

        if remove_ephemeris:
            ephemeris_path = cache_dir / DE440_FILE
            with suppress(OSError):
                ephemeris_path.unlink()

        if remove_empty_dir:
            with suppress(OSError):
                if not any(cache_dir.iterdir()):
                    cache_dir.rmdir()

    await hass.async_add_executor_job(_blocking_cleanup)


async def validate_ephemeris_file(
    hass: HomeAssistant,
    *,
    remove_on_invalid: bool = False,
) -> bool:
    """Validate the integrity of the ephemeris file.

    The validation checks:
    - file existence
    - conservative minimum size threshold
    - ability to load it with Skyfield
    - expected kernel type
    - presence of essential bodies

    Args:
        hass: Home Assistant instance.
        remove_on_invalid: If True, delete the ephemeris file when invalid.

    Returns:
        True if file is valid and usable, False if invalid or missing.
    """

    def _blocking_validation() -> bool:
        """Execute validation in executor to avoid blocking."""
        if not SKYFIELD_AVAILABLE:
            return False

        cache_dir = get_cache_dir(hass)
        ephemeris_path = cache_dir / DE440_FILE

        if not ephemeris_path.exists():
            return False

        min_expected_size = 100 * 1024 * 1024
        try:
            if ephemeris_path.stat().st_size < min_expected_size:
                if remove_on_invalid:
                    with suppress(OSError):
                        ephemeris_path.unlink()
                return False
        except OSError:
            return False

        try:
            loader = Loader(str(cache_dir))
            eph = loader(DE440_FILE)
        except (OSError, AttributeError, RuntimeError):
            if remove_on_invalid:
                with suppress(OSError):
                    ephemeris_path.unlink()
            return False

        if not isinstance(eph, SpiceKernel):
            if remove_on_invalid:
                with suppress(OSError):
                    ephemeris_path.unlink()
            return False

        if 0 not in eph or 3 not in eph or 301 not in eph:
            if remove_on_invalid:
                with suppress(OSError):
                    ephemeris_path.unlink()
            return False

        return True

    return await hass.async_add_executor_job(_blocking_validation)


async def ensure_valid_ephemeris(hass: HomeAssistant) -> bool:
    """Ensure a valid ephemeris file exists, downloading only if necessary.

    Args:
        hass: Home Assistant instance.

    Returns:
        True if a valid ephemeris is available, False otherwise.
    """
    if not SKYFIELD_AVAILABLE:
        return False

    await cleanup_cache_dir(hass)

    if await validate_ephemeris_file(hass, remove_on_invalid=False):
        return True

    def _blocking_download() -> bool:
        """Download ephemeris file in executor.

        Returns:
            True if the download attempt did not raise and the expected file exists.
        """
        cache_dir = get_cache_dir(hass)
        cache_dir.mkdir(parents=True, exist_ok=True)
        ephemeris_path = cache_dir / DE440_FILE

        try:
            loader = Loader(str(cache_dir))
            loader(DE440_FILE)
        except (OSError, RuntimeError, ValueError) as exc:
            _LOGGER.info("Ephemeris download failed: %r", exc)
            return False

        # Ensure the file exists where the integration expects it.
        if not ephemeris_path.exists():
            # Some environments may redirect cache paths; log to diagnose quickly.
            try:
                actual = None
                if hasattr(loader, "path_to"):
                    actual = loader.path_to(DE440_FILE)
                _LOGGER.info(
                    "Ephemeris download completed but file not found at expected path: expected=%s actual=%s",
                    str(ephemeris_path),
                    str(actual) if actual is not None else "unknown",
                )
            except Exception as exc:  # pragma: no cover
                _LOGGER.info(
                    "Ephemeris download completed but file not found at expected path: expected=%s (error=%r)",
                    str(ephemeris_path),
                    exc,
                )
            return False

        _LOGGER.info("Ephemeris file ready in cache: %s", str(ephemeris_path))
        return True

    if not await hass.async_add_executor_job(_blocking_download):
        return False

    return await validate_ephemeris_file(hass, remove_on_invalid=False)
