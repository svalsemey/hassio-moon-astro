"""Utility functions for Moon Astro integration.

This module centralizes helpers for ephemeris lifecycle management.
All filesystem and Skyfield operations are executed in an executor to avoid
blocking the event loop.
"""

from __future__ import annotations

from contextlib import suppress
from pathlib import Path

from homeassistant.core import HomeAssistant

from .const import CACHE_DIR_NAME, DE440_FILE

try:
    from skyfield.api import Loader
    from skyfield.jpllib import SpiceKernel

    SKYFIELD_AVAILABLE = True
except ImportError:  # pragma: no cover
    SKYFIELD_AVAILABLE = False


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
        cache_dir = Path(hass.config.path(CACHE_DIR_NAME))
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
                remaining_files = list(cache_dir.iterdir())
                if not remaining_files:
                    cache_dir.rmdir()

    await hass.async_add_executor_job(_blocking_cleanup)


async def validate_ephemeris_file(hass: HomeAssistant) -> bool:
    """Validate the integrity of the ephemeris file.

    The validation is intentionally strict:
    - checks file existence
    - checks a conservative minimum size threshold
    - attempts to load the file with Skyfield
    - checks that essential bodies are present

    If the file is found invalid, it is removed so a later download can start
    from a clean state.

    Args:
        hass: Home Assistant instance.

    Returns:
        True if file is valid and usable, False if invalid or missing.
    """

    def _blocking_validation() -> bool:
        """Execute validation in executor to avoid blocking."""
        if not SKYFIELD_AVAILABLE:
            return False

        cache_dir = Path(hass.config.path(CACHE_DIR_NAME))
        ephemeris_path = cache_dir / DE440_FILE

        if not ephemeris_path.exists():
            return False

        # de440.bsp is large; use a conservative minimum threshold.
        min_expected_size = 100 * 1024 * 1024
        try:
            file_size = ephemeris_path.stat().st_size
            if file_size < min_expected_size:
                with suppress(OSError):
                    ephemeris_path.unlink()
                return False
        except OSError:
            return False

        try:
            loader = Loader(str(cache_dir))
            eph = loader(DE440_FILE)
        except (OSError, AttributeError, RuntimeError):
            with suppress(OSError):
                ephemeris_path.unlink()
            return False

        if not isinstance(eph, SpiceKernel):
            with suppress(OSError):
                ephemeris_path.unlink()
            return False

        # Ensure essential bodies exist: solar system barycenter, Earth, Moon.
        if 0 not in eph or 3 not in eph or 301 not in eph:
            with suppress(OSError):
                ephemeris_path.unlink()
            return False

        return True

    return await hass.async_add_executor_job(_blocking_validation)


async def ensure_valid_ephemeris(hass: HomeAssistant) -> bool:
    """Ensure a valid ephemeris file exists, downloading only if necessary.

    This function is designed to be idempotent:
    - If a valid file already exists, it returns True without downloading.
    - If the file is missing or invalid, it tries to download it and then validates it.

    A best-effort cleanup of temporary download artifacts is done before checking
    the local file, to avoid leaving stale partial downloads around.

    Args:
        hass: Home Assistant instance.

    Returns:
        True if a valid ephemeris is available, False otherwise.
    """
    if not SKYFIELD_AVAILABLE:
        return False

    # Remove stale temporary download files first to keep the directory consistent.
    await cleanup_cache_dir(hass)

    # Fast path: if the file is already present and valid, do not download.
    if await validate_ephemeris_file(hass):
        return True

    def _blocking_download() -> bool:
        """Download ephemeris file in executor.

        Returns:
            True if the download attempt did not raise, False otherwise.
        """
        cache_dir = Path(hass.config.path(CACHE_DIR_NAME))
        cache_dir.mkdir(parents=True, exist_ok=True)

        try:
            loader = Loader(str(cache_dir))
            loader(DE440_FILE)
        except (OSError, RuntimeError, ValueError):
            return False

        return True

    download_ok = await hass.async_add_executor_job(_blocking_download)
    if not download_ok:
        return False

    # Final authoritative validation (also deletes corrupted files).
    return await validate_ephemeris_file(hass)
