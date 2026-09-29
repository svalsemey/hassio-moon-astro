"""Coordinators for Moon Astro.

Two coordinators share the loaded DE440 kernel: one recomputes the current Moon
position on the configured scan interval, the other recomputes the values that only
change at astronomical events and schedules its own refresh after the next event.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, tzinfo
import logging
from typing import Any

from skyfield.api import wgs84
from skyfield.jpllib import SpiceKernel
from skyfield.timelib import Time, Timescale

from homeassistant.config_entries import ConfigEntry
from homeassistant.core import CALLBACK_TYPE, HomeAssistant
from homeassistant.helpers.event import async_track_point_in_time
from homeassistant.helpers.update_coordinator import DataUpdateCoordinator, UpdateFailed
from homeassistant.util import dt as dt_util

from .astronomy import (
    compute_current_payload,
    compute_events_payload,
    round_to_minute_utc,
)
from .const import (
    CONF_ALT,
    CONF_LAT,
    CONF_LON,
    KEY_NEXT_APOGEE,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_LAST_QUARTER,
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_PERIGEE,
)

_LOGGER = logging.getLogger(__name__)

# Errors turning a refresh into a failure: a kernel read error, or a Skyfield search or
# numeric failure such as an instant outside the ephemeris range.
_RECOVERABLE_UPDATE_ERRORS: tuple[type[Exception], ...] = (
    OSError,
    ValueError,
    ArithmeticError,
    RuntimeError,
)

# Margin after an event boundary before the event-based values are recomputed.
_BOUNDARY_REFRESH_DELAY = timedelta(minutes=2)

# Upcoming event keys driving the boundary refresh of the events coordinator.
_NEXT_EVENT_KEYS: tuple[str, ...] = (
    KEY_NEXT_NEW_MOON,
    KEY_NEXT_FIRST_QUARTER,
    KEY_NEXT_FULL_MOON,
    KEY_NEXT_LAST_QUARTER,
    KEY_NEXT_APOGEE,
    KEY_NEXT_PERIGEE,
)


class _MoonAstroBaseCoordinator(DataUpdateCoordinator[dict[str, Any]]):
    """Coordinator running a Skyfield computation in the executor at each refresh."""

    def __init__(
        self,
        hass: HomeAssistant,
        entry: MoonAstroConfigEntry,
        *,
        name: str,
        eph: SpiceKernel,
        ts: Timescale,
        interval: timedelta,
    ) -> None:
        """Initialize the coordinator.

        Args:
            hass: Home Assistant instance.
            entry: Config entry owning the coordinator.
            name: Coordinator name used in logs.
            eph: Loaded ephemeris shared by the coordinators.
            ts: Loaded timescale shared by the coordinators.
            interval: Update interval.
        """
        super().__init__(
            hass,
            logger=_LOGGER,
            config_entry=entry,
            name=name,
            update_interval=interval,
        )
        self._eph = eph
        self._ts = ts

    def _compute_payload(self, t: Time) -> dict[str, Any]:
        """Compute the payload for the reference time t; runs in the executor."""
        raise NotImplementedError

    async def _async_update_data(self) -> dict[str, Any]:
        """Compute the payload for the current minute.

        Returns:
            The payload exposed by the entities.

        Raises:
            UpdateFailed: If the computation fails.
        """
        t = self._ts.from_datetime(round_to_minute_utc(dt_util.utcnow()))
        try:
            return await self.hass.async_add_executor_job(self._compute_payload, t)
        except _RECOVERABLE_UPDATE_ERRORS as err:
            raise UpdateFailed(f"{self.name} computation failed: {err}") from err


class MoonAstroCoordinator(_MoonAstroBaseCoordinator):
    """Coordinator computing the current Moon position and the surrounding rises and sets."""

    def __init__(
        self,
        hass: HomeAssistant,
        entry: MoonAstroConfigEntry,
        *,
        eph: SpiceKernel,
        ts: Timescale,
        interval: timedelta,
    ) -> None:
        """Initialize the coordinator with the observer location of the entry.

        Args:
            hass: Home Assistant instance.
            entry: Config entry providing the observer coordinates.
            eph: Loaded ephemeris shared by the coordinators.
            ts: Loaded timescale shared by the coordinators.
            interval: Update interval.
        """
        super().__init__(
            hass, entry, name="Moon Astro", eph=eph, ts=ts, interval=interval
        )
        self._observer = wgs84.latlon(
            latitude_degrees=float(entry.data.get(CONF_LAT, hass.config.latitude)),
            longitude_degrees=float(entry.data.get(CONF_LON, hass.config.longitude)),
            elevation_m=float(entry.data.get(CONF_ALT, hass.config.elevation)),
        )

    def _compute_payload(self, t: Time) -> dict[str, Any]:
        """Compute the current position payload; runs in the executor."""
        return compute_current_payload(self._eph, t, self._observer)


class MoonAstroEventsCoordinator(_MoonAstroBaseCoordinator):
    """Coordinator computing the values that only change at astronomical events.

    Besides the periodic fallback refresh, a refresh is scheduled shortly after the
    earliest upcoming event so that the previous and next values switch on time.
    """

    def __init__(
        self,
        hass: HomeAssistant,
        entry: MoonAstroConfigEntry,
        *,
        eph: SpiceKernel,
        ts: Timescale,
        interval: timedelta,
        tz: tzinfo,
    ) -> None:
        """Initialize the coordinator.

        Args:
            hass: Home Assistant instance.
            entry: Config entry owning the coordinator.
            eph: Loaded ephemeris shared by the coordinators.
            ts: Loaded timescale shared by the coordinators.
            interval: Fallback update interval.
            tz: Time zone defining the calendar months that name full moons.
        """
        super().__init__(
            hass, entry, name="Moon Astro Events", eph=eph, ts=ts, interval=interval
        )
        self._tz = tz
        self._unsub_boundary: CALLBACK_TYPE | None = None
        self.next_refresh_utc: datetime | None = None

    def _compute_payload(self, t: Time) -> dict[str, Any]:
        """Compute the event-based payload; runs in the executor."""
        return compute_events_payload(self._eph, t, self._tz)

    async def _async_update_data(self) -> dict[str, Any]:
        """Compute the payload and schedule the refresh following the next event."""
        data = await super()._async_update_data()
        self._schedule_boundary_refresh(data)
        return data

    def _schedule_boundary_refresh(self, data: dict[str, Any]) -> None:
        """Schedule a refresh shortly after the earliest upcoming event.

        Also records the instant of the next refresh, whether it comes from the event
        boundary or from the periodic fallback interval.

        Args:
            data: Freshly computed payload.
        """
        self._cancel_boundary_refresh()
        fallback = (
            None
            if (interval := self.update_interval) is None
            else dt_util.utcnow() + interval
        )
        upcoming = [dt for key in _NEXT_EVENT_KEYS if (dt := data[key]) is not None]
        boundary = None if not upcoming else min(upcoming) + _BOUNDARY_REFRESH_DELAY
        self.next_refresh_utc = min(
            (dt for dt in (fallback, boundary) if dt is not None), default=None
        )
        if boundary is not None:
            self._unsub_boundary = async_track_point_in_time(
                self.hass, self._async_refresh_at_boundary, boundary
            )

    async def _async_refresh_at_boundary(self, _now: datetime) -> None:
        """Request a refresh once the event boundary has passed."""
        self._unsub_boundary = None
        await self.async_request_refresh()

    def _cancel_boundary_refresh(self) -> None:
        """Cancel the pending boundary refresh, if any."""
        if self._unsub_boundary is not None:
            self._unsub_boundary()
            self._unsub_boundary = None

    async def async_shutdown(self) -> None:
        """Cancel the boundary refresh and stop the coordinator."""
        self._cancel_boundary_refresh()
        await super().async_shutdown()


@dataclass(frozen=True, slots=True)
class MoonAstroRuntimeData:
    """Runtime objects attached to a loaded config entry."""

    coordinator: MoonAstroCoordinator
    events_coordinator: MoonAstroEventsCoordinator


type MoonAstroConfigEntry = ConfigEntry[MoonAstroRuntimeData]
