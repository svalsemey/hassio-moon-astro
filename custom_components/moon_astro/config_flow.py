"""Config and options flows for Moon Astro."""

from __future__ import annotations

import asyncio
from collections.abc import Mapping
from typing import Any

import voluptuous as vol

from homeassistant.config_entries import (
    ConfigEntry,
    ConfigFlow,
    ConfigFlowResult,
    OptionsFlow,
)
from homeassistant.core import callback
from homeassistant.helpers import config_validation as cv
from homeassistant.helpers.selector import (
    SelectSelector,
    SelectSelectorConfig,
    SelectSelectorMode,
)

from .const import (
    CONF_ALT,
    CONF_EVENTS_REFRESH_FALLBACK,
    CONF_HIGH_PRECISION,
    CONF_LAT,
    CONF_LON,
    CONF_SCAN_INTERVAL,
    CONF_TIME_ZONE,
    CONF_USE_HA_TZ,
    DEFAULT_EVENTS_REFRESH_FALLBACK,
    DEFAULT_HIGH_PRECISION,
    DEFAULT_SCAN_INTERVAL,
    DEFAULT_USE_HA_TZ,
    DOMAIN,
    NAME,
)
from .utils import (
    EphemerisError,
    EphemerisKernel,
    async_get_ephemeris,
    available_time_zones,
)

LOCATION_SCHEMA = vol.Schema(
    {
        vol.Required(CONF_LAT): cv.latitude,
        vol.Required(CONF_LON): cv.longitude,
        vol.Required(CONF_ALT): vol.All(
            vol.Coerce(float), vol.Range(min=-500, max=10000)
        ),
    }
)
PRECISION_SCHEMA = vol.Schema(
    {vol.Required(CONF_HIGH_PRECISION, default=DEFAULT_HIGH_PRECISION): bool}
)


class MoonAstroConfigFlow(ConfigFlow, domain=DOMAIN):
    """Handle the initial configuration and the reconfiguration of Moon Astro."""

    VERSION = 1

    def __init__(self) -> None:
        """Initialize the flow state."""
        self._location: dict[str, Any] = {}
        self._download_task: asyncio.Task[EphemerisKernel] | None = None

    def _show_location_form(
        self,
        step_id: str,
        suggested: Mapping[str, Any],
        errors: dict[str, str] | None = None,
    ) -> ConfigFlowResult:
        """Show the observer location form pre-filled with suggested values.

        Args:
            step_id: Step handling the submitted form.
            suggested: Values pre-filling the fields.
            errors: Form errors to display.

        Returns:
            The form result.
        """
        return self.async_show_form(
            step_id=step_id,
            data_schema=self.add_suggested_values_to_schema(LOCATION_SCHEMA, suggested),
            errors=errors,
        )

    async def async_step_user(
        self, user_input: dict[str, Any] | None = None
    ) -> ConfigFlowResult:
        """Collect the observer location, then move on to the ephemeris download."""
        await self.async_set_unique_id(DOMAIN)
        self._abort_if_unique_id_configured()
        if user_input is None:
            return self._show_location_form(
                "user",
                {
                    CONF_LAT: self.hass.config.latitude,
                    CONF_LON: self.hass.config.longitude,
                    CONF_ALT: self.hass.config.elevation,
                },
            )
        self._location = user_input
        return await self.async_step_download()

    async def async_step_download(
        self, user_input: dict[str, Any] | None = None
    ) -> ConfigFlowResult:
        """Prepare the ephemeris behind a progress dialog.

        The flow manager calls this step again once the task completes.
        """
        if self._download_task is None:
            # A deferred start guarantees that the progress dialog is shown at
            # least once, even when the ephemeris is already available.
            self._download_task = self.hass.async_create_task(
                async_get_ephemeris(self.hass), eager_start=False
            )
        if not self._download_task.done():
            return self.async_show_progress(
                step_id="download",
                progress_action="download_ephemeris",
                progress_task=self._download_task,
            )
        try:
            await self._download_task
        except EphemerisError:
            next_step_id = "download_failed"
        else:
            next_step_id = "precision"
        finally:
            self._download_task = None
        return self.async_show_progress_done(next_step_id=next_step_id)

    async def async_step_download_failed(
        self, user_input: dict[str, Any] | None = None
    ) -> ConfigFlowResult:
        """Show the location form again with the download error."""
        return self._show_location_form(
            "user", self._location, errors={"base": "download_failed"}
        )

    async def async_step_precision(
        self, user_input: dict[str, Any] | None = None
    ) -> ConfigFlowResult:
        """Collect the precision option and create the entry."""
        if user_input is None:
            return self.async_show_form(
                step_id="precision", data_schema=PRECISION_SCHEMA
            )
        return self.async_create_entry(
            title=NAME, data=self._location, options=user_input
        )

    async def async_step_reconfigure(
        self, user_input: dict[str, Any] | None = None
    ) -> ConfigFlowResult:
        """Update the observer location of the existing entry."""
        entry = self._get_reconfigure_entry()
        if user_input is None:
            return self._show_location_form("reconfigure", entry.data)
        return self.async_update_reload_and_abort(entry, data_updates=user_input)

    @staticmethod
    @callback
    def async_get_options_flow(config_entry: ConfigEntry) -> OptionsFlow:
        """Return the options flow handler."""
        return MoonAstroOptionsFlow()


class MoonAstroOptionsFlow(OptionsFlow):
    """Handle the Moon Astro options."""

    async def async_step_init(
        self, user_input: dict[str, Any] | None = None
    ) -> ConfigFlowResult:
        """Show and store the options."""
        if user_input is not None:
            return self.async_create_entry(data=user_input)

        options = self.config_entry.options
        time_zones = await self.hass.async_add_executor_job(available_time_zones)
        schema = vol.Schema(
            {
                vol.Required(
                    CONF_SCAN_INTERVAL,
                    default=options.get(CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL),
                ): vol.All(vol.Coerce(int), vol.Range(min=60, max=21600)),
                vol.Required(
                    CONF_USE_HA_TZ,
                    default=options.get(CONF_USE_HA_TZ, DEFAULT_USE_HA_TZ),
                ): bool,
                vol.Required(
                    CONF_TIME_ZONE,
                    default=options.get(CONF_TIME_ZONE, self.hass.config.time_zone),
                ): SelectSelector(
                    SelectSelectorConfig(
                        options=time_zones, mode=SelectSelectorMode.DROPDOWN
                    )
                ),
                vol.Required(
                    CONF_HIGH_PRECISION,
                    default=options.get(CONF_HIGH_PRECISION, DEFAULT_HIGH_PRECISION),
                ): bool,
                vol.Required(
                    CONF_EVENTS_REFRESH_FALLBACK,
                    default=options.get(
                        CONF_EVENTS_REFRESH_FALLBACK, DEFAULT_EVENTS_REFRESH_FALLBACK
                    ),
                ): vol.All(vol.Coerce(int), vol.Range(min=3600, max=604800)),
            }
        )
        return self.async_show_form(step_id="init", data_schema=schema)
