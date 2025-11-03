"""Config and Options flow for Moon Astro."""

from __future__ import annotations

import voluptuous as vol

from homeassistant import config_entries
from homeassistant.core import callback
from homeassistant.data_entry_flow import FlowResult

from .const import (
    CONF_ALT,
    CONF_LAT,
    CONF_LON,
    CONF_SCAN_INTERVAL,
    CONF_USE_HA_TZ,
    DEFAULT_SCAN_INTERVAL,
    DOMAIN,
)


class MoonAstroConfigFlow(config_entries.ConfigFlow, domain=DOMAIN):
    """Handle a config flow for Moon Astro."""

    VERSION = 1

    async def async_step_user(self, user_input: dict | None = None) -> FlowResult:
        """Handle the initial step for GPS coordinates and altitude."""
        errors: dict[str, str] = {}

        if user_input is not None:
            await self.async_set_unique_id(DOMAIN)
            self._abort_if_unique_id_configured()
            return self.async_create_entry(title="Moon Astro", data=user_input)

        hass_lat = self.hass.config.latitude
        hass_lon = self.hass.config.longitude
        hass_elev = self.hass.config.elevation or 0

        schema = vol.Schema(
            {
                vol.Required(CONF_LAT, default=hass_lat): vol.Coerce(float),
                vol.Required(CONF_LON, default=hass_lon): vol.Coerce(float),
                vol.Optional(CONF_ALT, default=hass_elev): vol.Coerce(float),
            }
        )
        return self.async_show_form(step_id="user", data_schema=schema, errors=errors)

    async def async_step_import(self, user_input: dict) -> FlowResult:
        """Support YAML import if needed in the future."""
        return await self.async_step_user(user_input)

    async def async_step_reconfigure(
        self, user_input: dict | None = None
    ) -> FlowResult:
        """Handle reconfiguration (reserved for future)."""
        return await self.async_step_user(user_input)

    @staticmethod
    @callback
    def async_get_options_flow(config_entry: config_entries.ConfigEntry):
        """Return the options flow handler."""
        return MoonAstroOptionsFlow(config_entry)


class MoonAstroOptionsFlow(config_entries.OptionsFlow):
    """Handle options for Moon Astro."""

    def __init__(self, entry: config_entries.ConfigEntry) -> None:
        self._entry = entry

    async def async_step_init(self, user_input: dict | None = None) -> FlowResult:
        """First step of options flow."""
        if user_input is not None:
            return self.async_create_entry(title="", data=user_input)

        options = self._entry.options
        schema = vol.Schema(
            {
                vol.Optional(
                    CONF_SCAN_INTERVAL,
                    default=options.get(CONF_SCAN_INTERVAL, DEFAULT_SCAN_INTERVAL),
                ): vol.All(vol.Coerce(int), vol.Range(min=30, max=21600)),
                vol.Optional(
                    CONF_USE_HA_TZ,
                    default=options.get(CONF_USE_HA_TZ, True),
                ): bool,
            }
        )
        return self.async_show_form(step_id="init", data_schema=schema)
