# Moon Astro

High-precision Moon ephemeris integration for Home Assistant, powered by Skyfield (DE440). Provides current topocentric/geocentric ecliptic coordinates, illumination, distance, parallax, moonrise/set, next lunation timestamps, apogee/perigee, and zodiac information.

- Accurate ecliptic-of-date conversion (IAU 1980 nutation, true obliquity)
- Topocentric elevation/azimuth from your configured location
- Fully localized entities and state translations
- Config flow with options for scan interval and time zone handling

## Features

- Current:
  - Phase, Azimuth, Elevation, Illumination (%), Distance (km), Parallax (°)
  - Ecliptic longitude/latitude (topocentric and geocentric)
  - Zodiac sign and zodiac degree (current moon position)
- Next events:
  - Moonrise, Moonset
  - Apogee, Perigee
  - New Moon, First Quarter, Full Moon, Last Quarter
- Lunation context:
  - Geocentric ecliptic longitude/latitude at next new/full moon (true-of-date)
  - Zodiac sign and zodiac degree at next new/full moon
  - Next full moon name and usual alternative names (Gregorian calendar)
- Binary sensor:
  - Moon above horizon (on/off)

## Installation

HACS (recommended):
1. In HACS → Integrations → Custom repositories, add `https://github.com/svalsemey/hassio-moon-astro/` as type “Integration”.
2. Install “Moon Astro”.
3. Restart Home Assistant and add the integration from Devices & Services.

Manual install:
1. Download or clone the repository.
2. Copy the folder `custom_components/moon_astro` to your Home Assistant config directory:
   - Final path: `<config>/custom_components/moon_astro`
3. Restart Home Assistant.
4. Go to Settings → Devices & Services → Add Integration → search for “Moon Astro”.

## Configuration

During setup, the integration pre-fills your Home Assistant latitude, longitude, and elevation. You can adjust them if needed.

Options (via Configure on the integration):
- Scan interval (seconds). Default: 300 (min 30, max 21600)
- Use Home Assistant time zone. Default: true

## Localization

Translations are included for multiple languages. Entity names and some state values (for example zodiac signs and full moon names) are localized via Home Assistant’s translation system.

If you notice that a translation is incomplete or inaccurate, contributions are very welcome.

## Requirements

- Home Assistant 2023.6+ (recommended)
- Python dependencies are installed automatically by HA:
  - skyfield>=1.53
  - timezonefinder>=6.2.0

Note: Skyfield downloads ephemeris/timescale data to `<config>/.skyfield` on first run. These files are cached for later use.

## Privacy and Network

- No external network calls during runtime except Skyfield’s initial ephemeris download (cached locally).
- Computations are local on your Home Assistant instance.

## Troubleshooting

- Entities show “unavailable”:
  - Check Logs → “Moon Astro” for errors.
  - Ensure your HA location and elevation are set.
  - Verify ephemeris download completed (see `<config>/.skyfield` content).
- Time zone issues:
  - Toggle the “Use Home Assistant time zone” option or verify your HA system time zone.
- Data looks stale:
  - Reduce the scan interval (respect CPU usage).
  - Manually reload the integration from the UI.

## Known Notes

- Ecliptic-of-date values use IAU 1980 nutation with true obliquity.
- If you modify locally translation files, reload the integration to apply changes.

## Contributing

Issues and PRs are welcome!
- Issues: https://github.com/svalsemey/hassio-moon-astro/issues
- Pull requests: please lint and follow HA integration best practices.

## License

MIT © Sébastien VALSEMEY

## Credits

- Skyfield by Brandon Rhodes
- Home Assistant core team and community
