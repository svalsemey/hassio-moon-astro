# Moon Astro

High-precision Moon ephemeris integration for Home Assistant, powered by Skyfield (DE440). Provides current topocentric/geocentric ecliptic coordinates, illumination, distance, parallax, moonrise/set, next lunation timestamps, apogee/perigee, and zodiac information at the next new/full moon.

- Accurate ecliptic-of-date conversion (IAU 1980 nutation, true obliquity)
- Topocentric Alt/Az from your configured location
- Fully localized entities and state translations
- Config flow with options for scan interval and time zone handling

## Features

- Current:
  - Phase, Azimuth, Elevation, Illumination (%), Distance (km), Parallax (°)
  - Ecliptic longitude/latitude (topocentric and geocentric)
- Next events:
  - Moonrise, Moonset
  - Apogee, Perigee
  - New Moon, First Quarter, Full Moon, Last Quarter
- Lunation context:
  - Geocentric ecliptic longitude/latitude at next New/Full Moon (true-of-date)
  - Zodiac sign and degree at next New/Full Moon
- Binary sensor:
  - Moon above horizon (on/off)

## Installation

Manual install:
1. Download or clone the repository.
2. Copy the folder `custom_components/moon_astro` to your Home Assistant config directory:
   - Final path: `<config>/custom_components/moon_astro`
3. Restart Home Assistant.
4. Go to Settings → Devices & Services → Add Integration → search for “Moon Astro”.

HACS (recommended):
1. In HACS → Integrations → Custom repositories, add `https://github.com/svalsemey/hassio-moon-astro/` as type “Integration”.
2. Install “Moon Astro”.
3. Restart Home Assistant and add the integration from Devices & Services.

## Configuration

During setup, the integration pre-fills your Home Assistant latitude, longitude, and elevation. You can adjust them if needed.

Options (via Configure on the integration):
- Scan interval (seconds). Default: 300 (min 30, max 21600)
- Use Home Assistant time zone. Default: true

## Entities

Sensors (examples of entity names; actual names are localized):
- Phase
- Azimuth (°)
- Elevation (°)
- Illumination (%)
- Distance (km)
- Parallax (°)
- Ecliptic longitude/latitude (topocentric)
- Ecliptic longitude/latitude (geocentric)
- Next moonrise / moonset (timestamp)
- Next apogee / perigee (timestamp)
- Next new moon / first quarter / full moon / last quarter (timestamp)
- Ecliptic lon/lat at next new moon
- Ecliptic lon/lat at next full moon
- Zodiac sign at next new/full moon
- Zodiac degree at next new/full moon

Binary sensors:
- Moon above horizon

## Localization

Translations included for multiple languages (e.g., en, sv, tr, uk, zh_CN, sl, etc.). Entity names and certain state texts are localized via Home Assistant’s translation system.

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
- If you modify translation files, reload the integration to apply changes.
- Ensure the key name for “ecliptic_latitude_next_new_moon” is consistently spelled across code and translations.

## Contributing

Issues and PRs are welcome!
- Issues: https://github.com/svalsemey/hassio-moon-astro/issues
- Pull requests: please lint and follow HA integration best practices.

## License

MIT © Sébastien VALSEMEY

## Credits

- Skyfield by Brandon Rhodes
- Home Assistant core team and community
