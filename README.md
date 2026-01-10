# Moon Astro

High-precision Moon ephemeris integration for Home Assistant, powered by Skyfield (DE440). Provides current topocentric/geocentric ecliptic coordinates, illumination, distance, parallax, moonrise/set, lunation timestamps (next and previous), apogee/perigee, and zodiac information.

- Accurate ecliptic-of-date conversion (IAU 1980 nutation, true obliquity)
- Topocentric elevation/azimuth from your configured location
- Fully localized entities and state translations
- Config flow with options for scan interval, time zone handling, and high-precision calculations

## Features

- Current:
  - Phase, Azimuth, Elevation, Illumination (%), Distance (km), Parallax (°)
  - Ecliptic longitude/latitude (topocentric and geocentric)
  - Zodiac sign and zodiac degree (current moon position)

- Next events:
  - Moonrise, Moonset
  - Apogee, Perigee
  - New Moon, First Quarter, Full Moon, Last Quarter
  - Next full moon name and usual alternative names (Gregorian calendar)

- Previous events:
  - Previous Moonrise, Previous Moonset
  - Previous Apogee, Previous Perigee
  - Previous New Moon, Previous First Quarter, Previous Full Moon, Previous Last Quarter
  - Previous full moon name and usual alternative names (Gregorian calendar)

- Lunation context:
  - Geocentric ecliptic longitude/latitude at:
    - next new moon and next full moon (true-of-date)
    - previous new moon and previous full moon (true-of-date)
  - Zodiac sign and zodiac degree from geocentric coordinates at:
    - current moon position
    - next new moon / next full moon
    - previous new moon / previous full moon

- Binary sensor:
  - Moon above horizon (on/off)

## Precision notes (raw values)

Some intermediate computations use non-rounded (“raw”) values internally to avoid boundary artifacts, especially when values are close to zodiac sign cusps (0°, 30°, 60°…). The values exposed as sensor states remain rounded for readability, but zodiac sign/degree calculations can rely on raw longitudes for higher precision.

## High precision mode

Moon Astro can run in a high precision mode designed to reduce timestamp variability for event sensors (lunation phases, apogee/perigee, moonrise/moonset) by using a finer sampling step and a wider refinement bracket.

In this mode, additional refinement is applied to apsides computations (apogee/perigee) to improve the stability of the computed timestamps:
- the search uses a finer coarse sampling step to identify candidate extrema more reliably
- the extremum bracket is refined more aggressively to converge to a stable solution
- after convergence, the resulting instant is validated against neighboring minute instants to select the most extreme minute-aligned solution

**CPU note:** high precision mode requires more computations and can significantly increase CPU usage. It is generally not recommended on low-power hardware (for example Raspberry Pi models with limited resources). If you enable it, prefer using a longer scan interval and monitor system load.

## Update frequency and minimum granularity

To keep computations stable and avoid unnecessary CPU load, Moon Astro updates on a **minute-based** granularity. Update intervals shorter than one minute are not supported and are not useful for astronomical values that are typically consumed at human-scale cadence in Home Assistant.

When configuring the integration, choose an interval that matches your needs:
- For dashboards and general automations: a few minutes is usually enough.
- For time-sensitive automations (moonrise/set, phase boundaries): use a shorter interval while staying at **one minute or above**, and consider the CPU impact (especially in high precision mode).

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

## Options

Options are available via **Configure** on the integration:

- **Sensors update interval (seconds)**
  Default: 300 (min 60, max 21600)

- **Use Home Assistant time zone**
  Default: true
  When enabled, the integration uses Home Assistant’s configured time zone for timestamps and event calculations.

- **High precision mode**
  Default: false
  Reduces timestamp variability but increases CPU usage due to finer sampling and wider refinement.

## Localization

Translations are included for multiple languages. Entity names and some state values (for example zodiac signs and full moon names) are localized via Home Assistant’s translation system.

If you notice that a translation is incomplete or inaccurate, contributions are very welcome.

## Requirements

- Home Assistant 2024.4+ (implies Python 3.12)
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

- CPU usage is high:
  - Disable “High precision mode”.
  - Increase the sensors update interval.
  - Avoid running with a one-minute interval on low-power hardware.

- Data looks stale:
  - Reduce the scan interval (respect CPU usage and the one-minute minimum granularity).
  - Manually reload the integration from the UI.

## Known Notes

- Ecliptic-of-date values use IAU 1980 nutation with true obliquity.
- If you modify local translation files, reload the integration to apply changes.

## Contributing

Issues and PRs are welcome!
- Issues: https://github.com/svalsemey/hassio-moon-astro/issues
- Pull requests: please lint and follow HA integration best practices.

## License

MIT © Sébastien VALSEMEY

## Credits

- Skyfield by Brandon Rhodes
- Home Assistant core team and community
