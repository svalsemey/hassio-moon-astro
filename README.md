[![HACS Default](https://img.shields.io/badge/HACS-Default-blue?style=flat&logo=homeassistantcommunitystore&logoSize=auto)](https://my.home-assistant.io/redirect/hacs_repository/?owner=svalsemey&repository=hassio-moon-astro&category=integration)
[![HACS Passing](https://github.com/svalsemey/hassio-moon-astro/actions/workflows/validate.yml/badge.svg)](https://github.com/svalsemey/hassio-moon-astro/actions/workflows/validate.yml)
[![Total Downloads](https://img.shields.io/github/downloads/svalsemey/hassio-moon-astro/total.svg)](https://github.com/svalsemey/hassio-moon-astro/releases)
[![Latest Release Downloads](https://img.shields.io/github/downloads/svalsemey/hassio-moon-astro/latest/total.svg)](https://github.com/svalsemey/hassio-moon-astro/releases/latest)
[![Buy Me a Coffee](https://img.shields.io/badge/Buy%20Me%20a%20Coffee-orange?logo=buy-me-a-coffee)](https://www.buymeacoffee.com/sebastien_valsemey)

# Moon Astro

High-precision Moon ephemeris integration for Home Assistant, powered by Skyfield (DE440). Provides current topocentric/geocentric ecliptic coordinates, illumination, distance, parallax, moonrise/set, lunation timestamps (next and previous), apogee/perigee, and zodiac information.

- Coordinates in the true ecliptic and equinox of date (IAU 2006/2000A precession-nutation)
- Topocentric elevation/azimuth from your configured location
- Fully localized entities and state translations
- Config flow with options for scan interval, time zone handling and event refresh resilience

## Features

- Current:
  - Phase, Azimuth, Elevation, Illumination (%), Geocentric distance (km), Parallax (°)
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

- Diagnostics:
  - Downloadable report with the refresh state, the computation times and the latest values, the location being rounded to about 10 km (see [Diagnostics](#diagnostics))

## Sensor update model

Moon Astro uses two coordinators:

- **Periodic sensors** (position, illumination, distance, parallax, ecliptic coordinates, current zodiac position) are recomputed on the configured scan interval.
- **Event-based sensors** (moonrise/moonset, lunation phases, apogee/perigee, full moon names, and the ecliptic/zodiac position at lunations) are recomputed shortly after the earliest upcoming event, and at least once per **events refresh fallback interval**.

A full refresh takes a fraction of a second in a background thread, so both coordinators are refreshed during setup and every entity starts with data.

### Scheduled refresh around the next event

After each refresh, the events coordinator schedules the next one two minutes after the earliest of the next moonrise, moonset, new moon, first quarter, full moon, last quarter, apogee and perigee. The margin keeps the computation safely past the boundary.

### `next_update` attribute for event-based sensors

All **event-based sensors** expose an extra state attribute named `next_update`: the UTC instant of the next event-based recomputation, i.e. the earlier of the scheduled boundary refresh and the fallback refresh.

### Events refresh fallback interval

The **Events refresh fallback interval (seconds)** forces a periodic recomputation of the event-based sensors even when no boundary refresh fires (for example after a time change). The default of 24 hours is plenty.

## Precision notes (raw values)

Some intermediate computations use non-rounded (“raw”) values internally to avoid boundary artifacts, especially when values are close to zodiac sign cusps (0°, 30°, 60°…). The values exposed as sensor states remain rounded for readability, but zodiac sign/degree calculations can rely on raw longitudes for higher precision.

## Accuracy

- Phase instants (new moon, quarters, full moon) come from Skyfield's discrete event search, accurate to the millisecond before the instants are rounded to the minute.
- Moonrise and moonset are found at any latitude, poles included, however briefly the Moon stays above or below the horizon: the culminations (altitude extremes) are located first, then the single crossing of the horizon between two consecutive culminations is bisected to the second, then interpolated linearly to the millisecond before being rounded to the minute.
- Apogee and perigee are the extremes of the geometric geocentric distance, located with Skyfield's vectorized extremum search to within ten seconds; the published minute is the one whose boundary is closest to the extreme, so it does not depend on the sampling grid.
- Ecliptic coordinates are measured in the true ecliptic and equinox of date: Skyfield rotates the ICRS positions with its IAU 2006/2000A precession-nutation model and the true obliquity. Zodiac signs and degrees derive from the unrounded geocentric longitude.
- Rise, set and the above-horizon state use the Moon center at −34′ of altitude (standard refraction).
- The phase sensor reports new moon, first quarter, full moon or last quarter within about six hours (3° of phase angle) of the exact instant.

## Update frequency and minimum granularity

To keep computations stable and avoid unnecessary CPU load, Moon Astro updates sensor timestamps values on a **minute-based** granularity. Update intervals shorter than one minute are not supported and are not useful for astronomical values that are typically consumed at human-scale cadence in Home Assistant.

When configuring the integration, choose an interval that matches your needs:
- For dashboards and general automations: a few minutes is usually enough.
- For time-sensitive automations: use a shorter interval while staying at **one minute or above**, and consider the CPU impact.

## Installation

### HACS (recommended)

Use this link to directly go to the repository in HACS

[![Add this integration to Home Assistant](https://my.home-assistant.io/badges/hacs_repository.svg)](https://my.home-assistant.io/redirect/hacs_repository/?owner=svalsemey&repository=hassio-moon-astro&category=integration)

or

1. In HACS → Integrations → Custom repositories, add `https://github.com/svalsemey/hassio-moon-astro/` as type “Integration”.
2. Install “Moon Astro”.
3. Restart Home Assistant and add the integration from Devices & Services.

### Manual install
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

  Timestamps are always published as UTC instants and rendered by Home Assistant in your locale. The time zone only defines the calendar month used to name full moons and to detect blue moons.

- **Time zone**

  Used only if "Use Home Assistant time zone" is disabled

- **Events refresh fallback interval (seconds)**

  Default: 86400 (24 hours)

  Forces a periodic refresh of event-based sensors if the scheduled refresh around the next event is missed (restart, time change, delay). Lower values refresh more often; higher values reduce CPU usage.

## Localization

Translations are included for multiple languages. Entity names and some state values (for example zodiac signs and full moon names) are localized via Home Assistant’s translation system.

If you notice that a translation is incomplete or inaccurate, contributions are very welcome.

If you modify local translation files, reload the integration to apply changes.

## Requirements

- Home Assistant 2026.9+ (Python 3.14)
- Python dependencies are installed automatically by Home Assistant:
  - `skyfield==1.55`

Note: Moon Astro downloads the JPL DE440 ephemeris (about 115 MB) from ssd.jpl.nasa.gov into <config>/.skyfield during the initial setup. The file is validated at each start-up; a damaged copy is downloaded again automatically.

## Privacy and Network

- No external network calls except the one-time download of the DE440 ephemeris from ssd.jpl.nasa.gov (JPL).
- Computations are local on your Home Assistant instance.
- The diagnostics file is only generated when you download it and only leaves your machine if you share it; it holds a location rounded to about 10 km (see [Privacy of the diagnostics](#privacy-of-the-diagnostics)).

## Troubleshooting

- Entities show “unavailable”:
  - Check Logs → “Moon Astro” for errors.
  - Ensure your HA location and elevation are set.
  - Verify ephemeris download completed (see `<config>/.skyfield` content).
  - If the download fails at start-up, the setup is retried automatically; check that Home Assistant can reach ssd.jpl.nasa.gov.

- Time zone issues:
  - Toggle the “Use Home Assistant time zone” option or verify your HA system time zone.

- CPU usage is high:
  - Check the computation times in the [diagnostics](#computation-times).
  - Increase the sensors update interval; keep the events refresh fallback interval at its 24-hour default.

- Data looks stale:
  - Reduce the scan interval (respect CPU usage and the one-minute minimum granularity).
  - If only event-based sensors look stale, reduce the events refresh fallback interval.
  - Manually reload the integration from the UI.

## Diagnostics

Moon Astro supports the Home Assistant diagnostics: a JSON file gathering what is needed to analyze a problem, namely the configuration, the refresh state, the computation times and the latest computed values.

### Downloading the diagnostics

1. Go to **Settings → Devices & services → Moon Astro**.
2. Open the **⋮** menu of the integration entry and select **Download diagnostics**.
3. Home Assistant saves a `.json` file; attach it unchanged to your issue (drag and drop).

The diagnostics are available once the integration is loaded.

### Content of the file

Home Assistant adds its own sections at the top of the file:

| Section | Content |
| --- | --- |
| `home_assistant` | Home Assistant version, installation type, operating system, architecture, Python version, time zone |
| `custom_components` | Installed custom integrations and their versions |
| `integration_manifest` | Moon Astro manifest, including its version and requirements |
| `setup_times` | Time spent setting up the integration, which includes the validation of the ephemeris and the first two computations |

The Moon Astro content sits under `data`:

| Key | Content |
| --- | --- |
| `entry.data` | Observer latitude, longitude and elevation, rounded (see [Privacy of the diagnostics](#privacy-of-the-diagnostics)) |
| `entry.options` | Options as stored; an empty object means that the defaults apply |
| `coordinators.current` | Computation of the periodic sensors (position, illumination, distance, parallax, ecliptic coordinates, current zodiac position) |
| `coordinators.events` | Computation of the event-based sensors (moonrise/moonset, lunation phases, apogee/perigee, full moon names, lunation coordinates) |

Each coordinator section holds:

| Key | Content |
| --- | --- |
| `last_update_success` | `false` when the latest refresh failed; the entities of this coordinator are then unavailable |
| `last_exception` | Error raised by the latest failed refresh, or `null` |
| `update_interval_seconds` | Effective periodic interval: the sensors update interval for `current`, the events refresh fallback interval for `events` |
| `computation` | Statistics of the computation times (see below) |
| `data` | Latest values published to the entities, as raw state codes and UTC instants |
| `next_refresh` | `events` only: UTC instant of the next recomputation, as exposed by the `next_update` attribute |

### Computation times

The `computation` object of each coordinator holds:

| Key | Content |
| --- | --- |
| `reference_time` | UTC minute the latest values were computed for; together with the rounded location, it allows the computation to be reproduced |
| `count` | Number of successful computations |
| `last_ms` | Duration of the latest computation, in milliseconds |
| `mean_ms` | Mean duration, in milliseconds |
| `min_ms` | Shortest duration, in milliseconds |
| `max_ms` | Longest duration, in milliseconds |

The durations are measured around the Skyfield computation itself, in the background thread running it: they exclude the wait for a free worker thread but include any slowdown caused by other activity on the machine. The statistics start over each time the integration is loaded (Home Assistant start, reload, options or location change).

As an indication, on an ARM64 development machine a computation takes about 3 ms for `current`, and for `events` about 50 ms at mid-latitudes and up to about 200 ms in polar regions, where the Moon can stay above or below the horizon for days and the moonrise/moonset search window widens. A small single-board computer can be several times slower.

Reading the statistics:

- The CPU time spent per day is about `mean_ms` multiplied by the number of computations per day: with the default 300-second interval, `current` computes 288 times a day, i.e. about one second of CPU per day at 3 ms.
- `events` recomputes shortly after each moonrise, moonset, phase and apsis, and at least once per events refresh fallback interval: a few times a day. A `count` growing much faster points to a scheduling problem worth reporting.
- A `max_ms` far above `mean_ms` usually comes from the first computation, which reads the ephemeris from the disk, or from a temporary load of the machine (start-up, backup, other integrations).

### Privacy of the diagnostics

The integration never sends the diagnostics anywhere: the file is generated when you download it and only leaves your machine if you share it.

Astronomical values depend on the observer location, so a problem cannot be reproduced without it. The location is therefore rounded rather than removed, coarsely enough not to identify a home:

| Value | Rounding | Remaining precision |
| --- | --- | --- |
| Latitude | 0.1° | ±5.5 km |
| Longitude | 0.1° | ±5.5 km at the equator, ±3.9 km at 45° of latitude |
| Elevation | 100 m | ±50 m |
| Topocentric Moon azimuth, elevation and ecliptic longitude/latitude (`coordinators.current.data`) | 0.1° | ±5.5 km at best once converted into a location |

The rounded latitude and longitude designate an area of about 100 km², the size of a town or of a district of a city.

The topocentric angles are rounded too because they depend on the observer location: given the time of the computation, which the geocentric values reveal anyway, the Moon azimuth and elevation published to the ten-thousandth of a degree would locate the observer to within about ten meters, on the principle of celestial navigation. Rounded to 0.1°, they reveal no more than the rounded coordinates.

The other values are either the same everywhere on Earth (phase, illumination, distance, parallax, geocentric coordinates, lunation and apsis instants) or published to the minute (moonrise and moonset), which matches about a quarter of a degree of longitude.

The file holds no exact coordinates, no account and no network information. The sections added by Home Assistant describe the installation (versions, operating system, time zone, user running Home Assistant, other custom integrations): review the file before posting it publicly.

## Reporting an issue

Questions, help requests and ideas to discuss belong in [Discussions](https://github.com/svalsemey/hassio-moon-astro/discussions). Issues are opened through forms that ask for the information needed to analyze them:

- [Bug report](https://github.com/svalsemey/hassio-moon-astro/issues/new?template=bug_report.yml): attach the [diagnostics file](#downloading-the-diagnostics) and, when relevant, debug logs and reference values (IMCCE, USNO, timeanddate.com…) with the date and time of the comparison.
- [Feature request](https://github.com/svalsemey/hassio-moon-astro/issues/new?template=feature_request.yml): describe the use case and, for a new astronomical value, its definition and a reference source.

To collect debug logs:

1. Go to **Settings → Devices & services → Moon Astro**, open the **⋮** menu of the integration entry and select **Enable debug logging**.
2. Reproduce the problem (a reload of the integration triggers both computations).
3. Open the **⋮** menu again and select **Disable debug logging**: Home Assistant downloads the log file.

## Contributing

Issues and PRs are welcome!
- Issues: https://github.com/svalsemey/hassio-moon-astro/issues/new/choose (see [Reporting an issue](#reporting-an-issue))
- Pull requests: please lint, follow HA integration best practices and keep the [tests](#tests) passing.

## Tests

The test suite runs on every push and pull request (workflow **Tests**) and requires a coverage of at least 95% of the integration code.

### Running the tests

From the repository root, with [uv](https://docs.astral.sh/uv/) and `jq`:

```shell
uv venv --python 3.14
uv pip install --requirement requirements_test.txt $(jq -r '.requirements[]' custom_components/moon_astro/manifest.json)
uv run --no-project pytest --cov
```

The test dependencies are listed in `requirements_test.txt`, the integration requirements are read from its manifest so that both always match.

### Test ephemeris

The tests never download anything. They use `tests/fixtures/de440_2026_2027.bsp`, a 280 KB excerpt of DE440 holding the very coefficients of the full kernel for 2026 and 2027, so that results are identical to those of the 115 MB file. It was produced from the full kernel with:

```shell
python -m jplephem excerpt 2026/1/1 2028/1/1 de440.bsp tests/fixtures/de440_2026_2027.bsp
```

The clock is frozen at 2026-11-02 22:52 UTC, so the tests keep passing after 2027.

### Snapshots

Entity states and the diagnostics are compared with [Syrupy](https://github.com/syrupy-project/syrupy) snapshots stored in `tests/snapshots`. After an intended change of the values or of the entity definitions, regenerate them and review the differences before committing:

```shell
uv run --no-project pytest --snapshot-update
```

### What is tested

| File | Scope |
| --- | --- |
| `test_config_flow.py` | Setup from the Home Assistant location, download failure and recovery, invalid locations, single instance, reconfiguration, options and their limits |
| `test_init.py` | Loading, unloading, setup retry when the ephemeris or a first computation fails, reload on options change, ephemeris deletion on removal |
| `test_coordinator.py` | Recomputation after the next event, fallback refresh, no refresh after unloading, computation time statistics |
| `test_sensor.py`, `test_binary_sensor.py` | Every entity and its state (snapshots), sensors disabled by default, unavailability while a computation fails |
| `test_diagnostics.py` | Diagnostics content (snapshot), rounding of the location and of the topocentric angles |
| `test_utils.py` | Ephemeris cache, download, damaged files, network errors, deletion, time zone resolution |
| `test_astronomy.py` | Accuracy and edge cases of the computations (see below) |

The astronomical computations are checked against independent references and against their definitions:

- Phase instants agree within one minute with those published by the [US Naval Observatory](https://aa.usno.navy.mil/data/MoonPhases), over five lunations.
- Moonrises and moonsets in Paris agree within three minutes with those of the USNO, which refer to the upper limb of the Moon rather than its center.
- At six locations from the equator to both poles, the Moon center crosses the −34′ horizon within the published minute of every rise and set, and a previous and next instant is found for every kind of event.
- The phase angle is within 0.005° of 0°, 90°, 180° or 270° at the published phase instants, and the apogees and perigees are the distance extremes over three days.
- Blue moons, the time zone defining the month of a full moon, phase names near the principal phases, zodiac cusps and the rounding to the minute.

## License

MIT © Sébastien VALSEMEY

## Credits

- Skyfield by Brandon Rhodes
- Home Assistant core team and community
