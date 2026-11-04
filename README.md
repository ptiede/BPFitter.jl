# BPFitter.jl

Joint bandpass calibration of HOPS-fringed VLBI tracks, built on
[Gustavo](https://github.com/ptiede/Gustavo.jl).

The station bandpass is a property of the instrument, not of the source, so
BPFitter solves one bandpass per station and feed from every source of a track
together. A station's bandpass is cut into segments only where the data show it
changed abruptly. The breaks are set in the campaign file. `bpfit breaks`
suggests candidates: it tests each gap between a station's scans for a step,
refits the track with the strongest candidates, and writes a table and plots
of each for you to review before adding it. Slow changes are left to the
per-scan terms.

## Install

BPFitter needs Julia 1.12 or later. Install the `bpfit` app with

```
pkg> app add https://github.com/ptiede/BPFitter.jl
```

or call it from a clone with `julia --project=path/to/BPFitter -m BPFitter <command> …`.

Tracks without an ANTAB file take simulated SEFDs from
[ngehtsim](https://github.com/Smithsonian/ngehtsim) (its `refactor` branch),
which BPFitter runs with [uv](https://docs.astral.sh/uv/); `uv` must be on the
`PATH`.

## Use

Describe the campaign in a TOML file (see
[`examples/campaign.toml`](examples/campaign.toml)), then

```
bpfit inspect campaign.toml            # sources, stations, timeline, dead channels, sub-array overlaps, odd SEFDs
bpfit sefd    campaign.toml            # simulate ANTAB files for tracks that have none
bpfit breaks  campaign.toml            # suggest bandpass breaks; writes a table and plots for review
bpfit solve   campaign.toml            # solve; writes solutions, χ² tables and plots
bpfit apply   campaign.toml            # write corrected uvfits from the saved solutions
bpfit run     campaign.toml [track…]   # sefd where needed, solve, apply, per track
```

Each command runs on the named tracks, or on all of them.

### What a solve does

1. Load every source's uvfits for the track, put the visibilities in janskys
   from the ANTAB system temperatures and gain curves, and apply the configured
   flags.
2. Solve one bandpass per station and feed (cut at its breaks) over all
   sources, with an explicit per-scan, per-baseline source term. The phase and
   log-amplitude of each segment are Ornstein–Uhlenbeck processes along
   frequency, with their correlation bandwidth and amplitude fit to the data.
3. Fit a delay per station and scan on the bandpass-corrected data, apply it,
   and solve again, with per-scan terms for any station configured with
   `scan_terms`.
4. Apply the solution and report the goodness of fit against the data without
   it.

### Breaks

`bpfit breaks` starts from the breaks in the campaign file, if any, and:

1. solves the track and, for each station, splits every baseline to it at each
   gap between its scans. The F ratio of a two-sided fit against one mean
   spectrum peaks at a step; each station's `peaks` highest local maxima are
   its candidates, ranked by that F;
2. refits the track with each of the best `tries` candidates added, reusing the
   solve's per-scan delays, and keeps the one with the largest
   F = (Δχ²/Δp)/(χ²/ν) if it reaches `min_F` (and the 1 − α quantile of
   F(Δp, ν)). The residuals carry structure from scan to scan, so
   almost any cut passes the quantile alone; `min_F` sets a higher floor;
3. repeats until no candidate passes, then retests every break against the
   solve with all of them and drops any that no longer pass.

It writes `<output>/break_suggestions/`: `suggestions.txt`, one row per break
or candidate with its F test, the size of the bandpass change across it per
feed and the sources observed around it, and one plot per row of the phase and
amplitude change against frequency. At high signal-to-noise a change of a
degree passes the F test, so judge each by the size of its change and whether
it lines up with a change of sources. To adopt breaks, add them to the track's
`breaks` in the campaign file; the table ends with the entry that adopts all
the suggestions.

## Workflow

A new campaign usually goes through these steps, one track at a time.

1. **Write the campaign file.** Start from
   [`examples/campaign.toml`](examples/campaign.toml): the `data`, `antab` and
   `output` paths, and a `[tracks.<name>]` table per track. The
   [configuration reference](#configuration-reference) lists every key.
2. **Inspect the raw data.** `bpfit inspect campaign.toml <track>` reads the
   uvfits and prints what the campaign file needs:
   - channels that read no signal on a station or feed: flag them with
     `feed_flags` or `station_flags`;
   - sources observed at the same time by two sub-arrays: list one in
     `skip_sources`;
   - scans whose ANTAB SEFD is more than 3× above or below that of the
     station's nearby scans of the same source: often a failed measurement,
     but also a source rising or setting fast or a step in a station's
     calibration, so look before flagging them with `scan_flags`;
   - the scan timeline and whether the ANTAB file is there.
3. **SEFDs.** A track whose `antab` file is missing takes simulated SEFDs:
   `bpfit sefd` writes `<output>/simulated.AN` with ngehtsim (`run` does this
   when needed). These set the amplitude scale and the relative weight of the
   baselines; they hardly change the bandpass shape.
4. **Solve without breaks.** `bpfit solve campaign.toml <track>`, then read
   `chi2.txt` and the plots (see [Outputs](#outputs)). A station whose
   `tss/dof` stays well above 1 after the solve is the first place to look
   for a break, a bad channel or a bad scan.
5. **Look for breaks.** `bpfit breaks campaign.toml <track>` and review
   `<output>/break_suggestions/` (see [Reviewing break
   suggestions](#reviewing-break-suggestions)). Add the breaks you accept to the
   track's `breaks` and solve again.
6. **Apply.** `bpfit apply campaign.toml <track>` writes one corrected uvfits
   per source to `corrected`, from the solution `solve` saved.

Once a campaign file is settled, `bpfit run campaign.toml` repeats steps 3, 4
and 6 for every track.

## Outputs

`solve` writes to the track's `output` folder:

| file | what it holds |
|---|---|
| `chi2.txt` | the goodness of fit of the scan averages, with no bandpass and with the solved one: each scan's spectrum against its own band mean, as a sum of squares per degree of freedom, in total, per station and per feed pair. Near 1 is a spectrum the model explains down to the noise. |
| `timeline.txt` | every scan with its source, time span and stations, and the gaps a break could go in |
| `phase_<A>-<B>.png`, `amp_<A>-<B>.png` | each scan's spectrum on one baseline before and after the correction, for the busiest baselines of `watch_stations` and two others |
| `solution_phase.png`, `solution_amp.png` | the solved bandpass of every station and feed, one trace per segment (with breaks) |
| `break_delta.txt` | how much each station's bandpass moves across each break (with breaks) |
| `bandpass.zarr`, `scan_delays.zarr` | the solutions, which `apply` reads |

`breaks` writes `break_suggestions/` there, and `sefd` writes `simulated.AN`.

## Reviewing break suggestions

`suggestions.txt` has one row per break or candidate:

- `status`: `configured` (in the campaign file), `suggested`, `dropped` (kept
  during the search, then removed when retested against the solve with all the
  others), `passed` (it passed its test, but another candidate of the same
  round explained more, often the other station of a baseline that shows the
  step) or `rejected`;
- `gap F`: the step test that proposed it, and `F` and `needs`: the refit F test
  and the value it had to reach;
- per feed, the size of the change across it: `shape`, the rms phase change
  beyond its delay and band mean, in degrees; `scatter`, the rms amplitude
  change about its band mean, in percent; and the `delay` the phase change
  carries;
- the sources the station observed in the hour before and after.

Each row has a plot of the phase and amplitude change against frequency. Things
to weigh:

- **Size.** On a station with very high signal-to-noise, a change of a degree
  or a percent passes the F test easily. A real bandpass change usually has a
  shape you can see in the plot: a slope, a curve, or a jump at the edge of a
  block of channels. Channel-to-channel scatter of a degree or two with no
  shape is close to the noise of the solution.
- **Sources.** A change that comes with a change of the sources observed can
  come from the sources rather than the instrument, for example from source
  structure on long baselines or from noisy delays on weak scans.
- **Pairs of stations.** A change seen only on the baseline between two
  stations can be given to either one, and the search lists both. If the two
  are at the same site, the change may belong to that baseline rather than to
  either station's bandpass.
- **One channel.** A change confined to one channel is often better flagged
  than broken.

The table ends with a `breaks = { … }` line that adopts every suggestion; copy
the parts you accept into the track's table. The search does not change the
campaign file or the solve.

Each round of the search fits the per-scan delays once and refits every
candidate on the same delay-corrected data, so it takes from a few minutes to
about twenty minutes per track.

## Memory

A track's visibilities, weights and flags are held in memory with one
delay-corrected copy during a solve or a search round. A track of about 100
scans, 10 stations and 32 channels takes about 2 GB, and its solve or search
peaks near 7 GB, Julia itself included.

## Configuration reference

The campaign file has a `[defaults]` table, one `[tracks.<name>]` table per
track, and optional `[breaks]` and `[sefd]` tables. A track's keys replace the
default of the same name whole. Paths are relative to the campaign file and may
contain `{track}` and `{label}`. Unknown keys are an error.

Track keys:

| key | meaning |
|---|---|
| `data` | folder of the track's HOPS uvfits, one file per source |
| `antab` | the track's ANTAB file; when the file is missing, the SEFDs are simulated |
| `output` | where solutions, tables and plots go |
| `corrected` | where `apply` writes the corrected uvfits; leave out to write none |
| `export_name` | name of each corrected file, with `{track}` and `{source}` |
| `label` | the value of `{label}` in paths (default: the track name) |
| `ref_ant` | the station the bandpass gauge is pinned to (default: the station on the most scans among those without a break) |
| `breaks` | station => hours at which its bandpass is cut, measured from UT midnight before the track's first sample |
| `feed_flags` | station => feed => channels to flag on that feed |
| `station_flags` | station => channels to flag on every feed |
| `scan_flags` | list of `{ stations, scan }`: one station flags all its baselines in the scan, two flag the baseline between them |
| `skip_sources` | sources left out of the track |
| `apply_only_sources`, `apply_only_scans` | left out of the solve, but corrected with the solved bandpass and exported |
| `scan_delays` | fit a delay per station and scan between two bandpass solves (default `true`) |
| `scan_terms` | station => channel groups over which a phase level and a delay are fit per scan beside its bandpass, for a station whose phase slope restarts at fixed channel blocks and changes from scan to scan |
| `quantization_efficiency` | a number, or station => number, for the a-priori amplitude calibration (default 1) |
| `priors` | `{ bandwidth, phase_sigma, logamp_sigma }`: the centers of the hyperpriors on the bandpass's correlation bandwidth along frequency (Hz) and on the size of its phase (rad) and log-amplitude variations |
| `watch_stations` | stations whose baselines get before/after plots |

`[breaks]`, the settings of `bpfit breaks`: `min_F` (40) and `alpha` (10⁻³),
the refit F a suggestion must reach and the false-alarm rate it must also
pass; `peaks` (2), the candidates each station proposes per round; `tries` (4),
the candidates refit per round; `min_gap_minutes` (3), the shortest gap between
a station's scans a break may go in.

`[sefd]`, for tracks without ANTAB files: `array`, ngehtsim's array settings;
`stations`, ngehtsim's site name for every station code in the data (there is no
default table, because ngehtsim's two-letter codes are not the same as the
data's); `weather` (`"exact"`).

## Julia interface

```julia
using BPFitter
campaign = read_campaign("campaign.toml")
track = prepare_track(campaign["track1"])
write_suggestions("suggestions", track, find_breaks(track; search = campaign.search))
s = solve_track(track)
write_reports(track, s)
export_track(track.cfg, s)
```
