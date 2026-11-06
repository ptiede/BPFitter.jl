# BPFitter.jl

**Note this README was written with the aid of Claude. If there is an issue please file it on GitHub.**

Joint bandpass calibration of HOPS-fringed VLBI tracks, built on
[Gustavo](https://github.com/ptiede/Gustavo.jl).

The goal of this package is to correct residual issues in the EHT-HOPS pipeline,
specifically for EHT data that has both ALMA pre-polconvert and NOEMA, which
often causes issues with the data. The main things we try to correct for are 
residual delays/bandpass at sites, or relative delay breaks between two feeds. 
In addition, we do a very small fringe fit on the data to detect any residual
delays in the data. 

This repo contains various scripts built upon `Gustavo` which provides the basic
functions and interface. Solutions are saved as zarr files. 

## Install

BPFitter needs Julia 1.12 or later thanks to its *app* functionality. 

Install the `bpfit` app with
```
pkg> app add https://github.com/ptiede/BPFitter.jl
```

Julia puts the app in `~/.julia/bin`, which must be on the `PATH`. Update it
with `pkg> app update BPFitter`.

## Functionality

This repo provides functionality through several subfunctions in the `bpfit` function.
To run the function, specify a configuration TOML file. An example is provided in `examples/campaign.toml`



The current functions are:
```
bpfit inspect campaign.toml            # sources, stations, timeline, dead channels, sub-array overlaps, odd SEFDs
bpfit sefd    campaign.toml            # simulate ANTAB files for tracks that have none
bpfit breaks  campaign.toml            # suggest bandpass breaks; writes a table and plots for review
bpfit solve   campaign.toml            # solve; writes solutions, χ² tables and plots
bpfit apply   campaign.toml            # write corrected uvfits from the saved solutions
bpfit run     campaign.toml [track…]   # sefd where needed, solve, apply, per track
```
Each command runs on the named tracks, or on all of them, one after another in
one Julia process: the code is compiled for the first track only, which on a
small track is most of the time a command takes. Run a single track per call
to keep a failure on one track from stopping the others.


Some additional notes are below

## Apriori Calibration

There are two ways to provide a priori amplitude information i.e., SEFDs, which include gain dpfu's, system temperatures, etc.
The first preferred way is to provide a standard EHT metadata package. We assume you the user have access to this. 

The second way if the metadata package is unavailable is to simulate it with [ngehtsim](https://github.com/Smithsonian/ngehtsim) (its `refactor` branch),
which BPFitter runs with [uv](https://docs.astral.sh/uv/). This assumes that a user has `uv` installed and in its `PATH`.

### File layout

We assume that tracks are organized into distinct folders. That is, all the sources in a
single EHT track will be assumed to be in a single folder. This is important because we 
assume a priori (i.e., unless breaks are put in) that the bandpass is constant across all
sources. We do this to transfer bandpass solutions from strong calibration targets to weak ones. 
If you really want an independent bandpass per scan, you can manually put breaks in between each scan.


In more detail, `bpfit` reads one folder of HOPS uvfits per track and one ANTAB file per
track, wherever the campaign file's `data` and `antab` patterns point. Paths in
the campaign file are relative to it. A layout that matches
[`examples/campaign.toml`](examples/campaign.toml):

```
campaign.toml
data/
  hops/
    track1/                       # data = "data/hops/{track}"
      hops_track1_SRC1.uvfits     # one file per source
      hops_track1_SRC2.uvfits
    track2/
  antab/
    t1.AN                         # antab = "data/antab/{label}.AN"
    t2.AN
```

- A track's folder holds one uvfits per source, we assume that we have not band 
  averaged the data, e.g., there are 32 EHT IFs (or SPWs or HOPS channels depending on your preferred language)
  The source is the file's object name, and its scans are numbered within the file. Files whose names
  contain `+avg`, `_corr` or `_bpcal` are left out, and so are those ending in
  `_<source>.uvfits` for a source in `skip_sources`.
- A track whose ANTAB file is missing takes the SEFDs `bpfit sefd` simulates,
  written to its `output` folder as `simulated.AN`.

Each track's outputs go to its `output` and `corrected` folders as specified in the TOML file:

```
bpfit/track1/                     # output = "bpfit/{track}"
  chi2.txt  chi2_*.png            # fit tables and plots (solve)
  timeline.txt
  phase_<A>-<B>.png  amp_<A>-<B>.png
  solution_phase.png  solution_amp.png  break_delta.txt     # with breaks
  bandpass.zarr  scan_delays.zarr # solutions, read by apply
  break_suggestions/              # breaks: suggestions.txt and a plot per step
  simulated.AN                    # sefd
data/bpcal/track1/                # corrected = "data/bpcal/{track}"
  hops_track1_SRC1_bpcal.uvfits   # apply and run, named by export_name
```

### What `bpfit solve` does

1. Load every source's uvfits for the track, put the visibilities in janskys
   from the ANTAB system temperatures and gain curves, and apply the configured
   flags.
2. Solve one bandpass per station and feed (cut at its breaks) over all
   sources, with an explicit per-scan, per-baseline source term. The phase and
   log-amplitude of each segment are Ornstein–Uhlenbeck processes along
   frequency, with their correlation bandwidth and amplitude fit to the data.
3. Fit, on the bandpass-corrected data, a delay and a rate per station and
   scan and each feed's delay offset from the first feed per scan, apply them,
   and solve again, with per-scan terms for any station configured with
   `scan_terms`. A delay that jumps from one scan to the next, on one feed or
   both, needs no flag or break.
4. Apply the solution and report the goodness of fit against the data without
   it.

### Breaks

`bpfit breaks` starts from the breaks in the campaign file, if any, and:

1. solves the track and, for each station, splits every baseline to it at each
   gap between its scans. The F ratio of a two-sided fit against one mean
   spectrum peaks at a step, and the ratio of the two means is the change
   across it. Each station's local maxima of that F whose change, the median
   over its baselines, reaches `min_shape_deg` in phase or
   `min_scatter_percent` in amplitude are its candidates, ranked by that F.
   A step of another station shows on one baseline only, so its median change
   stays small;
2. refits the track with each of the best `tries` candidates added and keeps
   the one with the largest F = (Δχ²/Δp)/(χ²/ν) if it reaches `min_F` (and the
   1 − α quantile of F(Δp, ν)) and its change in the refit reaches one of the
   sizes. The residuals carry structure from scan to scan, so almost any cut
   passes the quantile alone; `min_F` sets a higher floor. At high
   signal-to-noise a change of a degree passes it too; the sizes keep the
   search to large changes;
3. repeats until no candidate passes, then retests every break against the
   solve with all of them and drops any that no longer pass.

Every solve of the search uses the per-scan delays of one solve with the
configured breaks.

It writes `<output>/break_suggestions/`: `suggestions.txt`, one row per break
or candidate with its F test, the size of the bandpass change across it per
feed and the sources observed around it, and one plot per row of the phase and
amplitude change against frequency. Judge each by the size of its change and
whether it lines up with a change of sources. To adopt breaks, add them to the track's
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
   for a break, a bad channel or a bad scan. Below the totals, `chi2.txt` has
   three tables that point to the usual fixes:
   - **per source**: how well each source's scans fit. The sum of squares
     grows with signal-to-noise, so bright sources stand highest; compare a
     source with others of similar brightness. `chi2_sources.png` shows each
     source's mean corrected residual spectrum. The bandpass is common to all
     sources, so a source with a spectral feature of its own (an absorption
     line along its line of sight, for example), or one that changes within
     its scans, pulls the solved bandpass toward it: it shows
     the feature with one sign and every other source shows it with the
     opposite sign. List such a source in `apply_only_sources`: it is
     corrected with the solved bandpass but does not enter the solve;
   - **scans that fit worst on a station**: a scan far above the station's
     typical scan, with the baseline carrying most of it. If one baseline
     carries it, flag that baseline (`scan_flags` with both stations); if it
     is spread over the station's baselines, the station's data in the scan
     are corrupt: flag the station. If the station's bandpass
     changed during it, `bpfit breaks` may suggest breaks on either side,
     which give the scan a bandpass of its own;
   - **channel edges** where a baseline's phase steps differently from scan
     to scan: when the baselines of one station show the same edges, its
     phase slope restarts at fixed channel blocks; fit per-scan terms over
     those blocks with `scan_terms`. The groups must cover the band from
     channel 1; a group of one channel has no slope to fit its delay to, so
     flag that channel too.

   A baseline with very high signal-to-noise, such as one between two
   stations at the same site, can stay many times above 1 after the
   correction from residuals of a degree or a percent that change from scan
   to scan. `chi2_baselines.png` shows such a baseline standing alone; no
   break or flag removes it, and residuals of that size are within the
   systematic errors usually assigned to VLBI data.

   Make these changes and solve again before looking for breaks: a source or
   scan the model cannot fit raises the χ² every break is tested against.
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
| `chi2.txt` | the goodness of fit of the scan averages, with no bandpass and with the solved one: each scan's spectrum against its own band mean, as a sum of squares per degree of freedom, in total, per station and per feed pair. Near 1 is a spectrum the model explains down to the noise. Then the fit per source and per baseline, the scans that fit worst on each station, and the channel edges where a baseline's phase steps differently from scan to scan (see [Workflow](#workflow), step 4). |
| `chi2_scans.png` | each station's fit (tss/dof, 1 is the noise) in each scan, scans in time order with the sources marked; the scans furthest above their station's median are numbered and listed under the plot, with the baseline that carries most of each |
| `chi2_baselines.png` | every baseline at once, as station-by-station grids of the fit before and after the correction and of how much the bandpass changes it (also a table in `chi2.txt`) |
| `chi2_sources.png` | each source's mean corrected residual spectrum, with its correlation `r` with the other sources'; the title names any source whose residual is the opposite of the others' |
| `chi2_edges.png` | the scan-to-scan spread of the phase step at each channel edge of each baseline; the title names any station whose baselines share edges |
| `timeline.txt` | every scan with its source, time span and stations, and the gaps a break could go in |
| `phase_<A>-<B>.png`, `amp_<A>-<B>.png` | each scan's spectrum on one baseline before and after the correction, for the busiest baselines of `watch_stations`, the four baselines that fit worst after it, and the two the solved bandpass changes most |
| `solution_phase.png`, `solution_amp.png` | the solved bandpass of every station and feed, one trace per segment (with breaks) |
| `break_delta.txt` | how much each station's bandpass moves across each break (with breaks): span, delay and its ramp across the band, shape and amplitude scatter per feed |
| `bandpass.zarr`, `scan_delays.zarr` | the solutions, which `apply` reads |

`breaks` writes `break_suggestions/` there, and `sefd` writes `simulated.AN`.

## Reviewing break suggestions

`suggestions.txt` has one row per break or candidate:

- `status`: `configured` (in the campaign file), `suggested`, `dropped` (kept
  during the search, then removed when retested against the solve with all the
  others), `passed` (it passed its test, but another candidate of the same
  round explained more, often the other station of a baseline that shows the
  step), `small` (it passed the F test, but its change in the refit reaches
  neither size) or `rejected`;
- `gap F`: the step test that proposed it, and `F` and `needs`: the refit F test
  and the value it had to reach;
- per feed, the size of the change across it: `span`, the phase change of the
  highest channel less the lowest, in degrees; `ramp`, the part of it a delay
  across the band makes, which the per-scan delays could take up; `shape`, the
  rms phase change beyond that delay and the band mean, in degrees; and
  `scatter`, the rms amplitude change about its band mean, in percent. A
  candidate is refit and kept on its `shape` or `scatter`, so a change that is
  mostly a delay can have a large `span` and a small `shape`;
- the sources the station observed in the hour before and after.

Each row has a plot of the phase and amplitude change against frequency. Things
to weigh:

- **Size.** A real bandpass change usually has a shape you can see in the
  plot: a slope, a curve, or a jump at the edge of a block of channels.
  Channel-to-channel scatter of a degree or two with no shape is close to the
  noise of the solution.
- **Sources.** A change that comes with a change of the sources observed can
  come from the sources rather than the instrument, for example from source
  structure on long baselines or from noisy delays on weak scans.
- **Pairs of stations.** A change seen only on the baseline between two
  stations can be given to either one, and the search lists both. If the two
  are at the same site, the change may belong to that baseline rather than to
  either station's bandpass.
- **One channel.** A change confined to one channel is often better flagged
  than broken.
- **Breaks close together.** Breaks on both sides of one scan give that scan
  a bandpass of its own, fitted from its baselines alone: the station's
  bandpass changed during or around it. Solve with them and look at the scan
  in `chi2.txt`'s worst scans: if it still fits badly, the bandpass moved
  within the scan, and flagging it is better.

The table ends with a `breaks = { … }` line that adopts every suggestion; copy
the parts you accept into the track's table. The search does not change the
campaign file or the solve.

The search fits the per-scan delays once and refits only the candidates with
large changes, so it takes a few minutes per track.

## Memory

A track's visibilities, weights and flags are held in memory with one
delay-corrected copy during a solve or a search. A track of about 100 scans,
10 stations and 32 channels takes about 2 GB, and its solve or search peaks at
6–10 GB, Julia itself included. Several tracks in one call run one after
another, and the peak stays that of the largest.

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
| `scan_terms` | station => channel groups over which a phase level and a delay are fit per scan beside its bandpass, for a station whose phase slope restarts at fixed channel blocks and changes from scan to scan; the groups cover the band from channel 1, and a one-channel group's channel is best flagged |
| `quantization_efficiency` | a number, or station => number, for the a-priori amplitude calibration (default 1) |
| `priors` | `{ bandwidth, phase_sigma, logamp_sigma }`: the centers of the hyperpriors on the bandpass's correlation bandwidth along frequency (Hz) and on the size of its phase (rad) and log-amplitude variations |
| `watch_stations` | stations whose baselines get before/after plots |

`[breaks]`, the settings of `bpfit breaks`: `min_shape_deg` (2) and
`min_scatter_percent` (5), the sizes of the change of which a candidate must
reach one to be refit and a suggestion one to be kept; `min_F` (40) and
`alpha` (10⁻³), the refit F a suggestion must reach and the false-alarm rate it
must also pass; `tries` (4), the candidates refit per round;
`min_gap_minutes` (3), the shortest gap between a station's scans a break may
go in.

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
