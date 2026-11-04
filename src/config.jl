"""
    BreakSearch(; min_F = 40, alpha = 1e-3, peaks = 2, min_gap_minutes = 3, tries = 4)

Settings of the automatic break search ([`find_breaks`](@ref)), which looks
for large steps in a station's bandpass and leaves slow changes to the per-scan
terms.

- `min_F`: the refit F a break must reach. The F test treats what a solve
  leaves as noise, but the residuals carry structure from scan to scan, so
  nearly any cut passes its `1 - alpha` quantile, however small `alpha` is
  (with 128 parameters the quantile stays below about 2); this sets a higher
  floor. A break must pass both. At high signal-to-noise a small change passes
  it too, so the size of each change is reported for review.
- `alpha`: the false-alarm rate of the refit F test.
- `peaks`: the candidates each station proposes per round, its highest local
  maxima of the gap F ratio ([`Gustavo.Fring.step_tests`](@ref)) along the
  track. Scan-to-scan scatter keeps the gap F of a real step low on some
  stations, so its value ranks candidates but does not decide them.
- `min_gap_minutes`: the shortest gap between a station's scans a break may be
  placed in.
- `tries`: the candidates refit per round, highest gap F first; the one with
  the largest refit F that passes is kept.
"""
struct BreakSearch{T}
    min_F::T
    alpha::T
    min_gap_minutes::T
    peaks::Int
    tries::Int
end
function BreakSearch(; min_F = 40.0, alpha = 1.0e-3, peaks = 2, min_gap_minutes = 3.0, tries = 4)
    return BreakSearch(promote(min_F, alpha, min_gap_minutes)..., peaks, tries)
end

"""
    BandpassPriors(; bandwidth = 96e6, phase_sigma = 0.1, logamp_sigma = 0.05)

The hyperprior centers of the bandpass shape. Each station's phase (radians) and
log-amplitude (nepers) is an Ornstein–Uhlenbeck process along frequency whose
correlation bandwidth (Hz) and amplitude are fit to the data under log-normal
hyperpriors centered on these values.
"""
struct BandpassPriors{T}
    bandwidth::T
    phase_sigma::T
    logamp_sigma::T
end
function BandpassPriors(; bandwidth = 96.0e6, phase_sigma = 0.1, logamp_sigma = 0.05)
    return BandpassPriors(promote(bandwidth, phase_sigma, logamp_sigma)...)
end

"""
    SimulatedSEFD(; array, stations, weather = "exact")

How [`simulate_antab`](@ref) runs ngehtsim for a track with no ANTAB file.
`array` is the ngehtsim array whose receiver and efficiency settings apply.
`stations` maps each station code in the data to its ngehtsim site name
(`"KT" => "KP"`); ngehtsim's own two-letter codes are not the EHT's, so every
station needs an entry. `weather` is ngehtsim's weather form: `"exact"` takes
each UT day's record, `"median"` and the others take climatology for the month.
"""
struct SimulatedSEFD
    array::String
    stations::Dict{String, String}
    weather::String
end
function SimulatedSEFD(; array, stations, weather = "exact")
    return SimulatedSEFD(String(array), Dict(String(k) => String(v) for (k, v) in stations), String(weather))
end

"""
    TrackConfig

The settings of one track's reduction, from a campaign file
([`read_campaign`](@ref)).

Paths:
- `data`: the folder of the track's HOPS uvfits, one file per source.
- `antab`: its ANTAB file, or `nothing` when the track has none and takes the
  SEFDs [`simulate_antab`](@ref) writes to `simulated_antab(cfg)`.
- `output`: where solutions, tables and plots go.
- `corrected`: where the corrected uvfits go, named by `export_name` with
  `{track}` and `{source}` filled in; `nothing` writes none.

Reduction:
- `ref_ant`: the station the bandpass gauge is pinned to, or `nothing` for the
  station on the most scans among those without a break.
- `breaks`: station => hours at which its bandpass is cut, measured from UT
  midnight before the track's first sample. A station not named holds one
  bandpass over the track.
- `feed_flags`: station => feed => channels to flag on that feed.
- `station_flags`: station => channels to flag on every feed.
- `scan_flags`: `(stations, scan)` pairs; one station flags all its baselines
  in the scan, two flag the baseline between them.
- `skip_sources`: sources left out of the track entirely.
- `apply_only_sources`, `apply_only_scans`: left out of the solve, but
  corrected with the solved bandpass and exported.
- `scan_delays`: fit a delay per station and scan after a first solve and solve
  again on the delay-corrected data.
- `scan_terms`: station => channel groups over which the final solve fits a
  phase level and a delay per scan beside its bandpass.
- `quantization_efficiency`: a number, or station => number, for
  `Gustavo.apriori_calibrate!`.
- `priors`: [`BandpassPriors`](@ref).
- `watch_stations`: stations whose baselines get before/after plots.
"""
struct TrackConfig{B, Q, P}
    name::String
    data::String
    antab::Union{Nothing, String}
    output::String
    corrected::Union{Nothing, String}
    export_name::String
    ref_ant::Union{Nothing, String}
    breaks::B
    feed_flags::Dict{String, Dict{Int, Vector{Int}}}
    station_flags::Dict{String, Vector{Int}}
    scan_flags::Vector{Tuple{Vector{String}, String}}
    skip_sources::Vector{String}
    apply_only_sources::Vector{String}
    apply_only_scans::Vector{String}
    scan_delays::Bool
    scan_terms::Dict{String, Vector{UnitRange{Int}}}
    quantization_efficiency::Q
    priors::P
    watch_stations::Vector{String}
end

"The ANTAB file [`simulate_antab`](@ref) writes for a track without one."
simulated_antab(cfg::TrackConfig) = joinpath(cfg.output, "simulated.AN")

"The ANTAB file a track's a-priori calibration reads."
antab_path(cfg::TrackConfig) = something(cfg.antab, simulated_antab(cfg))

"""
    Campaign

Every track of a campaign file, with the break-search and SEFD-simulation
settings they share. Index it by track name: `campaign["track1"]`.
"""
struct Campaign{S}
    path::String
    tracks::Dict{String, TrackConfig}
    order::Vector{String}
    search::S
    sefd::Union{Nothing, SimulatedSEFD}
end
Base.getindex(c::Campaign, name::AbstractString) = haskey(c.tracks, name) ? c.tracks[name] :
    error("no track $(repr(name)) in $(c.path); its tracks are $(c.order)")
track_names(c::Campaign) = c.order

const TRACK_KEYS = (
    "label", "data", "antab", "output", "corrected", "export_name", "ref_ant", "breaks",
    "feed_flags", "station_flags", "scan_flags", "skip_sources", "apply_only_sources",
    "apply_only_scans", "scan_delays", "scan_terms", "quantization_efficiency", "priors",
    "watch_stations",
)
const TOP_KEYS = ("defaults", "tracks", "breaks", "sefd")

function _check_keys(tbl, valid, where)
    bad = setdiff(keys(tbl), valid)
    isempty(bad) || error("unknown key(s) $(sort(collect(bad))) in $where; valid keys are $(collect(valid))")
    return tbl
end

"""
    read_campaign(path) -> Campaign

Read a campaign file. The file has a `[defaults]` table, one `[tracks.<name>]`
table per track, and optional `[breaks]` ([`BreakSearch`](@ref)) and `[sefd]`
([`SimulatedSEFD`](@ref)) tables. A track's keys replace the defaults' keys of
the same name whole: a track's `station_flags` is the complete list for that
track. Path values may contain `{track}` and `{label}`, and relative paths are
taken from the campaign file's folder.

A track's `breaks` are the ones the file gives, or none; `bpfit breaks`
([`find_breaks`](@ref)) suggests breaks to add.

See `examples/` for complete files.
"""
function read_campaign(path)
    tbl = _check_keys(TOML.parsefile(path), TOP_KEYS, path)
    haskey(tbl, "tracks") || error("$path has no [tracks] table")
    defaults = _check_keys(get(tbl, "defaults", Dict{String, Any}()), TRACK_KEYS, "[defaults] of $path")
    root = dirname(abspath(path))
    order = sort!(collect(keys(tbl["tracks"])))
    tracks = Dict{String, TrackConfig}()
    for name in order
        t = _check_keys(tbl["tracks"][name], TRACK_KEYS, "[tracks.$name] of $path")
        tracks[name] = _track_config(name, merge(defaults, t), root)
    end
    search = BreakSearch(; _symbol_keys(_check_keys(get(tbl, "breaks", Dict{String, Any}()), ("min_F", "alpha", "peaks", "min_gap_minutes", "tries"), "[breaks] of $path"))...)
    sefd = haskey(tbl, "sefd") ? SimulatedSEFD(; _symbol_keys(_check_keys(tbl["sefd"], ("array", "stations", "weather"), "[sefd] of $path"))...) : nothing
    return Campaign(abspath(path), tracks, order, search, sefd)
end

_symbol_keys(d) = (Symbol(k) => v for (k, v) in d)

function _track_config(name, t, root)
    label = string(get(t, "label", name))
    fill_in(s) = replace(s, "{track}" => name, "{label}" => label)
    place(s) = normpath(joinpath(root, fill_in(s)))
    required(k) = haskey(t, k) ? t[k] : error("track $name has no `$k`; give it in [defaults] or [tracks.$name]")
    output = place(required("output"))
    breaks = _breaks(get(t, "breaks", Dict{String, Any}()))
    priors = BandpassPriors(; _symbol_keys(_check_keys(get(t, "priors", Dict{String, Any}()), ("bandwidth", "phase_sigma", "logamp_sigma"), "priors of track $name"))...)
    q = get(t, "quantization_efficiency", 1.0)
    return TrackConfig(
        name, place(required("data")), haskey(t, "antab") ? place(t["antab"]) : nothing, output,
        haskey(t, "corrected") ? place(t["corrected"]) : nothing,
        fill_in(get(t, "export_name", "hops_{track}_{source}_bpcal.uvfits")),
        get(t, "ref_ant", nothing),
        breaks,
        Dict{String, Dict{Int, Vector{Int}}}(String(st) => Dict{Int, Vector{Int}}(parse(Int, string(f)) => Int.(c) for (f, c) in feeds) for (st, feeds) in get(t, "feed_flags", Dict())),
        Dict{String, Vector{Int}}(String(st) => Int.(c) for (st, c) in get(t, "station_flags", Dict())),
        Tuple{Vector{String}, String}[(String.(f["stations"]), String(f["scan"])) for f in get(t, "scan_flags", [])],
        String.(get(t, "skip_sources", String[])),
        String.(get(t, "apply_only_sources", String[])),
        String.(get(t, "apply_only_scans", String[])),
        get(t, "scan_delays", true),
        Dict{String, Vector{UnitRange{Int}}}(String(st) => [g[1]:g[2] for g in groups] for (st, groups) in get(t, "scan_terms", Dict())),
        q isa AbstractDict ? Dict(String(k) => v for (k, v) in q) : q,
        priors,
        String.(get(t, "watch_stations", String[])),
    )
end

_breaks(tbl) = Dict{String, Vector{Float64}}(String(st) => sort!(Float64.(hrs)) for (st, hrs) in tbl if !isempty(hrs))

"""
    with_breaks(cfg, breaks) -> TrackConfig

`cfg` with its breaks replaced by `breaks` (station => hours).
"""
function with_breaks(cfg::TrackConfig, breaks)
    fields = map(fieldnames(TrackConfig)) do f
        f === :breaks ? _breaks(breaks) : getfield(cfg, f)
    end
    return TrackConfig(fields...)
end
