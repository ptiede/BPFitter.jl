# Name fragments of files that are not raw source data: time-averaged copies and
# the output of an earlier correction.
const DERIVED_TAGS = ("+avg", "_corr", "_bpcal")

"""
    track_source_paths(data_dir; skip = String[]) -> Vector{String}

The raw per-source uvfits files under `data_dir`, in name order, leaving out
the sources named in `skip` (files ending in `_<source>.uvfits`) and files
whose names mark them as derived (`+avg`, `_corr`, `_bpcal`).
"""
function track_source_paths(data_dir; skip = String[])
    isdir(data_dir) || error("no such directory: $data_dir")
    paths = sort(
        filter(readdir(data_dir; join = true)) do p
            endswith(p, ".uvfits") || return false
            name = first(splitext(basename(p)))
            any(t -> occursin(t, name), DERIVED_TAGS) && return false
            return !any(s -> endswith(name, "_" * s), skip)
        end
    )
    isempty(paths) && error("no raw source uvfits under $data_dir")
    return paths
end

"""
    load_track(data_dir; skip = String[]) -> ProcessingSet

Every raw source file under `data_dir` in one processing set, with its
visibilities recorded as `uncalib` and each scan named `<source>_<scan>`.

Time samples are absolute epochs, so files written against different reference
dates share one axis. Each file numbers its scans from 1, and a track's scan
names must be unique, hence the source prefix. HOPS writes `BUNIT = 'JY'` on
visibilities that are still correlation coefficients; recording them as
`uncalib` is what lets `apriori_calibrate!` put them in janskys.

The files are read one at a time: reading on several threads is not safe in the
FITSFiles extension.
"""
function load_track(data_dir; skip = String[])
    paths = track_source_paths(data_dir; skip)
    sets = map(paths) do p
        ps = load_uvfits(p)
        for ms in values(ps)
            DimensionalData.metadata(ms[:visibility])[:units] = "uncalib"
            src = source_name(ms)
            ms[:scan_name] = map(n -> "$(src)_$(n)", ms[:scan_name])
        end
        ps
    end
    ps = merge(sets...)
    @info "sources joined" files = basename.(paths) measurement_sets = length(ps)
    return ps
end

"The source a Measurement Set observes; a HOPS file holds one."
source_name(ms) = String(only(unique(collect(ms[:field_name]))))

track_sources(ps) = sort!(unique(source_name(ms) for ms in values(ps)))

"""
    channel_ifs(ps) -> Dict{String, UnitRange{Int}}

The ANTAB IF of each channel of every spectral window. HOPS writes one channel
per IF, which `load_uvfits` reads as the channels of one window, so channel `k`
is IF `k`.
"""
channel_ifs(ps) = Dict(XRadio.spectralwindow(ms) => 1:length(XRadio.frequencies(ms)) for ms in values(ps))

"""
    usable_weights(ms) -> Array

The weights of `ms` with every flagged sample zeroed, laid out as its
visibilities. A flag leaves the weight at its magnitude, so anything that
accumulates by weight has to consult the flag as well.
"""
usable_weights(ms) = ifelse.(parent(XRadio.flags(ms)), zero(eltype(ms[:weight])), parent(ms[:weight]))

"""
    flag_channels!(ps, station, channels; feeds = nothing) -> ps

Flag `channels` (1-based, in frequency order within each spectral window) on
every product that reads one of `station`'s `feeds`, feed indices as
`Gustavo.feed_pairs` numbers them, or all of its feeds when `feeds` is
`nothing`. Visibilities and weights are left as they are; the flag is what the
solve, every weighted accumulation and the exported uvfits read.
"""
function flag_channels!(ps, station::AbstractString, channels; feeds = nothing)
    for ms in values(ps)
        nchan = length(XRadio.frequencies(ms))
        all(c -> 1 <= c <= nchan, channels) ||
            error("channels $channels fall outside the 1:$nchan channels of $(XRadio.spectralwindow(ms))")
        fp = feed_pairs(ms)
        F = XRadio.flags(ms)
        for (bi, (a, b)) in pairs(XRadio.baselines(ms))
            a == b && continue
            for p in axes(fp, 1)
                fa, fb = fp[p, bi]
                hit = (String(a) == station && (isnothing(feeds) || fa in feeds)) ||
                    (String(b) == station && (isnothing(feeds) || fb in feeds))
                hit && (view(F, Polarization(p), BaselineID(bi), Frequency(channels)) .= true)
            end
        end
    end
    return ps
end

"""
    flag_scan!(ps, stations, scan) -> ps

Flag every channel and product of `scan` on the baselines to `stations`: all of
a station's baselines when it names one, the baseline between them when it
names two.
"""
function flag_scan!(ps, stations, scan::AbstractString)
    found = false
    for ms in values(ps)
        rows = findall(==(scan), collect(ms[:scan_name]))
        isempty(rows) && continue
        found = true
        F = XRadio.flags(ms)
        for (bi, (a, b)) in pairs(XRadio.baselines(ms))
            a == b && continue
            ends = (String(a), String(b))
            hit = length(stations) == 1 ? only(stations) in ends : issetequal(stations, ends)
            hit && (view(F, Ti(rows), BaselineID(bi)) .= true)
        end
    end
    found || error("no scan named $scan")
    return ps
end

"""
    Track

One track loaded and prepared for the solve: the processing set `ps` with its
a-priori amplitudes and flags applied, its configuration, the epoch every hour
is measured from, and its scan timeline.
"""
struct Track{C, P, T, R}
    cfg::C
    ps::P
    origin::T
    rows::R
end

"""
    prepare_track(cfg::TrackConfig) -> Track

Load the track's raw uvfits, put them in janskys from the ANTAB file (the
track's own, or the simulated one), and apply the configured flags.
"""
function prepare_track(cfg::TrackConfig)
    ps = load_track(cfg.data; skip = cfg.skip_sources)
    origin = track_origin(ps)
    rows = scan_rows(ps, origin)
    overlap = overlapping_sources(rows)
    isempty(overlap) || error(
        "track $(cfg.name) has sources observed at the same time by two sub-arrays: ",
        join(("$a / $b" for (a, b) in overlap), ", "), "; list one of each pair in `skip_sources`"
    )
    antab = antab_path(cfg)
    isfile(antab) || error(
        isnothing(cfg.antab) ?
            "track $(cfg.name) has no ANTAB file and no simulated one at $antab; run `bpfit sefd` first" :
            "no ANTAB file at $antab"
    )
    XRadio.read_antab!(ps, antab; ifs = channel_ifs(ps))
    foreach(ms -> apriori_calibrate!(ms; cfg.quantization_efficiency), values(ps))
    stations = Gustavo.Calibration.DataGeometry(ps).antennas
    for (st, feeds) in cfg.feed_flags, (feed, channels) in feeds
        st in stations && flag_channels!(ps, st, channels; feeds = (feed,))
    end
    for (st, channels) in cfg.station_flags
        st in stations && flag_channels!(ps, st, channels)
    end
    for (sts, scan) in cfg.scan_flags
        flag_scan!(ps, sts, scan)
    end
    return Track(cfg, ps, origin, rows)
end

"Whether a Measurement Set's scan enters the bandpass solve."
in_solve(cfg::TrackConfig) = ms -> source_name(ms) ∉ cfg.apply_only_sources && only(XRadio.scans(ms)) ∉ cfg.apply_only_scans
