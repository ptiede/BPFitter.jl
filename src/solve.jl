"""
    TrackSolve

One joint bandpass solve of a track: the model and breaks it was solved under,
the reference station, the bandpass solution `sol` (with any per-scan terms),
the per-scan delays and the per-integration phases and amplitudes (`adhoc`),
each `nothing` without them, the corrected data, and the goodness of fit of the
scans in the solve ([`track_tss`](@ref)).
"""
struct TrackSolve{M, B, S, D, A, P, C}
    model::M
    breaks::B
    ref_ant::String
    sol::S
    delays::D
    adhoc::A
    corrected::P
    chi2::C
end

"""
    track_tss(ps) -> DimStack

The weighted total sum of squares (`Gustavo.Fring.tss`) of the scan averages
of `ps` on parallel feed pairs, per antenna pair and feed pair: each scan's
spectrum against its own band mean, the χ² of a source constant across the
band. The weights are taken as the inverse variances of each visibility
component.
"""
track_tss(ps) = scan_tss(scan_averages(ps))

"[`track_tss`](@ref) of data already averaged to one sample per scan."
function scan_tss(averages)
    all_pairs = tss(averages)
    q = all_pairs[FeedPair = At([fp for fp in lookup(all_pairs, FeedPair) if fp[1] == fp[2]])]
    total(A) = dropdims(sum(A; dims = Ti); dims = Ti)
    return DimStack((; tss = total(q.tss), dof = total(q.dof)))
end

"The total sum of squares of a [`track_tss`](@ref) result and its degrees of freedom."
scan_chi2(q) = (; chi2 = sum(q.tss), dof = sum(q.dof))

"The free parameters of a solve: every parameter of its bandpass solution."
nparams(s::TrackSolve) = sum(c -> length(c.params), s.sol.components)

"""
    reference_antenna(track, breaks) -> String

The configured reference station, or else the station on the most scans among
those with no break.
"""
function reference_antenna(track::Track, breaks)
    ref = track.cfg.ref_ant
    stations = Gustavo.Calibration.DataGeometry(track.ps).antennas
    if isnothing(ref)
        ref = first(first(filter(p -> !haskey(breaks, first(p)), station_scan_counts(track.rows))))
    end
    ref in stations || error("reference antenna $ref is not on track $(track.cfg.name), which has $stations")
    return ref
end

# Fringe-rate search window (Hz). The data are already fringe-fit, so their
# residual rates are a few mHz; a wider window lets a weak baseline peak on a
# rate alias, such as a sidelobe of a periodic modulation in one station's data.
const SCAN_RATE_WINDOW = (-0.05, 0.05)

# The bandpass is applied scan by scan, so only one scan's corrected copy exists
# at a time. Which feeds of each station the fringe searches detected
# ([`scan_detections`](@ref)) goes into the solution's `info`, which is saved
# with it.
function fit_scan_delays(sol, ps)
    search = FringeSearch(; rate_window = SCAN_RATE_WINDOW)
    searches = reduce(vcat, values(mapsets(g -> fit(search, calibrate(sol, g)), groupby(ps, ByScan()))))
    delays = fit(FringeClosure(; model = scan_delay_model()), searches)
    found = scan_detections(searches)
    info = merge(delays.info, (; detected_station = found.station, detected_scan = found.scan, detected_feeds = found.feeds))
    return CalibrationSolution(delays.geom, delays.name, delays.components, info; spec = delays.spec)
end

# The detection threshold of `FringeClosure`'s default `Stationization`: cells
# above it carry a thousandfold σ in the closure and barely enter the fit.
const DETECTION_PFA = 1.0e-4

"""
    scan_detections(searches) -> NamedTuple

Which feeds of each station have a fringe detection in each scan, from the
fringe searches of [`fit_scan_delays`](@ref): a cell (`valid`,
`pfa ≤ $(DETECTION_PFA)`) on baseline `a`–`b`, feeds `fa`–`fb`, counts for feed
`fa` of `a` and feed `fb` of `b`. One entry per station and scan: `station`,
`scan` and `feeds`, a matrix with one row per entry and one column per feed,
1 where the feed has a detection and 0 where it has none.
"""
function scan_detections(searches)
    station, scan, seen = String[], String[], Vector{Int}[]
    for s in searches
        m = DimensionalData.metadata(s)
        pairs_ = collect(lookup(s, :AntennaPair))
        feeds = collect(lookup(s, :FeedPair))
        names = sort!(unique(reduce(vcat, [[a, b] for (a, b) in pairs_]; init = String[])))
        nfeed = maximum(max(fa, fb) for (fa, fb) in feeds)
        idx = Dict(n => i for (i, n) in pairs(names))
        found = zeros(Int, length(names), nfeed)
        valid, pfa = parent(s[:valid]), parent(s[:pfa])
        for (i, (a, b)) in pairs(pairs_), (j, (fa, fb)) in pairs(feeds)
            valid[i, j] && pfa[i, j] <= DETECTION_PFA || continue
            found[idx[a], fa] = found[idx[b], fb] = 1
        end
        for (k, n) in pairs(names)
            push!(station, n)
            push!(scan, String(m.scan))
            push!(seen, found[k, :])
        end
    end
    nfeed = maximum(length, seen; init = 0)
    return (; station, scan, feeds = [f <= length(x) ? x[f] : 0 for x in seen, f in 1:nfeed])
end

"""
    solve_track(track; breaks = track.cfg.breaks) -> TrackSolve

The joint bandpass solve over every source of `track`:

1. with `scan_delays`, one bandpass per station (cut at its `breaks`) solved on
   the scans in the solve, the gauge pinned to the reference station, and a
   delay per station and scan fitted on the bandpass-corrected data and applied
   to the a-priori data;
2. with `adhoc`, a phase and a log-amplitude per station and integration fitted
   in every scan and applied;
3. the bandpass solved, with the configured per-scan terms;
4. the solution applied: in full to the scans in the solve, its shape alone
   (no per-scan terms) to the rest.

`track.ps` is left as it is; the corrected data are a copy.
"""
function solve_track(track::Track; breaks = track.cfg.breaks)
    setup = solve_setup(track, breaks)
    delays = track_delays(track, setup)
    keep = in_solve(track.cfg)
    corrected = isnothing(delays) ? deepcopy(track.ps) : calibrate(delays, track.ps)
    adhoc = track_adhoc(track, corrected)
    isnothing(adhoc) || calibrate!(adhoc, corrected)
    sol = fit(Bandpass(; setup.model, setup.gauge), filter(keep, corrected))
    calibrate!(sol, filter(keep, corrected))
    rest = filter(!keep, corrected)
    isempty(rest) || calibrate!(shape_only(sol), rest)
    return TrackSolve(setup.model, breaks, setup.ref_ant, sol, delays, adhoc, corrected, track_tss(filter(keep, corrected)))
end

# What a solve with `breaks` fits: the bandpass model with the configured
# per-scan terms (`model`) and without them (`shape`, which fits the per-scan
# delays), the gauge and its reference station.
function solve_setup(track::Track, breaks)
    cfg = track.cfg
    for hrs in values(breaks), h in hrs
        check_gap(track.rows, h)
    end
    all(scan -> any(r -> r.scan == scan, track.rows), cfg.apply_only_scans) ||
        error("apply_only_scans names a scan track $(cfg.name) does not have: $(cfg.apply_only_scans)")
    ref_ant = reference_antenna(track, breaks)
    shape = track_bandpass_model(breaks, track.origin, cfg.priors)
    stations = Gustavo.Calibration.DataGeometry(track.ps).antennas
    terms = Dict(st => g for (st, g) in cfg.scan_terms if st in stations)
    model = isempty(terms) ? shape : with_scan_terms(shape, terms)
    return (; shape, model, gauge = PinAntenna(ref_ant), ref_ant)
end

# The per-scan delays of a solve with `setup` (from `solve_setup`), or
# `nothing` for a track configured without them.
function track_delays(track::Track, setup)
    track.cfg.scan_delays || return nothing
    return fit_scan_delays(fit(Bandpass(; model = setup.shape, setup.gauge), filter(in_solve(track.cfg), track.ps)), track.ps)
end

# The per-integration phases and amplitudes of `ps`, or `nothing` for a track
# configured without them.
track_adhoc(track::Track, ps) = track.cfg.adhoc ? fit(Adhoc(), ps) : nothing

"""
    save_solve(dir, s::TrackSolve)

The bandpass solution, the per-scan delays and the per-integration phases and
amplitudes as solution stores under `dir`, replacing the stores an earlier
solve left there.
"""
function save_solve(dir, s::TrackSolve)
    mkpath(dir)
    bandpass, delays, adhoc = joinpath(dir, "bandpass.zarr"), joinpath(dir, "scan_delays.zarr"), joinpath(dir, "adhoc.zarr")
    foreach(p -> rm(p; force = true, recursive = true), (bandpass, delays, adhoc))
    save_solution(bandpass, s.sol)
    isnothing(s.delays) || save_solution(delays, s.delays)
    isnothing(s.adhoc) || save_solution(adhoc, s.adhoc)
    return dir
end

"""
    export_track(cfg, s::TrackSolve) -> Vector{String}

One corrected uvfits per source under `cfg.corrected`, named by
`cfg.export_name`. Splitting after the joint solve keeps each output the shape
its input had, while the bandpass applied is the one every source constrained.
"""
function export_track(cfg::TrackConfig, s::TrackSolve)
    isnothing(cfg.corrected) && error("track $(cfg.name) sets no `corrected` folder to export to")
    mkpath(cfg.corrected)
    return map(track_sources(s.corrected)) do src
        out = joinpath(cfg.corrected, replace(cfg.export_name, "{source}" => src))
        write_uvfits(out, filter(ms -> source_name(ms) == src, s.corrected); overwrite = true)
        @info "wrote" out
        out
    end
end
