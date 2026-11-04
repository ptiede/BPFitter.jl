"""
    TrackSolve

One joint bandpass solve of a track: the model and breaks it was solved under,
the reference station, the bandpass solution `sol` (with any per-scan terms),
the per-scan delays (`nothing` without them), the corrected data, and the
goodness of fit of the scans in the solve ([`track_tss`](@ref)).
"""
struct TrackSolve{M, B, S, D, P, C}
    model::M
    breaks::B
    ref_ant::String
    sol::S
    delays::D
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

# The bandpass is applied scan by scan, so only one scan's corrected copy exists
# at a time.
function fit_scan_delays(sol, ps)
    return fit(
        FringeClosure(; model = scan_delay_model()),
        reduce(vcat, values(mapsets(g -> fit(FringeSearch(), calibrate(sol, g)), groupby(ps, ByScan())))),
    )
end

"""
    solve_track(track; breaks = track.cfg.breaks) -> TrackSolve

The joint bandpass solve over every source of `track`:

1. with `scan_delays`, one bandpass per station (cut at its `breaks`) solved on
   the scans in the solve, the gauge pinned to the reference station, and a
   delay per station and scan fitted on the bandpass-corrected data and applied
   to the a-priori data;
2. the bandpass solved, with the configured per-scan terms;
3. the solution applied: in full to the scans in the solve, its shape alone
   (no per-scan terms) to the rest.

`track.ps` is left as it is; the corrected data are a copy.
"""
function solve_track(track::Track; breaks = track.cfg.breaks)
    setup = solve_setup(track, breaks)
    delays = track_delays(track, setup)
    keep = in_solve(track.cfg)
    corrected = isnothing(delays) ? deepcopy(track.ps) : calibrate(delays, track.ps)
    sol = fit(Bandpass(; setup.model, setup.gauge), filter(keep, corrected))
    calibrate!(sol, filter(keep, corrected))
    rest = filter(!keep, corrected)
    isempty(rest) || calibrate!(shape_only(sol), rest)
    return TrackSolve(setup.model, breaks, setup.ref_ant, sol, delays, corrected, track_tss(filter(keep, corrected)))
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

"""
    save_solve(dir, s::TrackSolve)

The bandpass solution and the per-scan delays as solution stores under `dir`,
replacing the stores an earlier solve left there.
"""
function save_solve(dir, s::TrackSolve)
    mkpath(dir)
    bandpass, delays = joinpath(dir, "bandpass.zarr"), joinpath(dir, "scan_delays.zarr")
    rm(bandpass; force = true, recursive = true)
    rm(delays; force = true, recursive = true)
    save_solution(bandpass, s.sol)
    isnothing(s.delays) || save_solution(delays, s.delays)
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
