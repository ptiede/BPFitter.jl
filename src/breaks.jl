"""
    BreakTest

A break tested by solving the track with and without it: the station, the hour,
the χ² rise `dchi2` from leaving it out, the parameters it carries `dp`, the
ratio `F = (dchi2/dp)/(χ²/ν)` against the solve with the break, the
`threshold` F had to pass, and whether the change across it in the solve with
the break reaches a size of the search (`large`). `gap_F` is the gap F ratio
that proposed it (`NaN` for a break that was not proposed by the search). The
threshold is the larger of `search.min_F` and the `1 - search.alpha` quantile
of F(dp, ν). A break passes when F exceeds the threshold and it is large.
"""
struct BreakTest{T}
    station::String
    hour::T
    gap_F::T
    dchi2::T
    dp::Int
    F::T
    threshold::T
    large::Bool
end
passed_F(t::BreakTest) = t.F > t.threshold
passed(t::BreakTest) = passed_F(t) && t.large

"""
    break_test(with, without, station, hour, search::BreakSearch; gap_F = NaN, large) -> BreakTest

The F test of a break between a solve that has it (`with`) and one that does
not. The misfit the model with the break leaves scales the χ² rise. `large`
says whether the change across the break in `with` reaches
`search.min_shape_deg` or `search.min_scatter_percent`.
"""
function break_test(with::TrackSolve, without::TrackSolve, station, hour, search::BreakSearch; gap_F = NaN, large)
    χa, ν = scan_chi2(with.chi2)
    χb = scan_chi2(without.chi2).chi2
    dp = nparams(with) - nparams(without)
    F = ((χb - χa) / dp) / (χa / ν)
    threshold = max(search.min_F, quantile(FDist(dp, ν), 1 - search.alpha))
    return BreakTest(String(station), promote(Float64(hour), gap_F, χb - χa)..., dp, F, threshold, large)
end

"""
    BreakCandidate

A place a station's bandpass may step: the station, the hour (inside a gap no
scan of the track straddles), the gap F ratio of the station's baselines split
there, and the change across it measured on the residuals
([`BreakChange`](@ref)), channel by channel the median over the baselines.
"""
struct BreakCandidate{T}
    station::String
    hour::T
    gap_F::T
    change::BreakChange{T}
end

"""
    break_candidates(track, s::TrackSolve, search::BreakSearch) -> Vector{BreakCandidate}

The places each station's bandpass may step in the corrected data of `s`, a
[`search_solve`](@ref) whose corrected data are scan averages, highest gap F
first: the local maxima of the gap F ratio (`Gustavo.Fring.step_tests` on the
scan averages' normalized residuals) over the station's gaps of at least
`search.min_gap_minutes`, leaving out the gaps that hold one of its breaks
already, where the change across the gap reaches `search.min_shape_deg` or
`search.min_scatter_percent`. A step raises the gap F at the gaps around it
too, so it is the local maxima that mark where steps may be. The change is
the median over the station's baselines, so a step that shows on one baseline
only, as another station's step does, stays small. Each candidate is placed at
the middle of the longest stretch of its gap that no scan of the track
straddles; a gap that every instant of which some scan covers yields no
candidate.
"""
function break_candidates(track::Track, s::TrackSolve, search::BreakSearch)
    residuals = normalized_residuals(s.corrected)
    scan_names = parent(residuals.scan_name)
    intervals = scan_intervals(filter(in_solve(track.cfg), track.ps))
    ν = collect(lookup(residuals.r, Frequency))
    ν = ν .- mean(ν)
    free = scan_gaps(track.rows)
    out = BreakCandidate{Float64}[]
    for st in Gustavo.Calibration.DataGeometry(s.corrected).antennas
        rows = filter(r -> isfinite(r.F) && r.gap >= 60 * search.min_gap_minutes, step_tests(residuals, intervals, st))
        span(r) = (track_hours(r.time - r.gap / 2, track.origin), track_hours(r.time + r.gap / 2, track.origin))
        local_max(i) = (i == firstindex(rows) || rows[i - 1].F <= rows[i].F) && (i == lastindex(rows) || rows[i + 1].F <= rows[i].F)
        holds_break(r) = ((lo, hi) = span(r); any(h -> lo <= h <= hi, get(s.breaks, st, Float64[])))
        series = _station_series(residuals, st)
        for r in (rows[i] for i in eachindex(rows) if local_max(i) && !holds_break(rows[i]))
            change = _step_change(series, scan_names, intervals, r.time, ν)
            _large(change, search) || continue
            lo, hi = span(r)
            stretches = [(max(g.t0, lo), min(g.t1, hi)) for g in free if g.t1 > lo && g.t0 < hi]
            isempty(stretches) && continue
            a, b = stretches[argmax([b - a for (a, b) in stretches])]
            push!(out, BreakCandidate(String(st), _round_inside((a + b) / 2, a, b), r.F, change))
        end
    end
    return sort!(out; by = c -> -c.gap_F)
end

_large(c::BreakChange, search) = any(>=(deg2rad(search.min_shape_deg)), c.shape) || any(>=(search.min_scatter_percent / 100), c.scatter)

# The residual spectra (channel × scan) and variances on each baseline to
# `station` and parallel feed pair, with the feed and the sign that turns the
# baseline's phase change into the station's.
function _station_series(residuals, station)
    r, var = residuals.r, residuals.var
    antenna_pairs, feeds = lookup(r, AntennaPair), lookup(r, FeedPair)
    cells = [
        (j, q) for j in axes(r, AntennaPair), q in axes(r, FeedPair)
            if station in antenna_pairs[j] && feeds[q][1] == feeds[q][2]
    ]
    series = map(cells) do (j, q)
        R = parent(permutedims(r[AntennaPair(j), FeedPair(q)], (Frequency, Ti)))
        V = parent(permutedims(var[AntennaPair(j), FeedPair(q)], (Frequency, Ti)))
        scans = [s for s in axes(V, 2) if any(isfinite, view(V, :, s))]
        (; feed = feeds[q][1], sign = first(antenna_pairs[j]) == station ? 1 : -1, R, V, scans)
    end
    return filter!(x -> !isempty(x.scans), series)
end

# The change of a station's bandpass across the time `T` (s) in its `series`:
# on each baseline and feed the ratio of the weighted mean spectra of the scans
# after and before, then channel by channel the median over the baselines.
function _step_change(series, scan_names, intervals, T, ν)
    dp = [Float64[] for _ in eachindex(ν), _ in 1:maximum(x -> x.feed, series)]
    da = [Float64[] for _ in eachindex(ν), _ in 1:maximum(x -> x.feed, series)]
    for x in series
        μb, μa = _split_means(x.R, x.V, x.scans, s -> first(intervals[scan_names[s]]) < T)
        for c in eachindex(ν)
            r = μa[c] / μb[c]
            isfinite(r) || continue
            push!(dp[c, x.feed], x.sign * angle(r))
            push!(da[c, x.feed], log(abs(r)))
        end
    end
    middle(v) = isempty(v) ? NaN : median(v)
    return _station_change(middle.(dp), middle.(da), ν)
end

# The weighted mean of each channel of `R` (channel × scan, variances `V`) over
# the `scans` before the split and over those after; NaN for a channel with
# fewer than two scans on a side.
function _split_means(R, V, scans, before)
    μb = fill!(similar(R, axes(R, 1)), complex(NaN, NaN))
    μa = copy(μb)
    for c in axes(R, 1)
        Wb = Wa = 0.0
        Sb = Sa = zero(eltype(R))
        nb = na = 0
        for s in scans
            isfinite(V[c, s]) || continue
            w = inv(V[c, s])
            if before(s)
                Wb += w
                Sb += w * R[c, s]
                nb += 1
            else
                Wa += w
                Sa += w * R[c, s]
                na += 1
            end
        end
        (nb >= 2 && na >= 2) || continue
        μb[c], μa[c] = Sb / Wb, Sa / Wa
    end
    return μb, μa
end

# The hour to two decimals when that stays inside (a, b), as a person would write it.
_round_inside(h, a, b) = (r = round(h; digits = 2); a < r < b ? r : h)

"""
    scan_averages(ps) -> ProcessingSet

`ps` averaged to one sample per scan. `effective_integration_time` is left out
of the average, which does not know how to combine it.
"""
function scan_averages(ps)
    drop(ms) = haskey(ms, :effective_integration_time) ? delete!(copy(ms), :effective_integration_time) : ms
    return average(XRadio.transform(drop, ps), ByScan())
end

"""
    search_data(track, breaks) -> NamedTuple

The data the break search refits with the per-scan delays of a solve with
`breaks`: the `delays` (`nothing` for a track configured without them), the
scans in the solve with those delays applied (`delayed`), and their scan
averages (`averages`). It is the one copy of the track a round of the search
holds.
"""
function search_data(track::Track, breaks)
    delays = track_delays(track, solve_setup(track, breaks))
    kept = filter(in_solve(track.cfg), track.ps)
    delayed = isnothing(delays) ? kept : calibrate(delays, kept)
    return (; delays, delayed, averages = scan_averages(delayed))
end

"""
    search_solve(track, breaks, data) -> TrackSolve

The bandpass solve with `breaks` on [`search_data`](@ref) `data`, as
[`solve_track`](@ref) solves it after its per-scan delays. The solution is
applied to the scan averages alone, which are its corrected data, and copies
nothing at full resolution. The bandpass and the per-scan terms are constant
over a scan, so its χ² is the one `solve_track` reports with the same delays.
"""
function search_solve(track::Track, breaks, data)
    setup = solve_setup(track, breaks)
    sol = fit(Bandpass(; setup.model, setup.gauge), data.delayed)
    corrected = calibrate(sol, data.averages)
    return TrackSolve(setup.model, breaks, setup.ref_ant, sol, data.delays, corrected, scan_tss(corrected))
end

add_break(breaks, st, h) = merge(breaks, Dict(st => sort!(vcat(get(breaks, st, Float64[]), h))))
function drop_break(breaks, st, h)
    out = Dict(k => (k == st ? filter(!=(h), v) : copy(v)) for (k, v) in breaks)
    return filter!(kv -> !isempty(kv[2]), out)
end

"""
    BreakSearchResult

What [`find_breaks`](@ref) suggests: the breaks (station => hours, the
configured ones included), the final test of each against the solve with all of
them, every candidate refit along the way, the breaks dropped in the final
tests, the final solve (without its corrected data), and the
[`BreakChange`](@ref) of each break or candidate keyed by (station, hour),
measured in the solve that had it.
"""
struct BreakSearchResult{B, T, S, C}
    breaks::B
    final::Vector{T}
    tried::Vector{T}
    dropped::Vector{T}
    solve::S
    changes::Dict{Tuple{String, Float64}, C}
end

"""
    find_breaks(track; search = BreakSearch(), io = stdout) -> BreakSearchResult

Bandpass breaks to suggest for review ([`write_suggestions`](@ref)), by greedy
forward selection starting from the track's configured breaks:

1. solve with the current breaks and list each station's candidate steps
   ([`break_candidates`](@ref));
2. refit with each of the best `search.tries` candidates added and keep the one
   with the largest F among those whose [`break_test`](@ref) passes, with a
   change in the refit that reaches a size of `search`; a step seen on two
   stations through the baseline between them goes to the station that explains
   more of it;
3. repeat until no candidate passes.

Then every break is tested again against the solve with all of them, and the
weakest one that fails is dropped and the track re-solved, until all pass; a
configured break is held to the F test alone. Every solve uses the per-scan
delays of a solve with the configured breaks ([`search_data`](@ref),
[`search_solve`](@ref)). Progress goes to `io`.
"""
function find_breaks(track::Track; search::BreakSearch = BreakSearch(), io::IO = stdout)
    breaks = Dict(k => copy(v) for (k, v) in track.cfg.breaks)
    data = search_data(track, breaks)
    current = search_solve(track, breaks, data)
    candidates = break_candidates(track, current, search)
    tried = BreakTest{Float64}[]
    dropped = BreakTest{Float64}[]
    changes = Dict{Tuple{String, Float64}, BreakChange{Float64}}()
    while true
        println(io, "round with breaks ", _show_breaks(breaks), ": ", length(candidates), " candidates")
        for c in candidates
            @printf(
                io, "  %-4s %7.2f h  gap F %8.2f  span %s deg, shape %s deg\n", c.station, c.hour, c.gap_F,
                join([@sprintf("%.1f", rad2deg(x)) for x in c.change.span], " "), join([@sprintf("%.2f", rad2deg(x)) for x in c.change.shape], " ")
            )
        end
        flush(io)
        best = nothing
        for c in first(candidates, search.tries)
            trial = search_solve(track, add_break(breaks, c.station, c.hour), data)
            change = break_changes(trial.sol, c.hour, track.origin)[c.station]
            t = break_test(trial, current, c.station, c.hour, search; gap_F = c.gap_F, large = _large(change, search))
            push!(tried, t)
            changes[(c.station, c.hour)] = change
            _print_test(io, t)
            passed(t) && (isnothing(best) || t.F > best.F) && (best = t)
        end
        isnothing(best) && break
        breaks = add_break(breaks, best.station, best.hour)
        current = search_solve(track, breaks, data)
        candidates = break_candidates(track, current, search)
    end
    configured(st, h) = h in get(track.cfg.breaks, st, Float64[])
    while true
        final = BreakTest{Float64}[]
        for (st, hrs) in breaks, h in hrs
            without = search_solve(track, drop_break(breaks, st, h), data)
            proposed = findfirst(t -> (t.station, t.hour) == (st, h), tried)
            change = changes[(st, h)] = break_changes(current.sol, h, track.origin)[st]
            push!(
                final, break_test(
                    current, without, st, h, search;
                    gap_F = isnothing(proposed) ? NaN : tried[proposed].gap_F, large = configured(st, h) || _large(change, search),
                )
            )
        end
        println(io, "final tests")
        foreach(t -> _print_test(io, t), final)
        failing = filter(!passed, final)
        isempty(failing) && return BreakSearchResult(breaks, final, tried, dropped, current, changes)
        weakest = failing[argmin([t.F / t.threshold for t in failing])]
        push!(dropped, weakest)
        println(io, "dropping ", weakest.station, " at ", weakest.hour, " h")
        breaks = drop_break(breaks, weakest.station, weakest.hour)
        current = search_solve(track, breaks, data)
    end
    return
end

_show_breaks(b) = isempty(b) ? "none" : join(("$st $(join(h, ", "))" for (st, h) in sort(collect(b))), "; ")

function _print_test(io, t::BreakTest)
    @printf(
        io, "  %-4s %7.2f h  Δχ² %10.1f  Δp %4d  F %9.2f  (needs %.2f)  %s\n",
        t.station, t.hour, t.dchi2, t.dp, t.F, t.threshold, passed(t) ? "kept" : passed_F(t) ? "small" : "rejected"
    )
    return flush(io)
end

"""
    write_suggestions(dir, track, r::BreakSearchResult)

The breaks [`find_breaks`](@ref) suggests, for review, under `dir` (replaced
whole):

- `suggestions.txt`: one row per break or candidate, with its status
  (`configured`, `suggested`, `dropped` in the final tests, `passed`: it passed
  its test but another candidate of the same round was kept, `small`: it passed
  the F test but its change reaches neither size, or `rejected`), its
  gap and refit F, the size of the change across it per feed (phase shape,
  amplitude scatter, delay; see [`BreakChange`](@ref)), and the sources the
  station observed in the hour before and after; then the `breaks` entry that
  adopts the suggestions, to copy into the campaign file;
- `<status>_<station>_<hour>.png`: the phase and amplitude change across each.

Nothing reads these files back: a break is adopted by adding it to the campaign
file.
"""
function write_suggestions(dir, track::Track, r::BreakSearchResult)
    rm(dir; recursive = true, force = true)
    mkpath(dir)
    configured = track.cfg.breaks
    status(t) = t.hour in get(configured, t.station, Float64[]) ? "configured" :
        any(d -> (d.station, d.hour) == (t.station, t.hour), r.dropped) ? "dropped" :
        any(f -> (f.station, f.hour) == (t.station, t.hour), r.final) ? "suggested" :
        passed(t) ? "passed" : passed_F(t) ? "small" : "rejected"
    latest = Dict((t.station, t.hour) => t for t in vcat(r.tried, r.dropped, r.final))
    order = Dict("configured" => 1, "suggested" => 2, "dropped" => 3, "passed" => 4, "small" => 5, "rejected" => 6)
    rows = sort!(collect(values(latest)); by = t -> (order[status(t)], order[status(t)] > 3 ? -t.F : t.hour))
    ν = r.solve.sol.geom.channel_freqs ./ 1.0e9
    degrees(x) = join([@sprintf("%.1f", rad2deg(v)) for v in x], " ")
    open(joinpath(dir, "suggestions.txt"), "w") do io
        println(io, "Bandpass breaks suggested by `bpfit breaks` for track $(track.cfg.name), for review.")
        println(io, "Hours are since UT midnight before the track's first sample. F is the refit F test")
        println(io, "(against the solve with all suggested breaks for configured, suggested and dropped")
        println(io, "rows; against the solve without it during the search for the others). `passed`")
        println(io, "rows passed their test, but another candidate of the same round explained more,")
        println(io, "often the other station of a baseline that shows the step. `small` rows passed")
        println(io, "the F test with a change below min_shape_deg and min_scatter_percent. The change")
        println(io, "across each break is per feed, with the across-station median removed: span is")
        println(io, "the phase change of the highest channel less the lowest; ramp the part of it a")
        println(io, "delay across the band makes, which the per-scan delays could take up; shape the")
        println(io, "rms phase change beyond that delay and the band mean; scatter the rms amplitude")
        println(io, "change about its band mean. Only candidates whose shape or scatter on the")
        println(io, "residuals reaches min_shape_deg or min_scatter_percent are refit.\n")
        @printf(
            io, "%-10s %-4s %7s %8s %9s %7s  %-11s  %-11s  %-11s  %-11s  %s\n", "status", "stn", "hour", "gap F", "F", "needs",
            "span (deg)", "ramp (deg)", "shape (deg)", "scatter (%)", "sources 1 h before | after"
        )
        for t in rows
            c = r.changes[(t.station, t.hour)]
            @printf(
                io, "%-10s %-4s %7.2f %8.2f %9.2f %7.2f  %-11s  %-11s  %-11s  %-11s  %s\n", status(t), t.station, t.hour, t.gap_F, t.F,
                t.threshold, degrees(c.span), degrees(c.ramp), degrees(c.shape), join([@sprintf("%.1f", 100x) for x in c.scatter], " "),
                _nearby_sources(track.rows, t.station, t.hour)
            )
            save(joinpath(dir, "$(status(t))_$(t.station)_$(t.hour).png"), plot_break_change(c, ν, t, status(t), _nearby_sources(track.rows, t.station, t.hour)))
        end
        println(io, "\nTo adopt the suggestions, set in [tracks.$(track.cfg.name)]:")
        println(io, "breaks = { ", join(("$st = [$(join(h, ", "))]" for (st, h) in sort!(collect(r.breaks); by = first)), ", "), " }")
    end
    return dir
end

# The sources `station` observed in the hour before and the hour after `hour`,
# with scan counts.
function _nearby_sources(rows, station, hour)
    counts(sel) = join(("$(src)×$(n)" for (src, n) in sort!(collect(_tally(r.source for r in rows if station in r.stations && sel(r))); by = x -> -x[2])), " ")
    return counts(r -> hour - 1 <= r.t1 <= hour) * " | " * counts(r -> hour <= r.t0 <= hour + 1)
end

function _tally(xs)
    out = Dict{String, Int}()
    for x in xs
        out[x] = get(out, x, 0) + 1
    end
    return out
end

"""
    plot_break_change(c::BreakChange, ν, t::BreakTest, status, sources) -> Figure

The phase (degrees) and amplitude (percent) change across one break or
candidate against frequency `ν` (GHz), one line per feed.
"""
function plot_break_change(c::BreakChange, ν, t::BreakTest, status, sources)
    feed_colors = ("#2a78d6", "#eb6834")
    fig = Figure(size = (900, 620), backgroundcolor = "#fcfcfb")
    grid = (; xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06), backgroundcolor = "#fcfcfb")
    phase = Axis(fig[1, 1]; ylabel = "phase change (deg)", grid...)
    amp = Axis(fig[2, 1]; ylabel = "amplitude change (%)", xlabel = "frequency (GHz)", grid...)
    for f in axes(c.dp, 2)
        color = feed_colors[mod1(f, length(feed_colors))]
        lines!(
            phase, ν, rad2deg.(view(c.dp, :, f)); color, linewidth = 2,
            label = @sprintf(
                "feed %d: span %.1f°, delay %.3f ns (ramp %.1f°), shape %.1f°", f, rad2deg(c.span[f]), 1.0e9 * c.delay[f],
                rad2deg(c.ramp[f]), rad2deg(c.shape[f])
            )
        )
        lines!(amp, ν, 100 .* view(c.da, :, f); color, linewidth = 2, label = @sprintf("feed %d: scatter %.2f%%", f, 100 * c.scatter[f]))
    end
    axislegend(phase; position = :lt, framevisible = false, labelsize = 11)
    axislegend(amp; position = :lt, framevisible = false, labelsize = 11)
    linkxaxes!(phase, amp)
    Label(fig[0, 1], @sprintf("%s: %s at %.2f h, refit F %.1f (needs %.1f), gap F %.1f", status, t.station, t.hour, t.F, t.threshold, t.gap_F); fontsize = 15)
    Label(fig[3, 1], "sources in the hour before | after: " * sources; fontsize = 11, color = "#5f5e5a")
    return fig
end
