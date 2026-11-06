"""
    write_reports(track, s::TrackSolve; dir = track.cfg.output)

The tables and plots of one solve, under `dir`:

- `timeline.txt`: the scans, the gaps a break could go in, station coverage;
- `chi2.txt`: the goodness of fit of the a-priori and the corrected data;
- `break_delta.txt` and `solution_phase.png`/`solution_amp.png`: how each
  station's bandpass changed across each break;
- `phase_<A>-<B>.png`/`amp_<A>-<B>.png`: scan spectra before and after the
  correction on the busiest baselines of the watched stations, plus two
  baselines that touch none of them.
"""
function write_reports(track::Track, s::TrackSolve; dir = track.cfg.output)
    mkpath(dir)
    cfg = track.cfg
    boundaries = sort!(unique(reduce(vcat, values(s.breaks); init = Float64[])))
    write_timeline(joinpath(dir, "timeline.txt"), track.rows, boundaries)
    keep = in_solve(cfg)
    before, after = scan_averages(track.ps), scan_averages(s.corrected)
    solved = scan_fits(filter(keep, after))
    edges = channel_edges(normalized_residuals(filter(keep, after)))
    open(joinpath(dir, "chi2.txt"), "w") do io
        println(io, "Goodness of fit of the scan averages, parallel hands: each scan's spectrum")
        println(io, "against its own band mean (the weighted total sum of squares per degree of")
        println(io, "freedom). `none` carries the flags and the a-priori scaling but no bandpass.\n")
        left_out = vcat(cfg.apply_only_sources, cfg.apply_only_scans)
        isempty(left_out) || println(io, "Left out of the solve and of these totals: ", join(left_out, ", "), "\n")
        println(io, "model none")
        print_fit(io, track_tss(filter(keep, track.ps)))
        println(io, "\nmodel bandpass")
        print_fit(io, s.chi2; nparams = nparams(s), station_params = station_bandpass_params(s.model, s.sol.geom))
        println(io)
        print_sources(io, scan_fits(before), scan_fits(after), track.rows, left_out)
        println(io)
        print_worst_scans(io, solved, track.rows)
        println(io)
        print_channel_edges(io, edges)
    end
    save(joinpath(dir, "chi2_scans.png"), plot_scan_fits(solved, track.rows))
    save(joinpath(dir, "chi2_sources.png"), plot_source_residuals(normalized_residuals(after), track.rows))
    save(joinpath(dir, "chi2_edges.png"), plot_channel_edges(edges))
    isempty(boundaries) || write_break_deltas(dir, s, boundaries, track.origin)
    selected = plotted_baselines(track, cfg.watch_stations)
    for bl in selected, (quantity, tag) in ((:phase, "phase"), (:amplitude, "amp"))
        save(joinpath(dir, "$(tag)_$(join(bl, "-")).png"), plot_stability(before, after, bl, track.origin; quantity))
    end
    return dir
end

"""
    print_fit(io, q; nparams = 0, station_params = Dict{String, Int}())

A [`track_tss`](@ref) result as tables of the sum of squares per degree of
freedom: the total, each station's baselines (a baseline counts under both its
stations), and each feed pair. The total's degrees of freedom are charged
`nparams`, and each station's its own `station_params[station]`, so that
models with different numbers of parameters compare.
"""
function print_fit(io::IO, q; nparams = 0, station_params = Dict{String, Int}())
    pairs_ = lookup(q.tss, AntennaPair)
    total, dof = sum(q.tss), sum(q.dof)
    @printf(io, "free parameters %d\n  total   %9.4f   (dof %d)\n", nparams, total / (dof - nparams), dof - nparams)
    println(io, "  per station (a baseline counts under both stations)")
    @printf(io, "    %-8s %8s %10s %10s\n", "station", "params", "dof", "tss/dof")
    for st in sort!(unique(reduce(vcat, collect.(pairs_))))
        mine = [p for p in pairs_ if st in p]
        x, d = sum(q.tss[AntennaPair = At(mine)]), sum(q.dof[AntennaPair = At(mine)])
        k = get(station_params, st, 0)
        @printf(io, "    %-8s %8d %10d %10.4f\n", st, k, d - k, x / (d - k))
    end
    println(io, "  per feed pair")
    for fp in lookup(q.tss, FeedPair)
        @printf(io, "    %-8s %10.4f\n", "$(fp[1])-$(fp[2])", sum(q.tss[FeedPair = At(fp)]) / sum(q.dof[FeedPair = At(fp)]))
    end
    return nothing
end

"""
    scan_fits(averages) -> NamedTuple

The parallel-hand sum of squares of each scan average of `averages` about its
band mean, and its degrees of freedom, summed over feed pairs: `tss` and `dof`
over (`Ti`, `AntennaPair`), and the scan name of each `Ti`.
"""
function scan_fits(averages)
    spectra = baseline_spectra(averages)
    q = tss(spectra)
    parallel = [fp for fp in lookup(q.tss, FeedPair) if fp[1] == fp[2]]
    total(A) = dropdims(sum(A[FeedPair = At(parallel)]; dims = FeedPair); dims = FeedPair)
    return (; tss = total(q.tss), dof = total(q.dof), scans = String.(parent(spectra.scan_name)))
end

"""
    print_sources(io, before, after, rows, left_out)

The sum of squares per degree of freedom of each source's scans before and after
the correction ([`scan_fits`](@ref)), not charged for the bandpass parameters.
It grows with a source's signal-to-noise, so bright sources stand highest.
`left_out` lists the sources and scans already left out of the solve.
"""
function print_sources(io::IO, before, after, rows, left_out)
    source = Dict(r.scan => r.source for r in rows)
    println(io, "per source (scan averages, every station; not charged for the bandpass parameters)")
    @printf(io, "    %-14s %6s %10s %10s %10s\n", "source", "scans", "dof", "none", "bandpass")
    for src in sort!(unique(values(source)))
        mine(f) = [k for k in eachindex(f.scans) if get(source, f.scans[k], "") == src]
        ratio(f) = (k = mine(f); sum(f.tss[Ti = k]) / sum(f.dof[Ti = k]))
        k = mine(after)
        isempty(k) && continue
        @printf(
            io, "    %-14s %6d %10d %10.4f %10.4f%s\n", src, length(k), sum(after.dof[Ti = k]), ratio(before), ratio(after),
            src in left_out ? "  (not in the solve)" : ""
        )
    end
    return nothing
end

"""
    print_worst_scans(io, f, rows; factor = 3, limit = 20)

The scans in which a station's baselines fit worst after the correction: the sum
of squares per degree of freedom of the station's baselines in the scan against
its median over the station's scans, for the `limit` scans furthest above it
and at least `factor` times it, with the baseline that carries the largest share
of the station's sum of squares in the scan. A scan whose excess sits on one
baseline is a candidate for flagging that baseline (`scan_flags` with both
stations); one spread over the station's baselines, for flagging the station.
A run of them may be a station whose bandpass changed within the scans, which
breaks on either side of each give a bandpass of its own.
"""
function print_worst_scans(io::IO, f, rows; factor = 3, limit = 20)
    hour = Dict(r.scan => r.t0 for r in rows)
    source = Dict(r.scan => r.source for r in rows)
    pairs_ = lookup(f.tss, AntennaPair)
    found = @NamedTuple{station::String, scan::String, ratio::Float64, typical::Float64, baseline::String, share::Float64}[]
    for st in sort!(unique(reduce(vcat, collect.(pairs_))))
        mine = [p for p in pairs_ if st in p]
        x = vec(parent(sum(f.tss[AntennaPair = At(mine)]; dims = AntennaPair)))
        d = vec(parent(sum(f.dof[AntennaPair = At(mine)]; dims = AntennaPair)))
        k = findall(>(0), d)
        length(k) >= 3 || continue
        typical = median(x[k] ./ d[k])
        for j in k
            x[j] / d[j] >= factor * typical || continue
            each = vec(parent(f.tss[Ti = j, AntennaPair = At(mine)]))
            worst = argmax(each)
            push!(found, (; station = String(st), scan = f.scans[j], ratio = x[j] / d[j], typical, baseline = join(mine[worst], "-"), share = each[worst] / x[j]))
        end
    end
    sort!(found; by = r -> -r.ratio / r.typical)
    println(io, "scans that fit worst on a station, against the station's median scan (at least $(factor)×)")
    isempty(found) && return println(io, "    none")
    @printf(
        io, "    %-8s %-18s %-14s %8s %10s %10s %8s  %s\n", "station", "scan", "source", "hour", "tss/dof", "median", "ratio",
        "worst baseline (share of the station's sum of squares)"
    )
    for r in first(found, limit)
        @printf(
            io, "    %-8s %-18s %-14s %8.3f %10.3f %10.3f %7.1f×  %s (%.0f%%)\n", r.station, r.scan, source[r.scan], hour[r.scan], r.ratio,
            r.typical, r.ratio / r.typical, r.baseline, 100 * r.share
        )
    end
    length(found) > limit && println(io, "    … and $(length(found) - limit) more")
    return nothing
end

"""
    channel_edges(residuals; min_scans = 5) -> NamedTuple

How much the corrected phase step from each channel to the next changes from
scan to scan, per baseline. On each baseline and parallel feed pair of the
[`normalized_residuals`](@ref) `residuals` of scan averages, the step is taken
in every scan with at least half its channels; its spread over the scans (the
median absolute deviation, over at least `min_scans` scans) is divided by the
median spread over the edges. `ratio` is (baseline × edge), the larger of the
feed pairs, `NaN` where unmeasured; `baselines` are the antenna pairs and edge
`c` lies between channels `c` and `c + 1`.
"""
function channel_edges(residuals; min_scans = 5)
    r = residuals.r
    pairs_, feeds = lookup(r, AntennaPair), lookup(r, FeedPair)
    ratio = fill(NaN, length(pairs_), size(r, Frequency) - 1)
    for (row, j) in enumerate(axes(r, AntennaPair)), q in axes(r, FeedPair)
        feeds[q][1] == feeds[q][2] || continue
        R = parent(permutedims(r[AntennaPair(j), FeedPair(q)], (Frequency, Ti)))
        scans = [k for k in axes(R, 2) if 2 * count(isfinite, view(R, :, k)) > size(R, 1)]
        length(scans) >= min_scans || continue
        spread = map(firstindex(R, 1):(lastindex(R, 1) - 1)) do c
            d = [rem2pi(angle(R[c + 1, k]) - angle(R[c, k]), RoundNearest) for k in scans if isfinite(R[c, k]) && isfinite(R[c + 1, k])]
            length(d) >= min_scans ? median(abs.(d .- median(d))) : NaN
        end
        known = filter(isfinite, spread)
        isempty(known) && continue
        for c in eachindex(spread)
            v = spread[c] / median(known)
            isfinite(v) && (ratio[row, c] = isfinite(ratio[row, c]) ? max(ratio[row, c], v) : v)
        end
    end
    return (; baselines = collect(pairs_), ratio)
end

"""
    print_channel_edges(io, edges; factor = 4)

The channel edges of [`channel_edges`](@ref) `edges` at least `factor` times
their baseline's median edge. A station whose phase slope restarts at fixed
channel blocks and changes from scan to scan shows the same edges on its
baselines; fit per-scan terms over those blocks (`scan_terms`).
"""
function print_channel_edges(io::IO, edges; factor = 4)
    println(io, "channel edges where a baseline's phase steps differently from scan to scan (at least $(factor)× its median edge)")
    found = false
    for (row, bl) in pairs(edges.baselines)
        hits = [c for c in axes(edges.ratio, 2) if edges.ratio[row, c] >= factor]
        isempty(hits) && continue
        found = true
        @printf(io, "    %-8s %s\n", join(bl, "-"), join(("$(c)|$(c + 1) $(@sprintf("%.1f", edges.ratio[row, c]))×" for c in hits), ", "))
    end
    found || println(io, "    none")
    return nothing
end

const DIVERGING = ["#104281", "#3987e5", "#9ec5f4", "#f0efec", "#f4b1b0", "#e34948", "#9a2626"]
const SEQUENTIAL = ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
const FEED_COLORS = ("#2a78d6", "#eb6834")
const SURFACE = "#fcfcfb"

"""
    plot_scan_fits(f, rows) -> Figure

Each station's fit in each scan ([`scan_fits`](@ref) `f`), scans in time order:
the sum of squares per degree of freedom of the station's baselines in the scan
over its median across the station's scans, on a log scale. The sources
observed are labelled where they change.
"""
function plot_scan_fits(f, rows)
    hour = Dict(r.scan => r.t0 for r in rows)
    source = Dict(r.scan => r.source for r in rows)
    order = sortperm(f.scans; by = sc -> get(hour, sc, Inf))
    scans = f.scans[order]
    pairs_ = lookup(f.tss, AntennaPair)
    stations = sort!(unique(reduce(vcat, collect.(pairs_))))
    z = fill(NaN, length(scans), length(stations))
    for (i, st) in pairs(stations)
        mine = [p for p in pairs_ if st in p]
        x = vec(parent(sum(f.tss[AntennaPair = At(mine)]; dims = AntennaPair)))[order]
        d = vec(parent(sum(f.dof[AntennaPair = At(mine)]; dims = AntennaPair)))[order]
        k = findall(>(0), d)
        isempty(k) && continue
        typical = median(x[k] ./ d[k])
        z[k, i] = log10.((x[k] ./ d[k]) ./ typical)
    end
    starts = [k for k in eachindex(scans) if k == firstindex(scans) || source[scans[k]] != source[scans[k - 1]]]
    fig = Figure(size = (max(900, 12 * length(scans) + 200), 140 + 28 * length(stations)), backgroundcolor = SURFACE)
    ax = Axis(
        fig[1, 1]; backgroundcolor = SURFACE, yticks = (eachindex(stations), stations),
        xticks = (starts, [source[scans[k]] for k in starts]), xticklabelrotation = π / 2, xticklabelsize = 10,
        xlabel = "scans in time order (labelled where the source changes)", title = "station fit per scan, against the station's median scan",
    )
    hm = heatmap!(ax, eachindex(scans), eachindex(stations), clamp.(z, -1.5, 1.5); colormap = DIVERGING, colorrange = (-1.5, 1.5), nan_color = SURFACE)
    ticks = [0.03, 0.1, 0.3, 1, 3, 10, 30]
    Colorbar(fig[1, 2], hm; label = "tss/dof ÷ station median", ticks = (log10.(ticks), string.(ticks)))
    return fig
end

"""
    plot_source_residuals(residuals, rows) -> Figure

The corrected residual spectrum of each source, one row per source: on each
parallel feed pair the inverse-variance weighted mean over the source's scans
and baselines of the [`normalized_residuals`](@ref) `residuals`, as amplitude
(percent of the band mean) and phase (degrees). The bandpass is common to all
sources, so a source with a spectral feature of its own pulls the solve toward
it: its mean residual shows the feature with one sign and the other sources'
with the opposite sign.
"""
function plot_source_residuals(residuals, rows)
    r, var = residuals.r, residuals.var
    source = Dict(rw.scan => rw.source for rw in rows)
    scans = parent(residuals.scan_name)
    ν = collect(lookup(r, Frequency)) ./ 1.0e9
    feeds = lookup(r, FeedPair)
    sources = sort!(unique([source[sc] for sc in scans if haskey(source, sc)]))
    fig = Figure(size = (1000, 60 + 170 * length(sources)), backgroundcolor = SURFACE)
    grid = (; xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06), backgroundcolor = SURFACE)
    for (i, src) in pairs(sources)
        ks = [k for k in eachindex(scans) if get(source, scans[k], "") == src]
        amp = Axis(fig[i, 1]; ylabel = "$src\namplitude (%)", xlabel = i == length(sources) ? "frequency (GHz)" : "", grid...)
        phase = Axis(fig[i, 2]; ylabel = "phase (deg)", xlabel = i == length(sources) ? "frequency (GHz)" : "", grid...)
        for q in axes(r, FeedPair)
            f = feeds[q][1]
            f == feeds[q][2] || continue
            m = map(eachindex(ν)) do c
                R = vec(parent(r[Ti = ks, FeedPair = q, Frequency = c]))
                W = inv.(vec(parent(var[Ti = ks, FeedPair = q, Frequency = c])))
                ok = [i for i in eachindex(R, W) if isfinite(R[i]) && isfinite(W[i])]
                isempty(ok) ? complex(NaN, NaN) : sum(W[ok] .* R[ok]) / sum(W[ok])
            end
            color = FEED_COLORS[mod1(f, length(FEED_COLORS))]
            lines!(amp, ν, 100 .* (abs.(m) .- 1); color, linewidth = 2, label = "feed $f")
            lines!(phase, ν, rad2deg.(angle.(m)); color, linewidth = 2, label = "feed $f")
        end
        i == 1 && axislegend(amp; position = :lt, framevisible = false, labelsize = 11)
    end
    Label(fig[0, 1:2], "mean corrected residual spectrum per source"; fontsize = 15)
    return fig
end

"""
    plot_channel_edges(edges) -> Figure

The [`channel_edges`](@ref) of every baseline: how much the phase step across
each channel edge changes from scan to scan, over the baseline's median edge.
"""
function plot_channel_edges(edges)
    rows_ = [row for row in axes(edges.ratio, 1) if any(isfinite, view(edges.ratio, row, :))]
    labels = [join(edges.baselines[row], "-") for row in rows_]
    fig = Figure(size = (900, 140 + 16 * length(rows_)), backgroundcolor = SURFACE)
    ax = Axis(
        fig[1, 1]; backgroundcolor = SURFACE, yticks = (eachindex(rows_), labels), yticklabelsize = 9,
        xticks = 1:2:size(edges.ratio, 2), xlabel = "channel edge c (between channels c and c + 1)",
        title = "scan-to-scan spread of the phase step at each channel edge, over the baseline's median edge",
    )
    z = permutedims(log10.(clamp.(edges.ratio[rows_, :], 1, 30)))
    hm = heatmap!(ax, axes(edges.ratio, 2), eachindex(rows_), z; colormap = SEQUENTIAL, colorrange = (0, log10(30)), nan_color = SURFACE)
    ticks = [1, 2, 4, 10, 30]
    Colorbar(fig[1, 2], hm; label = "spread ÷ median edge", ticks = (log10.(ticks), string.(ticks)))
    return fig
end

function plotted_baselines(track, watch; partners = 2, quiet = 2)
    selected = Tuple{String, String}[]
    counts = baseline_scan_counts(track.ps)
    stations = Gustavo.Calibration.DataGeometry(track.ps).antennas
    for st in watch
        st in stations || (@warn "watch station absent from the track" station = st; continue)
        mine = filter(p -> st in first(p), counts)
        append!(selected, first.(first(mine, partners)))
    end
    for (bl, _) in counts
        length(selected) >= partners * length(watch) + quiet && break
        any(in(watch), bl) || bl in selected || push!(selected, bl)
    end
    return selected
end

"""
    segment_picks(times, boundaries, origin) -> (picks, labels)

One representative index into `times` per interval of `boundaries` (hours since
`origin`), and the label each carries. An interval the track never reaches picks
`nothing`.
"""
function segment_picks(times, boundaries, origin)
    edges = vcat(-Inf, collect(boundaries), Inf)
    hours = track_hours(times, origin)
    picks = [findfirst(h -> edges[k] <= h < edges[k + 1], hours) for k in 1:(length(edges) - 1)]
    labels = [
        k == 1 ? "< $(edges[2]) h" : k == length(picks) ? "≥ $(edges[k]) h" : "$(edges[k])–$(edges[k + 1]) h"
            for k in eachindex(picks)
    ]
    return picks, labels
end

function gain_array(sol)
    G = gains(sol)
    return parent(permutedims(G, (Frequency, Ti, AntennaName, Feed))), String.(lookup(G, AntennaName)), collect(lookup(G, Ti))
end

"""
    plot_bandpass_solution(sol, boundaries, origin; quantity = :phase) -> Figure

The solved bandpass against frequency, one row per station and one column per
feed, one trace per interval of `boundaries`. A station whose traces lie on top
of one another did not change; a station with no break repeats one trace.
"""
function plot_bandpass_solution(sol, boundaries, origin; quantity::Symbol = :phase)
    G, names, times = gain_array(sol)
    picks, labels = segment_picks(times, boundaries, origin)
    keep = findall(!isnothing, picks)
    fghz = sol.geom.channel_freqs ./ 1.0e9
    reduce_y = quantity === :phase ? angle : abs
    ylab = quantity === :phase ? "phase (rad)" : "|gain|"
    nfeed = size(G, 4)
    fig = Figure(size = (520 * nfeed + 60, 200 * length(names) + 60))
    axs = Axis[]
    for (row, a) in pairs(names), fi in 1:nfeed
        ax = Axis(fig[row, fi]; xlabel = row == length(names) ? "frequency (GHz)" : "", ylabel = a, title = row == 1 ? "feed $(fi)  $ylab" : "")
        for k in keep
            scatter!(ax, fghz, reduce_y.(view(G, :, picks[k], row, fi)); markersize = 5, label = labels[k])
        end
        push!(axs, ax)
    end
    linkxaxes!(axs...)
    Legend(fig[0, 1:nfeed], first(axs); orientation = :horizontal, framevisible = false)
    return fig
end

function best_delay(v, w, ν)
    best = (s = -Inf, τ = 0.0)
    for τ in range(-0.5, 0.5; length = 2001) .* (length(ν) / (maximum(ν) - minimum(ν)))
        s = abs(sum(w .* v .* cis.(-2π * τ .* ν)))
        s > best.s && (best = (; s, τ))
    end
    return best.τ
end

"""
    BreakChange

How one station's bandpass moves across a break: per channel and feed the phase
change `dp` (rad) and log-amplitude change `da`; per feed the `delay` (s) the
phase change carries, its `shape` (rad), the rms of the phase change beyond that
delay and its band mean, and `scatter`, the rms of `da` about its band mean.
Channels with `NaN` changes are left out of the per-feed values.
"""
struct BreakChange{T}
    dp::Matrix{T}
    da::Matrix{T}
    delay::Vector{T}
    shape::Vector{T}
    scatter::Vector{T}
end

"""
    break_changes(sol, hour, origin) -> Dict{String, BreakChange}

Every station's [`BreakChange`](@ref) across `hour` (hours since `origin`) in
the bandpass shape of `sol`, from its last sample before `hour` to its first
after, with the across-station median removed channel by channel: that median
is the gauge difference between two independently solved segments, not a
change.
"""
function break_changes(sol, hour, origin)
    shape = shape_only(sol)
    after = searchsortedfirst(track_hours(shape.geom.times, origin), hour)
    G = gains(shape; Ti = [after - 1, after])
    g = parent(permutedims(G, (Frequency, Ti, AntennaName, Feed)))
    dp = strip_common(rem2pi.(angle.(view(g, :, 2, :, :)) .- angle.(view(g, :, 1, :, :)), RoundNearest))
    da = strip_common(log.(abs.(view(g, :, 2, :, :))) .- log.(abs.(view(g, :, 1, :, :))))
    ν = shape.geom.channel_freqs .- mean(shape.geom.channel_freqs)
    return Dict(
        String(st) => _station_change(view(dp, :, a, :), view(da, :, a, :), ν)
            for (a, st) in pairs(lookup(G, AntennaName))
    )
end

function _station_change(dp, da, ν)
    stats = map(f -> _feed_change(view(dp, :, f), view(da, :, f), ν), axes(dp, 2))
    return BreakChange(collect(dp), collect(da), [x.τ for x in stats], [x.shape for x in stats], [x.scatter for x in stats])
end

function _feed_change(dp, da, ν)
    ok = [c for c in eachindex(dp, da, ν) if isfinite(dp[c]) && isfinite(da[c])]
    length(ok) < 2 && return (; τ = NaN, shape = NaN, scatter = NaN)
    z = cis.(dp[ok])
    τ = best_delay(z, ones(length(ok)), ν[ok])
    beyond = z .* cis.(-2π * τ .* ν[ok])
    return (; τ, shape = rms(angle.(beyond .* conj(sum(beyond)))), scatter = std(da[ok]))
end

"""
    write_break_deltas(dir, s, boundaries, origin)

`break_delta.txt` and the per-segment solution plots: how far each station's
bandpass moved across each boundary ([`break_changes`](@ref)).
"""
function write_break_deltas(dir, s::TrackSolve, boundaries, origin)
    shape = shape_only(s.sol)
    for (quantity, tag) in ((:phase, "phase"), (:amplitude, "amp"))
        save(joinpath(dir, "solution_$(tag).png"), plot_bandpass_solution(shape, boundaries, origin; quantity))
    end
    edges = vcat(-Inf, collect(boundaries), Inf)
    label(k) = k == 1 ? "< $(edges[2]) h" : k == length(edges) - 1 ? "≥ $(edges[k]) h" : "$(edges[k])–$(edges[k + 1]) h"
    open(joinpath(dir, "break_delta.txt"), "w") do io
        println(io, "Change across each boundary. Boundaries (h): $(join(boundaries, ", "))")
        println(io, "The across-station median is removed channel by channel; it is the gauge")
        println(io, "difference between two independently solved segments, not a change.")
        println(io, "Per station and feed: the delay the phase change carries, the rms of the phase")
        println(io, "change beyond that delay and its band mean, the rms of the log-amplitude change")
        println(io, "about its band mean, and for each feed after the first the change of its band-mean")
        println(io, "phase from the first feed's.\n")
        for (k, h) in enumerate(boundaries)
            println(io, "$(label(k))  →  $(label(k + 1))")
            println(io, "station feed  delay(ns)  shape(deg)  dlogamp(%)  feed-phase(deg)")
            for (st, c) in sort!(collect(break_changes(s.sol, h, origin)); by = first), f in axes(c.dp, 2)
                feed = f == 1 ? 0.0 : angle(sum(cis.(view(c.dp, :, f))) * conj(sum(cis.(view(c.dp, :, 1)))))
                @printf(io, "%-7s %4d  %9.3f  %10.2f  %10.2f  %15.1f\n", st, f, 1.0e9 * c.delay[f], rad2deg(c.shape[f]), 100 * c.scatter[f], rad2deg(feed))
            end
            println(io)
        end
    end
    return nothing
end

# One scan's spectrum on baseline `bi`, product `p`, divided by its weighted band
# mean, with 1σ errors of amplitude and phase from its weights. `nothing` when its
# median channel signal-to-noise is below `min_snr`.
function normalized_spectrum(ms, bi, p; min_snr = 3)
    vis = ms[:visibility]
    ν = collect(lookup(vis, Frequency))
    v = zeros(ComplexF64, length(ν))
    w = zeros(length(ν))
    U = usable_weights(ms)
    order = dimnum(vis, (Frequency, Ti, BaselineID, Polarization))
    V = PermutedDimsArray(parent(vis), order)
    W = PermutedDimsArray(U, order)
    for t in axes(V, 2), c in axes(V, 1)
        wi = W[c, t, bi, p]
        vi = V[c, t, bi, p]
        (wi > 0 && isfinite(vi)) || continue
        v[c] += wi * vi
        w[c] += wi
    end
    ok = w .> 0
    count(ok) > 1 || return nothing
    v[ok] ./= w[ok]
    m = sum(w[ok] .* v[ok]) / sum(w[ok])
    σ = sqrt.(1 ./ w[ok])
    median(abs(m) ./ σ) >= min_snr || return nothing
    return (; ν = ν[ok], r = v[ok] ./ m, σamp = σ ./ abs(m), σph = σ ./ abs.(v[ok]), t = mean(lookup(vis, Ti)))
end

"""
    plot_stability(before, after, baseline, origin; quantity = :phase) -> Figure

Each scan's average spectrum on `baseline`, parallel hands, before (left) and
after (right) the correction (each a processing set, at any time resolution;
[`write_reports`](@ref) passes scan averages), colored by time, for scans whose median channel
signal-to-noise is at least 3. Each spectrum is divided by its weighted band
mean, which leaves the structure across the band: the bandpass shape, and any
delay the correction left, as a phase slope. The black line is the
inverse-variance weighted mean.
"""
function plot_stability(before, after, baseline, origin; quantity::Symbol = :phase)
    fig = Figure(size = (1100, 650))
    hours(t) = track_hours(t, origin)
    trange = extrema(hours(mean(lookup(ms[:visibility], Ti))) for ms in values(after))
    nfeed = maximum(ms -> maximum(first, feed_pairs(ms)), values(after))
    for f in 1:nfeed
        row = Axis[]
        for (col, (data, title)) in enumerate(((before, "before"), (after, "after")))
            ax = Axis(
                fig[f, col]; title = "$(join(baseline, "-"))  $title  feeds $f$f",
                xlabel = f == nfeed ? "frequency (GHz)" : "",
                ylabel = quantity === :phase ? "phase (deg)" : "amplitude / band mean",
            )
            push!(row, ax)
            sums = Dict{Float64, Tuple{ComplexF64, Float64, Float64}}()
            for ms in values(data)
                bi = findfirst(==(baseline), [(String(a), String(b)) for (a, b) in XRadio.baselines(ms)])
                isnothing(bi) && continue
                p = findfirst(==((f, f)), feed_pairs(ms)[:, bi])
                isnothing(p) && continue
                s = normalized_spectrum(ms, bi, p)
                isnothing(s) && continue
                x = s.ν ./ 1.0e9
                y, e = quantity === :phase ? (rad2deg.(angle.(s.r)), rad2deg.(s.σph)) : (abs.(s.r), s.σamp)
                c = fill(hours(s.t), length(x))
                errorbars!(ax, x, y, e; color = c, colormap = :viridis, colorrange = trange, linewidth = 0.6, alpha = 0.5)
                scatter!(ax, x, y; color = c, colormap = :viridis, colorrange = trange, markersize = 3)
                for k in eachindex(x)
                    w = 1 / s.σamp[k]^2
                    z, a, ww = get(sums, x[k], (0.0im, 0.0, 0.0))
                    sums[x[k]] = (z + w * s.r[k], a + w * abs(s.r[k]), ww + w)
                end
            end
            isempty(sums) && continue
            xs = sort(collect(keys(sums)))
            ym = [quantity === :phase ? rad2deg(angle(sums[x][1])) : sums[x][2] / sums[x][3] for x in xs]
            lines!(ax, xs, ym; color = :black, linewidth = 1.5)
        end
        linkyaxes!(row...)
    end
    Colorbar(fig[1:nfeed, 3]; colormap = :viridis, limits = trange, label = "hours from UT midnight, first day")
    return fig
end
