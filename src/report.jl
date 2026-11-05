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
    end
    isempty(boundaries) || write_break_deltas(dir, s, boundaries, track.origin)
    selected = plotted_baselines(track, cfg.watch_stations)
    before, after = scan_averages(track.ps), scan_averages(s.corrected)
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
