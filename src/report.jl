"""
    write_reports(track, s::TrackSolve; dir = track.cfg.output)

The tables and plots of one solve, under `dir`:

- `timeline.txt`: the scans, the gaps a break could go in, station coverage;
- `chi2.txt`: the goodness of fit of the a-priori and the corrected data;
- `break_delta.txt` and `solution_phase.png`/`solution_amp.png`: how each
  station's bandpass changed across each break;
- `phase_<A>-<B>.png`/`amp_<A>-<B>.png`: scan spectra before and after the
  correction on the busiest baselines of the watched stations, the baselines
  that fit worst after it, and those the bandpass changes most.
"""
function write_reports(track::Track, s::TrackSolve; dir = track.cfg.output)
    mkpath(dir)
    cfg = track.cfg
    boundaries = sort!(unique(reduce(vcat, values(s.breaks); init = Float64[])))
    write_timeline(joinpath(dir, "timeline.txt"), track.rows, boundaries)
    keep = in_solve(cfg)
    before, after = scan_averages(track.ps), scan_averages(s.corrected)
    solved = scan_fits(filter(keep, after))
    found = worst_scans(solved)
    baselines = baseline_summary(scan_fits(filter(keep, before)), solved, s.sol)
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
        print_baselines(io, baselines)
        println(io)
        print_worst_scans(io, found, track.rows)
        println(io)
        print_channel_edges(io, edges)
    end
    save(joinpath(dir, "chi2_scans.png"), plot_scan_fits(solved, found, track.rows))
    save(joinpath(dir, "chi2_baselines.png"), plot_baselines(baselines))
    save(joinpath(dir, "chi2_sources.png"), plot_source_residuals(normalized_residuals(after), track.rows))
    save(joinpath(dir, "chi2_edges.png"), plot_channel_edges(edges))
    isempty(boundaries) || write_break_deltas(dir, s, boundaries, track.origin)
    selected = plotted_baselines(track, cfg.watch_stations, baselines)
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
    baseline_summary(before, after, sol; min_scans = 3) -> Vector{NamedTuple}

Every baseline of the scan fits `after` ([`scan_fits`](@ref)) with data in at
least `min_scans` scans, worst fit first: its scans, its sum of squares per
degree of freedom `before` and `after` the correction, and how much the solved
bandpass `sol` changes it ([`correction_sizes`](@ref)).
"""
function baseline_summary(before, after, sol; min_scans = 3)
    ratio(f, bl) = sum(f.tss[AntennaPair = At(bl)]) / sum(f.dof[AntennaPair = At(bl)])
    pairs_ = [bl for bl in lookup(after.tss, AntennaPair) if count(>(0), parent(after.dof[AntennaPair = At(bl)])) >= min_scans]
    size = correction_sizes(sol, [(String(bl[1]), String(bl[2])) for bl in pairs_])
    rows = map(pairs_) do bl
        name = (String(bl[1]), String(bl[2]))
        (;
            baseline = name, scans = count(>(0), parent(after.dof[AntennaPair = At(bl)])),
            before = bl in lookup(before.tss, AntennaPair) ? ratio(before, bl) : NaN, after = ratio(after, bl), correction = size[name],
        )
    end
    return sort!(rows; by = r -> -r.after)
end

"""
    print_baselines(io, rows)

The [`baseline_summary`](@ref) `rows` as a table.
"""
function print_baselines(io::IO, rows)
    println(io, "per baseline (scans in the solve; the correction is the rms phase of the solved bandpass about its band mean)")
    @printf(io, "    %-8s %6s %12s %12s %16s\n", "baseline", "scans", "none", "bandpass", "correction (deg)")
    for r in rows
        @printf(io, "    %-8s %6d %12.3f %12.3f %16.1f\n", join(r.baseline, "-"), r.scans, r.before, r.after, rad2deg(r.correction))
    end
    return nothing
end

"""
    plot_baselines(rows) -> Figure

The [`baseline_summary`](@ref) `rows` as three station-by-station grids: the
fit before the correction, after it (tss/dof on a log scale from 1, the noise),
and the size of the correction (degrees).
"""
function plot_baselines(rows)
    stations = sort!(unique(reduce(vcat, [collect(r.baseline) for r in rows]; init = String[])))
    at = Dict(st => i for (i, st) in pairs(stations))
    grid(v) = (z = fill(NaN, length(stations), length(stations)); foreach(r -> (z[at[r.baseline[1]], at[r.baseline[2]]] = z[at[r.baseline[2]], at[r.baseline[1]]] = v(r)), rows); z)
    n = length(stations)
    fig = Figure(size = (360 * 3 + 120, 120 + 30 * n), backgroundcolor = SURFACE)
    ticks = (eachindex(stations), stations)
    top = 100.0
    fits = [1, 3, 10, 30, 100]
    for (col, (title, z)) in pairs(
            (
                ("fit before the correction (tss/dof)", grid(r -> log10(clamp(r.before, 1, top)))),
                ("fit after the correction (tss/dof)", grid(r -> log10(clamp(r.after, 1, top)))),
            )
        )
        ax = Axis(fig[1, col]; title, titlesize = 12, xticks = ticks, yticks = ticks, aspect = DataAspect(), backgroundcolor = SURFACE, xticklabelrotation = π / 2)
        hm = heatmap!(ax, 1:n, 1:n, z; colormap = SEQUENTIAL, colorrange = (0, log10(top)), nan_color = SURFACE)
        col == 2 && Colorbar(fig[1, 3], hm; label = "tss/dof", ticks = (log10.(fits), string.(fits)))
    end
    corr = grid(r -> rad2deg(r.correction))
    ax = Axis(fig[1, 4]; title = "bandpass correction (deg rms)", titlesize = 12, xticks = ticks, yticks = ticks, aspect = DataAspect(), backgroundcolor = SURFACE, xticklabelrotation = π / 2)
    hm = heatmap!(ax, 1:n, 1:n, corr; colormap = SEQUENTIAL, colorrange = (0, max(1.0, maximum(filter(isfinite, corr); init = 1.0))), nan_color = SURFACE)
    Colorbar(fig[1, 5], hm; label = "degrees")
    return fig
end

"""
    worst_scans(f; factor = 3) -> Vector{NamedTuple}

The scans in which a station's baselines fit worst after the correction
([`scan_fits`](@ref) `f`): every (station, scan) whose sum of squares per degree
of freedom over the station's baselines is at least `factor` times its median
over the station's scans, furthest above it first, with the baseline that
carries the largest share of the station's sum of squares in the scan.
"""
function worst_scans(f; factor = 3)
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
    return sort!(found; by = r -> -r.ratio / r.typical)
end

"""
    print_worst_scans(io, found, rows; factor = 3, limit = 20)

The first `limit` of the [`worst_scans`](@ref) `found` (at least `factor` times
their station's median). A scan whose excess sits on one baseline is a
candidate for flagging that baseline (`scan_flags` with both stations); one
spread over the station's baselines, for flagging the station. A run of them
may be a station whose bandpass changed within the scans, which breaks on
either side of each give a bandpass of its own.
"""
function print_worst_scans(io::IO, found, rows; factor = 3, limit = 20)
    hour = Dict(r.scan => r.t0 for r in rows)
    source = Dict(r.scan => r.source for r in rows)
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
    shared_edges(edges; factor = 4) -> Vector{NamedTuple}

The stations whose baselines share channel edges of [`channel_edges`](@ref)
`edges` at least `factor` times their baseline's median edge: per station with
such edges on two or more of its baselines, the edges flagged on two or more of
them and how many baselines flag any, most first.
"""
function shared_edges(edges; factor = 4)
    stations = sort!(unique(reduce(vcat, collect.(edges.baselines))))
    out = @NamedTuple{station::String, edges::Vector{Int}, baselines::Int}[]
    for st in stations
        rows_ = [row for (row, bl) in pairs(edges.baselines) if st in bl]
        flagged = [row for row in rows_ if any(>=(factor), filter(isfinite, view(edges.ratio, row, :)))]
        length(flagged) >= 2 || continue
        common = [c for c in axes(edges.ratio, 2) if count(row -> edges.ratio[row, c] >= factor, flagged) >= 2]
        isempty(common) || push!(out, (; station = String(st), edges = common, baselines = length(flagged)))
    end
    return sort!(out; by = x -> -x.baselines)
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
    for x in shared_edges(edges; factor)
        println(io, "  shared by $(x.baselines) baselines of $(x.station): ", join(("$(c)|$(c + 1)" for c in x.edges), ", "))
    end
    return nothing
end

const DIVERGING = ["#104281", "#3987e5", "#9ec5f4", "#f0efec", "#f4b1b0", "#e34948", "#9a2626"]
const SEQUENTIAL = ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
const FEED_COLORS = ("#2a78d6", "#eb6834")
const SURFACE = "#fcfcfb"

"""
    plot_scan_fits(f, found, rows; label = 6) -> Figure

Each station's fit in each scan ([`scan_fits`](@ref) `f`), scans in time order:
the sum of squares per degree of freedom of the station's baselines in the scan,
on a log scale from 1, the noise. The first `label` of the
[`worst_scans`](@ref) `found` are numbered on the plot and listed under it, and
the sources observed are labelled where they change.
"""
function plot_scan_fits(f, found, rows; label = 6)
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
        z[k, i] = log10.(max.(x[k] ./ d[k], 1))
    end
    starts = [k for k in eachindex(scans) if k == firstindex(scans) || source[scans[k]] != source[scans[k - 1]]]
    marked = first(found, label)
    fig = Figure(size = (max(900, 12 * length(scans) + 200), 200 + 28 * length(stations) + 18 * length(marked)), backgroundcolor = SURFACE)
    ax = Axis(
        fig[1, 1]; backgroundcolor = SURFACE, yticks = (eachindex(stations), stations),
        xticks = (starts, [source[scans[k]] for k in starts]), xticklabelrotation = π / 2, xticklabelsize = 10,
        xlabel = "scans in time order (labelled where the source changes)",
        title = "station fit per scan: tss/dof of the station's baselines, 1 is the noise; numbered, the scans furthest above their station's median",
    )
    top = 100.0
    hm = heatmap!(ax, eachindex(scans), eachindex(stations), clamp.(z, 0, log10(top)); colormap = SEQUENTIAL, colorrange = (0, log10(top)), nan_color = SURFACE)
    notes = String[]
    for (n, r) in pairs(marked)
        x, y = findfirst(==(r.scan), scans), findfirst(==(r.station), stations)
        (isnothing(x) || isnothing(y)) && continue
        text!(ax, x, y; text = string(n), align = (:center, :center), fontsize = 11, color = :black, font = :bold)
        push!(notes, @sprintf("%d  %s in %s (%s, %.2f h): %.1f× its median scan, %.0f%% on %s", n, r.station, r.scan, source[r.scan], hour[r.scan], r.ratio / r.typical, 100 * r.share, r.baseline))
    end
    ticks = [1, 3, 10, 30, 100]
    Colorbar(fig[1, 2], hm; label = "tss/dof", ticks = (log10.(ticks), string.(ticks)))
    isempty(notes) || Label(fig[2, 1:2], join(notes, "\n"); halign = :left, justification = :left, fontsize = 11)
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
    parallel = [q for q in axes(r, FeedPair) if feeds[q][1] == feeds[q][2]]
    mean_residual(ks, q) = map(eachindex(ν)) do c
        R = vec(parent(r[Ti = ks, FeedPair = q, Frequency = c]))
        W = inv.(vec(parent(var[Ti = ks, FeedPair = q, Frequency = c])))
        ok = [i for i in eachindex(R, W) if isfinite(R[i]) && isfinite(W[i])]
        isempty(ok) ? complex(NaN, NaN) : sum(W[ok] .* R[ok]) / sum(W[ok])
    end
    m = Dict(src => [mean_residual([k for k in eachindex(scans) if get(source, scans[k], "") == src], q) for q in parallel] for src in sources)
    amplitude(src) = [mean(q -> abs(m[src][q][c]) - 1, eachindex(parallel)) for c in eachindex(ν)]
    against = Dict(src => _correlation(amplitude(src), [mean(o -> amplitude(o)[c], filter(!=(src), sources)) for c in eachindex(ν)]) for src in sources if length(sources) > 1)
    opposite = sort!([src for (src, ρ) in against if ρ < -0.5]; by = src -> against[src])
    fig = Figure(size = (1000, 80 + 170 * length(sources)), backgroundcolor = SURFACE)
    grid = (; xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06), backgroundcolor = SURFACE)
    for (i, src) in pairs(sources)
        ρ = get(against, src, NaN)
        amp = Axis(
            fig[i, 1]; ylabel = "$src\namplitude (%)", xlabel = i == length(sources) ? "frequency (GHz)" : "",
            title = @sprintf("r = %.2f with the other sources", ρ), titlesize = 11, titlecolor = ρ < -0.5 ? "#9a2626" : "#52514e", grid...,
        )
        phase = Axis(fig[i, 2]; ylabel = "phase (deg)", xlabel = i == length(sources) ? "frequency (GHz)" : "", grid...)
        for (j, q) in pairs(parallel)
            f = feeds[q][1]
            color = FEED_COLORS[mod1(f, length(FEED_COLORS))]
            lines!(amp, ν, 100 .* (abs.(m[src][j]) .- 1); color, linewidth = 2, label = "feed $f")
            lines!(phase, ν, rad2deg.(angle.(m[src][j])); color, linewidth = 2, label = "feed $f")
        end
        i == 1 && axislegend(amp; position = :lt, framevisible = false, labelsize = 11)
    end
    finding = isempty(opposite) ? "no source's residual is the opposite of the others'" :
        join(opposite, ", ") * ": residual opposite to the other sources' (a spectrum of its own; see apply_only_sources)"
    Label(fig[-1, 1:2], "mean corrected residual spectrum per source"; fontsize = 15)
    Label(fig[0, 1:2], finding; fontsize = 12, color = isempty(opposite) ? "#52514e" : "#9a2626")
    return fig
end

function _correlation(a, b)
    ok = [c for c in eachindex(a, b) if isfinite(a[c]) && isfinite(b[c])]
    length(ok) < 3 && return NaN
    x, y = a[ok] .- mean(a[ok]), b[ok] .- mean(b[ok])
    return sum(x .* y) / sqrt(sum(abs2, x) * sum(abs2, y))
end

"""
    plot_channel_edges(edges) -> Figure

The [`channel_edges`](@ref) of every baseline: how much the phase step across
each channel edge changes from scan to scan, over the baseline's median edge.
"""
function plot_channel_edges(edges)
    shared = shared_edges(edges)
    finding = isempty(shared) ? "no station's baselines share a channel edge" :
        join(("$(x.station): edges " * join(("$(c)|$(c + 1)" for c in x.edges), ", ") * " on $(x.baselines) baselines" for x in shared), "; ") *
        " (see scan_terms)"
    rows_ = [row for row in axes(edges.ratio, 1) if any(isfinite, view(edges.ratio, row, :))]
    labels = [join(edges.baselines[row], "-") for row in rows_]
    fig = Figure(size = (900, 140 + 16 * length(rows_)), backgroundcolor = SURFACE)
    ax = Axis(
        fig[1, 1]; backgroundcolor = SURFACE, yticks = (eachindex(rows_), labels), yticklabelsize = 9,
        xticks = 1:2:size(edges.ratio, 2), xlabel = "channel edge c (between channels c and c + 1)",
        title = "scan-to-scan spread of the phase step at each channel edge, over the baseline's median edge", titlesize = 12,
    )
    z = permutedims(log10.(clamp.(edges.ratio[rows_, :], 1, 30)))
    hm = heatmap!(ax, axes(edges.ratio, 2), eachindex(rows_), z; colormap = SEQUENTIAL, colorrange = (0, log10(30)), nan_color = SURFACE)
    ticks = [1, 2, 4, 10, 30]
    Colorbar(fig[1, 2], hm; label = "spread ÷ median edge", ticks = (log10.(ticks), string.(ticks)))
    Label(fig[0, 1:2], finding; fontsize = 12, color = isempty(shared) ? "#52514e" : "#9a2626")
    return fig
end

# The baselines the before/after plots show: each watched station's busiest
# `partners`, then from the [`baseline_summary`](@ref) `rows` the `worst` fits
# after the correction and the `largest` corrections.
function plotted_baselines(track, watch, rows; partners = 2, worst = 4, largest = 2)
    selected = Tuple{String, String}[]
    counts = baseline_scan_counts(track.ps)
    stations = Gustavo.Calibration.DataGeometry(track.ps).antennas
    for st in watch
        st in stations || (@warn "watch station absent from the track" station = st; continue)
        mine = filter(p -> st in first(p), counts)
        append!(selected, first.(first(mine, partners)))
    end
    for (ranked, n) in ((sort(rows; by = r -> -r.after), worst), (sort(rows; by = r -> -r.correction), largest))
        added = 0
        for r in ranked
            added >= n && break
            r.baseline in selected || (push!(selected, r.baseline); added += 1)
        end
    end
    return selected
end

"""
    correction_sizes(sol, baselines) -> Dict

How much the bandpass shape of `sol` changes each of `baselines` (station-name
pairs): the rms over channels of the phase of `g_a conj(g_b)` about its band
mean, in radians, the largest over the solution's time samples and feeds.
"""
function correction_sizes(sol, baselines)
    G = gains(shape_only(sol))
    g = parent(permutedims(G, (Frequency, Ti, AntennaName, Feed)))
    index = Dict(String(st) => a for (a, st) in pairs(lookup(G, AntennaName)))
    out = Dict{Tuple{String, String}, Float64}()
    for bl in baselines
        a, b = get(index, bl[1], nothing), get(index, bl[2], nothing)
        (isnothing(a) || isnothing(b)) && (out[bl] = 0.0; continue)
        out[bl] = maximum(Iterators.product(axes(g, 2), axes(g, 4))) do (t, f)
            z = view(g, :, t, a, f) .* conj.(view(g, :, t, b, f))
            ok = filter(isfinite, z)
            isempty(ok) ? 0.0 : rms(angle.(ok .* conj(sum(ok))))
        end
    end
    return out
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
change `dp` (rad) and log-amplitude change `da`; per feed the `span` of the
phase change (rad, highest channel less lowest), the `delay` (s) it carries and
the `ramp` that delay makes across the band (rad), its `shape` (rad), the rms
of the phase change beyond that delay and its band mean, and `scatter`, the rms
of `da` about its band mean. Channels with `NaN` changes are left out of the
per-feed values.
"""
struct BreakChange{T}
    dp::Matrix{T}
    da::Matrix{T}
    span::Vector{T}
    delay::Vector{T}
    ramp::Vector{T}
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
    field(k) = [getfield(x, k) for x in stats]
    return BreakChange(collect(dp), collect(da), field(:span), field(:τ), field(:ramp), field(:shape), field(:scatter))
end

function _feed_change(dp, da, ν)
    ok = [c for c in eachindex(dp, da, ν) if isfinite(dp[c]) && isfinite(da[c])]
    length(ok) < 2 && return (; span = NaN, τ = NaN, ramp = NaN, shape = NaN, scatter = NaN)
    z = cis.(dp[ok])
    τ = best_delay(z, ones(length(ok)), ν[ok])
    beyond = z .* cis.(-2π * τ .* ν[ok])
    return (;
        span = maximum(dp[ok]) - minimum(dp[ok]), τ, ramp = 2π * abs(τ) * (maximum(ν[ok]) - minimum(ν[ok])),
        shape = rms(angle.(beyond .* conj(sum(beyond)))), scatter = std(da[ok]),
    )
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
        println(io, "Per station and feed: the span of the phase change (highest channel less lowest),")
        println(io, "the delay it carries and the ramp that delay makes across the band, the rms of the")
        println(io, "phase change beyond that delay and its band mean (shape), the rms of the")
        println(io, "log-amplitude change about its band mean, and for each feed after the first the")
        println(io, "change of its band-mean phase from the first feed's.\n")
        for (k, h) in enumerate(boundaries)
            println(io, "$(label(k))  →  $(label(k + 1))")
            println(io, "station feed  span(deg)  delay(ns)  ramp(deg)  shape(deg)  dlogamp(%)  feed-phase(deg)")
            for (st, c) in sort!(collect(break_changes(s.sol, h, origin)); by = first), f in axes(c.dp, 2)
                feed = f == 1 ? 0.0 : angle(sum(cis.(view(c.dp, :, f))) * conj(sum(cis.(view(c.dp, :, 1)))))
                @printf(
                    io, "%-7s %4d  %9.1f  %9.3f  %9.1f  %10.2f  %10.2f  %15.1f\n", st, f, rad2deg(c.span[f]), 1.0e9 * c.delay[f],
                    rad2deg(c.ramp[f]), rad2deg(c.shape[f]), 100 * c.scatter[f], rad2deg(feed)
                )
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
