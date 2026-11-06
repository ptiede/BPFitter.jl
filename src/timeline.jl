const SECONDS_PER_HOUR = 3600.0

rms(x) = sqrt(mean(abs2, x))

"""
    track_origin(ps) -> Float64

The epoch of UT midnight on the day of the track's first sample, in the
absolute seconds the time axis carries. Every hour in a configuration, a
timeline or a break is measured from here, so a track running through UT
midnight runs past 24.
"""
function track_origin(ps)
    t0 = minimum(minimum(XRadio.times(ms)) for ms in values(ps))
    return datetime2unix(DateTime(Date(unix2datetime(t0))))
end

"Absolute epoch(s) → hours since `origin`."
track_hours(t, origin) = (t .- origin) ./ SECONDS_PER_HOUR

"Hours since `origin` → absolute epoch(s)."
track_epoch(h, origin) = origin .+ SECONDS_PER_HOUR .* h

"The cross baselines of `ms`, as station-name pairs in stored order."
cross_baselines(ms) = [(String(a), String(b)) for (a, b) in XRadio.baselines(ms) if a != b]

"The stations carrying a cross baseline in `ms`."
ms_stations(ms) = sort!(unique(reduce(vcat, [[a, b] for (a, b) in cross_baselines(ms)]; init = String[])))

"""
    scan_rows(ps, origin) -> Vector{<:NamedTuple}

One row per scan, in time order: the source, the scan name, the window it spans
in hours since `origin`, and the stations on it.
"""
function scan_rows(ps, origin)
    rows = map(collect(values(groupby(ps, ByScan())))) do g
        ms = first(values(g))
        t = XRadio.times(ms)
        (;
            source = source_name(ms), scan = String(only(unique(collect(ms[:scan_name])))),
            t0 = track_hours(minimum(t), origin), t1 = track_hours(maximum(t), origin),
            stations = ms_stations(ms),
        )
    end
    return sort!(rows; by = r -> r.t0)
end

"""
    scan_gaps(rows; min_hours = 0.0) -> Vector{<:NamedTuple}

The intervals of at least `min_hours` in which the track observed nothing,
between the end of everything started so far and the start of the next scan.
"""
function scan_gaps(rows; min_hours = 0.0)
    gaps = map(r -> (; t0 = r.t0, t1 = r.t1, hours = r.t1 - r.t0), empty(rows))
    tmax = -Inf
    for r in rows
        if isfinite(tmax) && r.t0 - tmax >= min_hours
            push!(gaps, (; t0 = tmax, t1 = r.t0, hours = r.t0 - tmax))
        end
        tmax = max(tmax, r.t1)
    end
    return gaps
end

"""
    check_gap(rows, hour)

Raise unless `hour` falls between scans rather than inside one: a boundary
through a scan would split its data across two bandpasses.
"""
function check_gap(rows, hour)
    straddling = [(r.source, r.scan, r.t0, r.t1) for r in rows if r.t0 < hour < r.t1]
    isempty(straddling) ||
        error("boundary $hour h falls inside scans $straddling; move it into a gap")
    return nothing
end

"""
    overlapping_sources(rows) -> Vector{Tuple{String, String}}

Pairs of sources observed at the same time, by two sub-arrays. One time axis
cannot rank a sample into two scans, so one of each pair has to be skipped
(`skip_sources`) for the joint solve.
"""
function overlapping_sources(rows)
    found = Set{Tuple{String, String}}()
    for (i, a) in pairs(rows)
        for b in rows[(i + 1):end]
            b.t0 >= a.t1 && break
            a.source == b.source || push!(found, minmax(a.source, b.source))
        end
    end
    return sort!(collect(found))
end

"""
    station_scan_counts(rows) -> Vector{Pair{String, Int}}

How many scans each station appears on, most first.
"""
function station_scan_counts(rows)
    counts = Dict{String, Int}()
    for r in rows, st in r.stations
        counts[st] = get(counts, st, 0) + 1
    end
    return sort!(collect(counts); by = last, rev = true)
end

"""
    baseline_scan_counts(ps) -> Vector{Pair{Tuple{String, String}, Int}}

How many scans each baseline, in stored order, appears on, most first.
"""
function baseline_scan_counts(ps)
    counts = Dict{Tuple{String, String}, Int}()
    for g in values(groupby(ps, ByScan())), bl in cross_baselines(first(values(g)))
        counts[bl] = get(counts, bl, 0) + 1
    end
    return sort!(collect(counts); by = last, rev = true)
end

"""
    write_timeline(io, rows, boundaries; gap_min_hours = 0.1)

The scan table with the segment each scan falls in under `boundaries`, followed
by the gaps a break could be placed in and each station's coverage.
"""
function write_timeline(io::IO, rows, boundaries; gap_min_hours = 0.1)
    seg_of(t) = searchsortedlast(sort(collect(boundaries)), t) + 1
    println(io, "source          scan    t0(h)   t1(h)  seg  stations")
    for r in rows
        @printf(
            io, "%-14s %6s  %7.3f %7.3f  %3d  %s\n",
            r.source, r.scan, r.t0, r.t1, seg_of(r.t0), join(r.stations, " ")
        )
    end
    gaps = scan_gaps(rows; min_hours = gap_min_hours)
    @printf(io, "\ngaps of at least %.2f h (a break belongs in one)\n", gap_min_hours)
    for g in sort(gaps; by = g -> g.hours, rev = true)
        @printf(io, "  %7.3f – %7.3f   %5.2f h\n", g.t0, g.t1, g.hours)
    end
    println(io, "\nstation coverage")
    @printf(io, "  %-6s %6s %8s %8s\n", "st", "scans", "first(h)", "last(h)")
    for (st, n) in station_scan_counts(rows)
        ts = [r.t0 for r in rows if st in r.stations]
        @printf(io, "  %-6s %6d %8.3f %8.3f\n", st, n, minimum(ts), maximum(ts))
    end
    return nothing
end
write_timeline(path::AbstractString, args...; kw...) = open(io -> write_timeline(io, args...; kw...), path, "w")

"""
    strip_common(d; keep) -> Array

`d` with the across-station median removed channel by channel and feed by feed,
over the second axis. That median is the part every station shares, which is
what a gauge difference between two independently solved segments looks like;
what is left belongs to individual stations. `keep` is a (station, feed) mask of
the cells the median is taken over.
"""
function strip_common(d; keep = trues(size(d, 2), size(d, 3)))
    axes(keep) == (axes(d, 2), axes(d, 3)) ||
        throw(DimensionMismatch("keep has axes $(axes(keep)), the (station, feed) axes are $((axes(d, 2), axes(d, 3)))"))
    out = similar(d)
    for f in axes(d, 3), c in axes(d, 1)
        idx = findall(view(keep, :, f))
        out[c, :, f] = view(d, c, :, f) .- (isempty(idx) ? zero(eltype(d)) : median(view(d, c, idx, f)))
    end
    return out
end
