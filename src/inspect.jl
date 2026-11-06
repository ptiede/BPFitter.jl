"""
    dead_channels(ps; threshold = 0.2) -> Vector{NamedTuple}

Channels a station's feed carries no signal in: on each of the station's cross
baselines and scans, each channel's time-averaged amplitude against the median
channel's, and the median of that ratio over baselines and scans. A
`(station, feed, channel, ratio)` row is returned for every ratio below
`threshold`. The weights of such channels usually look sound, so they need a
`feed_flags` or `station_flags` entry.
"""
function dead_channels(ps; threshold = 0.2)
    ratios = Dict{Tuple{String, Int, Int}, Vector{Float64}}()
    for ms in values(ps)
        vis = ms[:visibility]
        order = dimnum(vis, (Frequency, Ti, BaselineID, Polarization))
        V = PermutedDimsArray(parent(vis), order)
        W = PermutedDimsArray(usable_weights(ms), order)
        fp = feed_pairs(ms)
        for (bi, (a, b)) in pairs(XRadio.baselines(ms)), p in axes(fp, 1)
            a == b && continue
            fa, fb = fp[p, bi]
            fa == fb || continue
            amp = map(axes(V, 1)) do c
                w = view(W, c, :, bi, p)
                sum(w) > 0 ? abs(sum(w .* view(V, c, :, bi, p)) / sum(w)) : NaN
            end
            ok = filter(isfinite, amp)
            length(ok) >= 4 || continue
            ref = median(ok)
            ref > 0 || continue
            for st in (String(a), String(b)), c in eachindex(amp)
                isfinite(amp[c]) && push!(get!(ratios, (st, fa, c), Float64[]), amp[c] / ref)
            end
        end
    end
    rows = [(; station = k[1], feed = k[2], channel = k[3], ratio = median(v)) for (k, v) in ratios]
    return sort!(filter(r -> r.ratio < threshold, rows); by = r -> (r.station, r.feed, r.channel))
end

"""
    sefd_outliers(ps; factor = 3.0, neighbors = 2) -> Vector{NamedTuple}

Scans whose ANTAB SEFD stands apart from the station's other scans of the same
source. `ps` must hold the system temperatures and gain curves
`XRadio.read_antab!` writes. A station's SEFD in a scan it has data in is the
median Tsys/DPFU of the ANTAB rows inside the scan, over receptors and IFs.
Opacity makes the SEFD depend strongly on elevation, so each scan is compared
only with scans of the same source close to it in time, which are at similar
elevations. A `(station, scan, hour, sefd, typical)` row is returned for every
scan more than `factor` from `typical`, the median of the `2·neighbors` such
scans nearest it (`neighbors` on either side away from the ends). Such a row is
usually a failed measurement; flag the scan with `scan_flags`.
"""
function sefd_outliers(ps; factor = 3.0, neighbors = 2)
    A, R = XRadio.AntennaName, XRadio.ReceptorLabel
    origin = track_origin(ps)
    series = Dict{Tuple{String, String}, Vector{@NamedTuple{scan::String, hour::Float64, sefd::Float64}}}()
    for g in values(groupby(ps, ByScan()))
        by_station = Dict{String, Vector{Float64}}()
        for ms in values(g)
            T = XRadio.system_temperatures(ms)
            dpfu = XRadio.gaincurves(ms)[:sensitivity]
            on_scan = ms_stations(ms)
            for a in intersect(lookup(T, A), lookup(dpfu, A)), r in lookup(T, R)
                String(a) in on_scan || continue
                v = vec(parent(T[A(At(a)), R(At(r))])) ./ dpfu[A(At(a)), R(At(r))]
                append!(get!(by_station, String(a), Float64[]), filter(x -> isfinite(x) && x > 0 && x != 999.0, v))
            end
        end
        ms = first(values(g))
        hour = track_hours(mean(XRadio.times(ms)), origin)
        scan, source = String(only(XRadio.scans(ms))), source_name(ms)
        for (st, v) in by_station
            isempty(v) || push!(get!(series, (st, source), valtype(series)()), (; scan, hour, sefd = median(v)))
        end
    end
    rows = @NamedTuple{station::String, scan::String, hour::Float64, sefd::Float64, typical::Float64}[]
    for ((st, _), s) in series
        sort!(s; by = x -> x.hour)
        for i in eachindex(s)
            lo = clamp(i - neighbors, firstindex(s), max(firstindex(s), lastindex(s) - 2 * neighbors))
            near = [s[j].sefd for j in lo:min(lastindex(s), lo + 2 * neighbors) if j != i]
            length(near) >= 2 || continue
            typical = median(near)
            ratio = s[i].sefd / typical
            max(ratio, 1 / ratio) > factor && push!(rows, (; station = st, s[i].scan, s[i].hour, s[i].sefd, typical))
        end
    end
    return sort!(rows; by = r -> (r.station, r.hour))
end

"""
    inspect_track(io, cfg::TrackConfig)

What a campaign file needs to know about a track, read from its raw uvfits:
sources, stations, channels, the scan timeline, sources observed at the same
time by two sub-arrays (one must be skipped), channels that read no signal,
whether its ANTAB file is there, and the scans whose ANTAB SEFDs stand apart
([`sefd_outliers`](@ref)).
"""
function inspect_track(io::IO, cfg::TrackConfig)
    ps = load_track(cfg.data; skip = cfg.skip_sources)
    origin = track_origin(ps)
    rows = scan_rows(ps, origin)
    println(io, "track ", cfg.name, "  (", cfg.data, ")")
    println(io, "sources:  ", join(track_sources(ps), ", "))
    println(io, "stations: ", join(sort!(unique(reduce(vcat, (r.stations for r in rows); init = String[]))), ", "))
    ms = first(values(ps))
    ν = XRadio.frequencies(ms)
    @printf(io, "channels: %d, %.4f–%.4f GHz\n", length(ν), minimum(ν) / 1.0e9, maximum(ν) / 1.0e9)
    antab = antab_path(cfg)
    println(io, "antab:    ", antab, isfile(antab) ? "" : isnothing(cfg.antab) ? "  (missing: run `bpfit sefd`)" : "  (MISSING)")
    overlap = overlapping_sources(rows)
    if !isempty(overlap)
        println(io, "\nsources observed at the same time (skip one of each pair with `skip_sources`):")
        foreach(p -> println(io, "  ", p[1], " / ", p[2]), overlap)
    end
    dead = dead_channels(ps)
    if !isempty(dead)
        println(io, "\nchannels reading under 20% of their baselines' median amplitude:")
        for r in dead
            @printf(io, "  %-4s feed %d  channel %3d  %.2f\n", r.station, r.feed, r.channel, r.ratio)
        end
    end
    if isfile(antab)
        XRadio.read_antab!(ps, antab; ifs = channel_ifs(ps))
        odd = sefd_outliers(ps)
        if !isempty(odd)
            println(io, "\nscans whose SEFD (Tsys/DPFU) is more than 3× above or below that of the station's nearby scans of the same source (look, then flag with `scan_flags`):")
            for r in odd
                @printf(io, "  %-4s %-16s %7.3f h  %10.1f Jy  (neighbors %.1f Jy)\n", r.station, r.scan, r.hour, r.sefd, r.typical)
            end
        end
    end
    println(io)
    write_timeline(io, rows, reduce(vcat, values(cfg.breaks); init = Float64[]))
    return nothing
end
