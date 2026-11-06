using BPFitter
using BPFitter: scan_gaps, check_gap, overlapping_sources, add_break, drop_break, _round_inside, track_hours
using Gustavo
using XRadio
using DimensionalData
using Random: Xoshiro
using Test

include(joinpath(pkgdir(Gustavo), "test", "synthetic_uvfits.jl"))

const T0 = 1.7125344e9

# Three stations near the north pole, so a source at +80° declination stays at
# about 80° elevation; four one-channel IFs, `nscans` scans ten minutes apart of
# three integrations each.
function write_track(dir; source = "SRC", vis = (r, p, c) -> complex(1.0, 0.0), nscans = 2)
    starts = T0 .+ 600.0 .* (0:(nscans - 1))
    times = Float64[]
    baselines = Tuple{Int, Int}[]
    for t0 in starts, t in t0 .+ (0, 10, 20), bl in ((1, 2), (1, 3), (2, 3))
        push!(times, t)
        push!(baselines, bl)
    end
    nrec = length(times)
    write_synthetic_uvfits(
        joinpath(dir, "hops_TEST_$(source).uvfits");
        times, baselines, uvw = zeros(nrec, 3),
        vis = ComplexF32[vis(r, p, c) for r in 1:nrec, p in 1:2, c in 1:4],
        weights = ones(Float32, nrec, 2, 4),
        ref_freq = 228.0e9, if_freqs = [0.0, 58.0e6, 116.0e6, 174.0e6], ch_width = fill(58.0e6, 4),
        antennas = ["AA", "BB", "CC"],
        positions = [0.0 1.0e5 -1.0e5; 1.0e5 0.0 -1.0e5; 6.356e6 6.356e6 6.356e6],
        scans = [(t0 - 5, t0 + 25) for t0 in starts],
        object = source, ra = 187.7, dec = 80.0, rdate = "2024-04-08",
    )
    return dir
end

# `sefd`: station => one SEFD, or one per scan. One TSYS row per scan, at the
# scan's start (scans ten minutes apart from midnight of day 99).
function write_antab(path, sefd; nscans = 2, nif = 4)
    open(path, "w") do io
        for st in sort(collect(keys(sefd)))
            println(io, "GAIN $st ELEV DPFU = 1.0, 1.0 POLY = 1.0 /")
        end
        index = join(("'L$i|R$i'" for i in 1:nif), ", ")
        for st in sort(collect(keys(sefd)))
            println(io, "TSYS $st FT = 1.0 TIMEOFF = 0 INDEX = $index /")
            for k in 1:nscans
                v = join(fill(string(sefd[st] isa Number ? sefd[st] : sefd[st][k]), nif), " ")
                h, m = divrem(10 * (k - 1), 60)
                println(io, " 99 $(lpad(h, 2, '0')):$(lpad(m, 2, '0')).2 $v")
            end
            println(io, "/")
        end
    end
    return path
end

# The bandpass of `STEP_STATIONS` across eight one-channel IFs: phase (rad) and
# log-amplitude per channel, and the change of BB's after scan `STEP_AFTER`.
const STEP_STATIONS = ["AA", "BB", "CC", "DD"]
const STEP_AFTER = 6
step_x(c) = (c - 4.5) / 3.5
step_phase(a, c) = (0.0, 0.3 * step_x(c)^2 - 0.1, 0.15 * cos(2 * step_x(c)), 0.25 * sin(3 * step_x(c)))[a]
step_logamp(a, c) = (0.05 * step_x(c)^2, -0.08 * step_x(c), 0.04 * cos(2 * step_x(c)), 0.0)[a]
step_dphase(c) = 0.35 * cos(2.5 * step_x(c))
step_dlogamp(c) = 0.1 * step_x(c)^2 - 0.03

# Four stations near the north pole observing one source at +80°, twelve scans
# of three integrations, with the bandpass above and noise of `σ` per
# visibility component (weights 1/σ²). The source is 1 Jy with SEFDs of 1 Jy.
function write_step_track(dir; σ = 0.05, nscans = 12)
    rng = Xoshiro(7)
    starts = T0 .+ 600.0 .* (0:(nscans - 1))
    times = Float64[]
    baselines = Tuple{Int, Int}[]
    scan_of = Int[]
    for (k, t0) in pairs(starts), t in t0 .+ (0, 10, 20), a in 1:4, b in (a + 1):4
        push!(times, t)
        push!(baselines, (a, b))
        push!(scan_of, k)
    end
    gain(a, c, k) = exp(
        complex(
            step_logamp(a, c) + (a == 2 && k > STEP_AFTER ? step_dlogamp(c) : 0.0),
            step_phase(a, c) + (a == 2 && k > STEP_AFTER ? step_dphase(c) : 0.0),
        )
    )
    nrec = length(times)
    vis = ComplexF32[
        gain(baselines[r][1], c, scan_of[r]) * conj(gain(baselines[r][2], c, scan_of[r])) + σ * complex(randn(rng), randn(rng))
            for r in 1:nrec, p in 1:2, c in 1:8
    ]
    write_synthetic_uvfits(
        joinpath(dir, "hops_TEST_SRC.uvfits");
        times, baselines, uvw = zeros(nrec, 3), vis,
        weights = fill(Float32(1 / σ^2), nrec, 2, 8),
        ref_freq = 228.0e9, if_freqs = 58.0e6 .* (0:7), ch_width = fill(58.0e6, 8),
        antennas = STEP_STATIONS,
        positions = [0.0 1.0e5 -1.0e5 0.0; 1.0e5 0.0 -1.0e5 -1.0e5; 6.356e6 6.356e6 6.356e6 6.356e6],
        scans = [(t0 - 5, t0 + 25) for t0 in starts],
        object = "SRC", ra = 187.7, dec = 80.0, rdate = "2024-04-08",
    )
    return dir
end

function write_campaign(dir, body)
    path = joinpath(dir, "campaign.toml")
    write(path, body)
    return path
end

@testset "BPFitter" begin
    @testset "campaign files" begin
        dir = mktempdir()
        path = write_campaign(
            dir, """
            [defaults]
            data = "data/{track}"
            antab = "antab/{label}.AN"
            output = "out/{track}"
            station_flags = { ST = [1] }
            feed_flags = { SG = { 1 = [30, 31] } }
            [breaks]
            alpha = 1e-4
            [sefd]
            array = "ARRAY"
            stations = { AA = "ALMA" }
            [tracks.A]
            label = "a1"
            breaks = { AA = [12.25, 3.0] }
            [tracks.B]
            station_flags = { ST = [1, 7, 8] }
            scan_flags = [{ stations = ["ST", "SP"], scan = "X_1" }]
            apply_only_sources = ["VARIABLE"]
            """
        )
        c = read_campaign(path)
        @test BPFitter.track_names(c) == ["A", "B"]
        a, b = c["A"], c["B"]
        @test a.data == joinpath(dir, "data", "A")
        @test a.antab == joinpath(dir, "antab", "a1.AN")
        @test b.antab == joinpath(dir, "antab", "B.AN")
        @test a.breaks == Dict("AA" => [3.0, 12.25])
        @test isempty(b.breaks)
        @test a.station_flags == Dict("ST" => [1])
        @test b.station_flags == Dict("ST" => [1, 7, 8])
        @test a.feed_flags == Dict("SG" => Dict(1 => [30, 31]))
        @test b.scan_flags == [(["ST", "SP"], "X_1")]
        @test b.apply_only_sources == ["VARIABLE"]
        @test c.search.alpha == 1.0e-4
        @test c.search.min_F == 40.0 && c.search.min_shape_deg == 2.0 && c.search.min_scatter_percent == 5.0 && c.search.tries == 4
        @test c.sefd.stations == Dict("AA" => "ALMA")
        @test c.sefd.array == "ARRAY"
        @test_throws "no track \"C\"" c["C"]

        bad = write_campaign(mktempdir(), "[defaults]\nouptut = \"x\"\n[tracks.A]\ndata = \"d\"\n")
        @test_throws "unknown key(s) [\"ouptut\"]" read_campaign(bad)
        missing_output = write_campaign(mktempdir(), "[tracks.A]\ndata = \"d\"\n")
        @test_throws "track A has no `output`" read_campaign(missing_output)

        mkpath(joinpath(dir, "out", "B"))
        write(joinpath(dir, "out", "B", "breaks.toml"), "[breaks]\nST = [20.5]\n")
        @test isempty(read_campaign(path)["B"].breaks)
        @test read_campaign(path)["A"].breaks == Dict("AA" => [3.0, 12.25])
        @test BPFitter.with_breaks(a, Dict("SP" => [1.0])).breaks == Dict("SP" => [1.0])
    end

    @testset "timeline" begin
        rows = [
            (; source = "A", scan = "A_1", t0 = 1.0, t1 = 1.2, stations = ["X", "Y"]),
            (; source = "B", scan = "B_1", t0 = 1.1, t1 = 1.3, stations = ["Z"]),
            (; source = "A", scan = "A_2", t0 = 2.0, t1 = 2.2, stations = ["X", "Y"]),
        ]
        @test scan_gaps(rows) == [(; t0 = 1.3, t1 = 2.0, hours = 2.0 - 1.3)]
        @test isnothing(check_gap(rows, 1.5))
        @test_throws "falls inside scans" check_gap(rows, 1.25)
        @test overlapping_sources(rows) == [("A", "B")]
        later = vcat(rows[3:3], [(; source = "C", scan = "C_1", t0 = 2.1, t1 = 2.3, stations = ["Z"])])
        @test overlapping_sources(vcat(rows[1:1], later)) == [("A", "C")]
    end

    @testset "break bookkeeping" begin
        b = Dict("AA" => [12.25])
        @test add_break(b, "AA", 3.0) == Dict("AA" => [3.0, 12.25])
        @test add_break(b, "ST", 1.0) == Dict("AA" => [12.25], "ST" => [1.0])
        @test drop_break(Dict("AA" => [3.0, 12.25], "ST" => [1.0]), "ST", 1.0) == Dict("AA" => [3.0, 12.25])
        @test drop_break(Dict("AA" => [3.0, 12.25]), "AA", 3.0) == Dict("AA" => [12.25])
        @test _round_inside(1.23456, 1.2, 1.3) == 1.23
        @test _round_inside(1.2049, 1.2, 1.2051) == 1.2049
    end

    @testset "loading and a-priori calibration" begin
        dir = mktempdir()
        write_track(mkpath(joinpath(dir, "data", "T")); source = "SRC")
        write_track(joinpath(dir, "data", "T"); source = "SKIPPED")
        cp(joinpath(dir, "data", "T", "hops_TEST_SRC.uvfits"), joinpath(dir, "data", "T", "hops_TEST_SRC+avg.uvfits"))
        sefd = Dict("AA" => 100.0, "BB" => 400.0, "CC" => 900.0)
        write_antab(mkpath(joinpath(dir, "antab")) * "/T.AN", sefd)
        path = write_campaign(
            dir, """
            [defaults]
            data = "data/{track}"
            antab = "antab/{track}.AN"
            output = "out/{track}"
            skip_sources = ["SKIPPED"]
            station_flags = { CC = [2] }
            [tracks.T]
            """
        )
        cfg = read_campaign(path)["T"]
        @test basename.(BPFitter.track_source_paths(cfg.data; skip = cfg.skip_sources)) == ["hops_TEST_SRC.uvfits"]
        track = prepare_track(cfg)
        @test BPFitter.track_sources(track.ps) == ["SRC"]
        @test sort(unique(r.scan for r in track.rows)) == ["SRC_1", "SRC_2"]
        @test track.rows[1].t0 ≈ track_hours(T0, track.origin)
        for ms in values(track.ps)
            V, F = ms[:visibility], XRadio.flags(ms)
            for (bi, (a, b)) in pairs(XRadio.baselines(ms))
                a == b && continue
                expected = sqrt(sefd[String(a)] * sefd[String(b)])
                cells = V[BaselineID(bi)]
                flagged = F[BaselineID(bi)]
                hit_cc = "CC" in (String(a), String(b))
                @test all(isapprox.(abs.(cells[.!flagged]), expected; rtol = 1.0e-5))
                @test all(flagged[Frequency = 2]) == hit_cc
            end
        end
    end

    @testset "SEFD outliers" begin
        dir = mktempdir()
        write_track(dir; nscans = 6)
        ps = BPFitter.load_track(dir)
        sefd = Dict("AA" => [100.0, 110.0, 105.0, 2000.0, 95.0, 100.0], "BB" => fill(400.0, 6), "CC" => [900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0])
        XRadio.read_antab!(ps, write_antab(joinpath(dir, "T.AN"), sefd; nscans = 6); ifs = BPFitter.channel_ifs(ps))
        odd = sefd_outliers(ps)
        @test [(r.station, r.scan) for r in odd] == [("AA", "SRC_4")]
        @test only(odd).sefd ≈ 2000.0
        @test only(odd).typical ≈ 102.5
        @test isempty(sefd_outliers(ps; factor = 25.0))
    end

    @testset "breaks, solve and apply on a synthetic track" begin
        dir = mktempdir()
        write_step_track(mkpath(joinpath(dir, "data", "T")))
        write_antab(mkpath(joinpath(dir, "antab")) * "/T.AN", Dict(st => 1.0 for st in STEP_STATIONS); nscans = 12, nif = 8)
        campaign_body(breaks) = """
        [defaults]
        data = "data/{track}"
        antab = "antab/{track}.AN"
        output = "out/{track}"
        corrected = "out/corrected/{track}"
        ref_ant = "AA"
        [tracks.T]
        $breaks
        """
        cfg = read_campaign(write_campaign(dir, campaign_body("")))["T"]
        track = prepare_track(cfg)
        gap = (track.rows[STEP_AFTER].t1, track.rows[STEP_AFTER + 1].t0)

        r = find_breaks(track; io = devnull)
        @test collect(keys(r.breaks)) == ["BB"]
        h = only(r.breaks["BB"])
        @test gap[1] < h < gap[2]
        @test all(BPFitter.passed, r.final)
        # A change below both sizes is neither refit nor suggested.
        @test isempty(find_breaks(track; search = BreakSearch(; min_shape_deg = 90, min_scatter_percent = 90), io = devnull).breaks)

        suggestions = write_suggestions(joinpath(cfg.output, "break_suggestions"), track, r)
        table = read(joinpath(suggestions, "suggestions.txt"), String)
        @test occursin("breaks = { BB = [$h] }", table)
        @test isfile(joinpath(suggestions, "suggested_BB_$h.png"))

        # The change across the break matches the injected one.
        ν = track.ps |> values |> first |> XRadio.frequencies
        ν = ν .- sum(ν) / length(ν)
        truth = BPFitter._station_change(hcat(step_dphase.(1:8), step_dphase.(1:8)), hcat(step_dlogamp.(1:8), step_dlogamp.(1:8)), ν)
        found = r.changes[("BB", h)]
        @test found.shape ≈ truth.shape rtol = 0.1
        @test found.scatter ≈ truth.scatter rtol = 0.15

        # With the break the corrected spectra are flat to the noise; without it they are not.
        χ²(s) = (q = BPFitter.scan_chi2(s.chi2); q.chi2 / q.dof)
        @test 0.8 < χ²(solve_track(track; breaks = r.breaks)) < 1.25
        @test χ²(solve_track(track; breaks = Dict{String, Vector{Float64}}())) > 3

        # `run` solves with the campaign file's breaks, writes the reports and exports;
        # `apply` writes the same files from the saved solution.
        campaign = read_campaign(write_campaign(dir, campaign_body("breaks = { BB = [$h] }")))
        redirect_stdout(devnull) do
            BPFitter.run_command("run", campaign, ["T"])
        end
        for f in ("chi2.txt", "timeline.txt", "break_delta.txt", "solution_phase.png", "bandpass.zarr", "scan_delays.zarr")
            @test ispath(joinpath(cfg.output, f))
        end
        # The corrected track fits to the noise: one source, no scan or channel edge stands out.
        chi2 = read(joinpath(cfg.output, "chi2.txt"), String)
        @test occursin(r"SRC +12 ", chi2)
        @test occursin(r"scans that fit worst[^\n]*\n    none", chi2)
        @test occursin(r"channel edges[^\n]*\n    none", chi2)
        exported = joinpath(cfg.corrected, "hops_T_SRC_bpcal.uvfits")
        @test 0.8 < BPFitter.scan_chi2(BPFitter.track_tss(load_uvfits(exported))) |> (q -> q.chi2 / q.dof) < 1.25
        rm(exported)
        redirect_stdout(devnull) do
            BPFitter.run_command("apply", campaign, ["T"])
        end
        @test 0.8 < BPFitter.scan_chi2(BPFitter.track_tss(load_uvfits(exported))) |> (q -> q.chi2 / q.dof) < 1.25
    end

    @testset "dead channels" begin
        dir = mktempdir()
        write_track(dir; vis = (r, p, c) -> c == 3 && r % 3 != 0 ? complex(0.01, 0.0) : complex(1.0, 0.0))
        ps = BPFitter.load_track(dir)
        dead = dead_channels(ps)
        # Records 1, 2 of every 3 are AA–BB and AA–CC: only AA is on all of the low ones.
        @test [(r.station, r.channel) for r in dead] == [("AA", 3), ("AA", 3)]
        @test sort(unique(r.feed for r in dead)) == [1, 2]
    end
end
