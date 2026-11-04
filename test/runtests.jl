using BPFitter
using BPFitter: scan_gaps, check_gap, overlapping_sources, add_break, drop_break, _round_inside, track_hours
using Gustavo
using XRadio
using DimensionalData
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

# `sefd`: station => one SEFD, or one per scan.
function write_antab(path, sefd; nscans = 2)
    open(path, "w") do io
        for st in sort(collect(keys(sefd)))
            println(io, "GAIN $st ELEV DPFU = 1.0, 1.0 POLY = 1.0 /")
        end
        for st in sort(collect(keys(sefd)))
            println(io, "TSYS $st FT = 1.0 TIMEOFF = 0 INDEX = 'L1|R1', 'L2|R2', 'L3|R3', 'L4|R4' /")
            for k in 1:nscans
                v = join(fill(string(sefd[st] isa Number ? sefd[st] : sefd[st][k]), 4), " ")
                println(io, " 99 00:$(lpad(10 * (k - 1), 2, '0')).2 $v")
            end
            println(io, "/")
        end
    end
    return path
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
        @test c.search.min_F == 40.0 && c.search.peaks == 2 && c.search.tries == 4
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
