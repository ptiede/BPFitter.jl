const USAGE = """
usage: bpfit <command> <campaign.toml> [track ...]

Every command runs on the named tracks, or on all the campaign's tracks.

  inspect   sources, stations, channels, timeline, sub-array overlaps and dead
            channels of each track, read from its raw uvfits
  sefd      simulate an ANTAB file with ngehtsim for each track that has none
  breaks    suggest bandpass breaks for review; writes a table and plots to
            <output>/break_suggestions (adopt them in the campaign file)
  solve     solve with the campaign file's breaks; writes solutions, tables
            and plots
  apply     write corrected uvfits from the solutions `solve` saved
  plot      redraw the solution plots from the solutions `solve` saved:
            solution_phase.png and solution_amp.png with the campaign file's
            breaks, and scan_delays.png
  run       sefd where needed, solve, apply
"""

"""
    apply_saved(track) -> TrackSolve

The solve `solve` saved under the track's output folder, applied to the track:
its per-scan delays and per-integration phases and amplitudes, then its
bandpass in full to the scans in the solve and its shape alone to the rest.
"""
function apply_saved(track::Track)
    dir = track.cfg.output
    sol = load_solution(joinpath(dir, "bandpass.zarr"))
    delays_path = joinpath(dir, "scan_delays.zarr")
    delays = isdir(delays_path) ? load_solution(delays_path) : nothing
    adhoc_path = joinpath(dir, "adhoc.zarr")
    adhoc = isdir(adhoc_path) ? load_solution(adhoc_path) : nothing
    corrected = isnothing(delays) ? deepcopy(track.ps) : calibrate(delays, track.ps)
    isnothing(adhoc) || calibrate!(adhoc, corrected)
    keep = in_solve(track.cfg)
    calibrate!(sol, filter(keep, corrected))
    rest = filter(!keep, corrected)
    isempty(rest) || calibrate!(shape_only(sol), rest)
    model = track_bandpass_model(track.cfg.breaks, track.origin, track.cfg.priors)
    return TrackSolve(model, track.cfg.breaks, something(track.cfg.ref_ant, ""), sol, delays, adhoc, corrected, track_tss(filter(keep, corrected)))
end

"""
    plot_saved(track) -> Vector{String}

The solution plots of the solve `solve` saved under the track's output folder:
`solution_phase.png` and `solution_amp.png`, with the track's configured
breaks, and `scan_delays.png` when per-scan delays were saved. Returns the paths
written.
"""
function plot_saved(track::Track)
    dir = track.cfg.output
    bandpass = joinpath(dir, "bandpass.zarr")
    isdir(bandpass) || error("no saved bandpass at $bandpass; run `bpfit solve` first")
    out = String[]
    write_solution_plots(dir, shape_only(load_solution(bandpass)), track.cfg.breaks, track.origin)
    append!(out, joinpath.(dir, ("solution_phase.png", "solution_amp.png")))
    delays = joinpath(dir, "scan_delays.zarr")
    if isdir(delays)
        push!(out, joinpath(dir, "scan_delays.png"))
        save(last(out), plot_scan_delays(load_solution(delays), track.rows, track.origin))
    end
    return out
end

function solve_and_report(track; breaks = track.cfg.breaks)
    s = solve_track(track; breaks)
    show(stdout, MIME"text/plain"(), s.sol)
    println()
    save_solve(track.cfg.output, s)
    write_reports(track, s)
    print(read(joinpath(track.cfg.output, "chi2.txt"), String))
    return s
end

function run_command(command, campaign::Campaign, names)
    for name in names
        cfg = campaign[name]
        @info "track $name: $command"
        if command == "inspect"
            inspect_track(stdout, cfg)
        elseif command == "sefd"
            cfg.apriori || (@info "track $name has `apriori = false`; nothing to simulate"; continue)
            isnothing(cfg.antab) || (@info "track $name has an ANTAB file; nothing to simulate"; continue)
            isnothing(campaign.sefd) && error("$(campaign.path) has no [sefd] table to simulate SEFDs with")
            simulate_antab(cfg, campaign.sefd)
        elseif command == "breaks"
            track = prepare_track(cfg)
            write_suggestions(joinpath(cfg.output, "break_suggestions"), track, find_breaks(track; search = campaign.search))
        elseif command == "solve"
            solve_and_report(prepare_track(cfg))
        elseif command == "apply"
            export_track(cfg, apply_saved(prepare_track(cfg)))
        elseif command == "plot"
            foreach(out -> @info("wrote", out), plot_saved(prepare_track(cfg)))
        elseif command == "run"
            if cfg.apriori && isnothing(cfg.antab) && !isfile(simulated_antab(cfg))
                isnothing(campaign.sefd) && error("track $name has no ANTAB file and $(campaign.path) has no [sefd] table")
                simulate_antab(cfg, campaign.sefd)
            end
            s = solve_and_report(prepare_track(cfg))
            isnothing(cfg.corrected) || export_track(cfg, s)
        else
            error("unknown command $(repr(command))\n\n$USAGE")
        end
    end
    return nothing
end

function (@main)(args)
    if length(args) < 2 || args[1] in ("-h", "--help", "help")
        print(USAGE)
        return length(args) < 2 && !(get(args, 1, "") in ("-h", "--help", "help")) ? 1 : 0
    end
    BLAS.set_num_threads(1)
    campaign = read_campaign(args[2])
    names = length(args) > 2 ? args[3:end] : track_names(campaign)
    run_command(args[1], campaign, names)
    return 0
end
