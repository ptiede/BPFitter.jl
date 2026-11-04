const PYTHON_PROJECT = joinpath(pkgdir(@__MODULE__), "python")

"""
    simulate_antab(cfg::TrackConfig, sefd::SimulatedSEFD; uv = "uv") -> String

Write `simulated_antab(cfg)`: an ANTAB file whose system temperatures are
ngehtsim's a-priori SEFDs per station, scan and IF, with a DPFU of 1 and a flat
gain curve, so `Gustavo.apriori_calibrate!` reads each value as the SEFD. The
SEFDs come from ngehtsim's weather, receiver and dish tables for the track's
dates, elevations and frequencies (`python/simulate_antab.py`, run with `uv`).
"""
function simulate_antab(cfg::TrackConfig, sefd::SimulatedSEFD; uv = "uv")
    out = simulated_antab(cfg)
    mkpath(dirname(out))
    paths = track_source_paths(cfg.data; skip = cfg.skip_sources)
    stations = join(("$k=$v" for (k, v) in sort(collect(sefd.stations))), ",")
    script = joinpath(PYTHON_PROJECT, "simulate_antab.py")
    run(`$uv run --project $PYTHON_PROJECT python $script --out $out --array $(sefd.array) --stations $stations --weather $(sefd.weather) $paths`)
    return out
end
