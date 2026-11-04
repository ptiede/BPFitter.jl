"""
    BPFitter

Joint bandpass calibration of HOPS-fringed VLBI tracks with Gustavo: every
source of a track constrains one bandpass per station, cut where the data show
it changed. Configure a campaign in a TOML file ([`read_campaign`](@ref)) and
run the `bpfit` app, or call [`prepare_track`](@ref), [`find_breaks`](@ref),
[`solve_track`](@ref) and [`export_track`](@ref) from Julia.
"""
module BPFitter

using Gustavo
using Gustavo.Fring: normalized_residuals, scan_intervals, step_tests, tss
using XRadio
using FITSFiles
using DimensionalData
using DimensionalData: dims, lookup, dimnum
using Distributions: LogNormal, FDist, quantile
using CairoMakie
using Dates: Date, DateTime, datetime2unix, unix2datetime
using LinearAlgebra: BLAS
using Printf: @printf, @sprintf
using Statistics
using TOML

export read_campaign, Campaign, TrackConfig, BreakSearch, BandpassPriors, SimulatedSEFD
export prepare_track, Track, solve_track, TrackSolve, find_breaks, BreakSearchResult, BreakTest, BreakChange, break_changes
export write_suggestions, write_reports, save_solve, export_track, apply_saved, simulate_antab
export inspect_track, dead_channels, sefd_outliers, track_bandpass_model, with_scan_terms, scan_delay_model

include("config.jl")
include("timeline.jl")
include("load.jl")
include("models.jl")
include("solve.jl")
include("report.jl")
include("breaks.jl")
include("inspect.jl")
include("sefd.jl")
include("cli.jl")

end
