# The joint solve carries one complex gain per (station, feed, time segment), so
# a station's phase and log-amplitude share a segmentation: a station given a
# break gets it in both observables.
#
# Leaving most stations whole is what makes a break measurable. Absolute
# bandpass phase and level are unobservable, so each connected piece of the
# array carries one arbitrary constant. Stations that hold across a boundary
# bridge it, and the ones that break are measured against them.

function shape_prior(bandwidth, σ)
    return OUPrior(; scale = LogNormal(log(bandwidth), 1.0), σ = LogNormal(log(σ), 1.5))
end

function _bandpass(ti, prior)
    return GainComponent(ConstantTerm(); Ti = ti, Frequency = ChannelBlocks(1), Feed = PerFeed(), prior = (; Frequency = prior))
end

function _both(ti, priors::BandpassPriors)
    return (;
        phase = (; bandpass = _bandpass(ti, shape_prior(priors.bandwidth, priors.phase_sigma))),
        logamp = (; bandpass = _bandpass(ti, shape_prior(priors.bandwidth, priors.logamp_sigma))),
    )
end

"""
    track_bandpass_model(breaks, origin, priors) -> GainModel

One bandpass per station and feed over the whole track, except for the stations
named in `breaks` (station => hours since `origin`), whose bandpass is cut at
those hours. Each piece is an Ornstein–Uhlenbeck process along frequency about a
free level under [`BandpassPriors`](@ref) `priors`.
"""
function track_bandpass_model(breaks, origin, priors::BandpassPriors)
    model = GainModel(; _both(GlobalTime(), priors)...)
    for (st, hrs) in pairs(breaks)
        isempty(hrs) && continue
        model = with_antenna(model, st; _both(InstrumentScans(track_epoch(sort(collect(hrs)), origin)), priors)...)
    end
    return model
end

"""
    scan_delay_model() -> GainModel

A delay per station and scan shared by its feeds, each other feed's offset
from the first, and a rate per station and scan, for a `FringeClosure` on
bandpass-corrected data. OU priors along time (10 ns and 0.1 Hz, one hour) let
a scan with no detection take its value from its neighbors instead of a flag;
at these widths a scan with detections is set by its own data.
"""
function scan_delay_model()
    delay = OUPrior(; scale = 3600.0, σ = 1.0e-8)
    rate = OUPrior(; scale = 3600.0, σ = 0.1)
    return GainModel(
        phase = (
            mbd = GainComponent(Delay(); Ti = PerScan(), Feed = SharedFeeds(), prior = delay),
            rel_delay = GainComponent(Delay(); Ti = PerScan(), Feed = ExceptFeed(1), prior = delay),
            rate = GainComponent(Rate(); Ti = PerScan(), Feed = SharedFeeds(), prior = rate),
        ),
    )
end

"""
    with_scan_terms(model, groups) -> GainModel

`model` with a phase level and a delay per scan and per channel group beside
the bandpass of each station `groups` names (station => channel ranges), both
common to the station's feeds. A phased array whose phase slope restarts at the
edges of its channel blocks and changes from scan to scan needs these; a delay
over the whole band cannot follow it.
"""
function with_scan_terms(model, groups)
    for (st, chans) in pairs(groups)
        fg = FreqGroups(chans)
        level = GainComponent(ConstantTerm(); Ti = PerScan(), Frequency = fg, Feed = SharedFeeds())
        delay = GainComponent(Delay(); Ti = PerScan(), Frequency = fg, Feed = SharedFeeds())
        tree = antenna_components(model, st)
        model = with_antenna(model, st; phase = merge(tree.phase, (; level, delay)), tree.logamp)
    end
    return model
end

"""
    station_bandpass_params(model, geom) -> Dict{String, Int}

The free parameters each station's bandpass carries under `model`: one complex
gain per channel, feed and occupied time segment, counted as the phase and
log-amplitude halves the solve stores.
"""
function station_bandpass_params(model, geom)
    nchan = length(geom.channel_freqs)
    return Dict(
        String(a) => nchan * 2 * geom.nfeed * last(Gustavo.Calibration.time_segment_ids(antenna_components(model, a).phase.bandpass.Ti, geom))
            for a in geom.antennas
    )
end

"The solution without its per-scan terms: what scans left out of the solve take."
shape_only(sol) = filter(c -> !(c.component.Ti isa PerScan), sol)
