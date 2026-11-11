"""Write an ANTAB file of ngehtsim's simulated a-priori SEFDs for one track.

For each station, scan and IF of the track's HOPS uvfits, the SEFD follows
ngehtsim's thermal model (``station_observation._station_terms_from_rows``):

    tau  = tau_zenith / sin(el)
    Tb   = T_CMB exp(-tau) + Tatm (1 - exp(-tau))
    Tsys = (T_R + eta_ff Tb + (1 - eta_ff) Tgnd) (1 + SSR)
    SEFD = 2 k Tsys / A_eff * wind factor

and the a-priori value written is SEFD exp(tau), corrected for the opacity
as ANTAB system temperatures are. The weather is the site's record for each
UT day the track covers (or the month's climatology, with --weather median and
the other ngehtsim forms), at each IF's center frequency.

The receivers are a named ngehtsim array's (--array) or ngehtsim's defaults for
a list of sites (--sites); ngehtsim drops a named array's settings when a site
list is given, so exactly one is accepted. --overrides names a TOML file of
ngehtsim's per-station overrides (T_R_overrides, sideband_ratio_overrides, ...),
passed to its obs_generator as given.

The file has a DPFU of 1 and a flat gain curve, so the TSYS values read as
SEFDs, one value per scan and IF, shared by both feeds.

    uv run --project python python python/simulate_antab.py \
        --out sim.AN (--array ARRAY | --sites ALMA,KP,...) [--overrides overrides.toml] \
        --stations AA=ALMA,KT=KP,... [--weather exact] file.uvfits ...
"""

import argparse
import datetime
import tomllib

import numpy as np

import ngehtsim.const_def as const
from ngehtsim.obs.obs_generator import obs_generator, windspeed_SEFD_modification
from ngehtsim.obs.observation_geometry import station_geometry_from_rows
from ngehtsim.obs.observation_template import ObservationTemplate

MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
MJD_UNIX_ZERO = 40587.0
OVERRIDES = (
    "T_R_overrides",
    "sideband_ratio_overrides",
    "D_overrides",
    "ap_eff_overrides",
    "surf_rms_overrides",
    "receiver_configuration_overrides",
    "bandwidth_overrides",
    "lo_freq_overrides",
    "hi_freq_overrides",
    "wind_loading_overrides",
    "custom_receivers",
)


def ut_date(mjd):
    return (datetime.datetime(1970, 1, 1) + datetime.timedelta(days=float(mjd) - MJD_UNIX_ZERO)).date()


def scans(time_mjd, gap_seconds=59.4):
    """Index arrays of the rows of each scan: runs of samples with no gap longer than `gap_seconds`."""
    order = np.argsort(time_mjd, kind="stable")
    t = time_mjd[order]
    cuts = np.flatnonzero(np.diff(t) * 86400.0 > gap_seconds) + 1
    return np.split(order, cuts)


def station_samples(path, stations):
    """Per station: (time_mjd, elevation_rad) of every row it is on, and the file's IF frequencies."""
    ds = ObservationTemplate.from_uvfits(path).dataset
    geom = station_geometry_from_rows(
        ds.stations.position_itrs_m, ds.time_mjd, ds.antenna1, ds.antenna2, ds.ra_hours, ds.dec_degrees
    )
    if geom is None:
        raise ValueError(f"{path}: a station has no position in the antenna table")
    out = {}
    for i, name in enumerate(ds.stations.names):
        if name not in stations:
            raise ValueError(f"{path}: station {name} has no ngehtsim site in --stations")
        first, second = ds.antenna1 == i, ds.antenna2 == i
        t = np.concatenate((ds.time_mjd[first], ds.time_mjd[second]))
        el = np.concatenate((geom.elevation1_rad[first], geom.elevation2_rad[second]))
        out[name] = (t, el)
    return out, np.asarray(ds.channel_frequency_hz), float(ds.channel_bandwidth_hz[0]), ds.ra_hours, ds.dec_degrees


def read_overrides(path):
    if path is None:
        return {}
    with open(path, "rb") as io:
        overrides = tomllib.load(io)
    bad = sorted(set(overrides) - set(OVERRIDES))
    if bad:
        raise ValueError(f"{path}: unknown ngehtsim override(s) {bad}; valid ones are {list(OVERRIDES)}")
    return overrides


def context(receivers, day, freq_hz, bandwidth_hz, ra_hours, dec_degrees, weather, sites):
    """ngehtsim's station context; `receivers` holds `array`, `sites` and the overrides."""
    gen = obs_generator(
        settings=dict(
            array=receivers["array"],
            sites=receivers["sites"],
            RA=ra_hours,
            DEC=dec_degrees,
            frequency=freq_hz / 1.0e9,
            bandwidth=bandwidth_hz / 1.0e9,
            month=MONTHS[day.month - 1],
            day=str(day.day),
            year=str(day.year),
            weather=weather,
            weather_year=day.year,
            weather_day=day.day,
            random_seed=0,
        ),
        **receivers["overrides"],
    )
    ctx = gen.station_context()
    missing = sorted(set(sites) - set(ctx["tau"]))
    if missing:
        source = f"array {receivers['array']}" if receivers["array"] else "--sites"
        raise ValueError(f"site(s) {missing} in --stations are not in {source}; ngehtsim has {sorted(map(str, ctx['tau']))}")
    return ctx


def apriori_sefd(ctx, site, el):
    """ngehtsim's SEFD times exp(tau) along the line of sight, at elevations `el` (radians)."""
    tau = ctx["tau"][site] / np.sin(el)
    Tb = const.T_CMB * np.exp(-tau) + ctx["Tatm"][site] * (1.0 - np.exp(-tau))
    Tsys = (ctx["receiver_temperature"][site] + const.eta_ff * Tb + (1.0 - const.eta_ff) * ctx["Tgnd"][site]) * (
        1.0 + ctx["sideband_ratio"][site]
    )
    sefd = 2.0 * const.k * Tsys / ctx["effective_area"][site]
    wl = ctx["wind_loading"][site]
    sefd = sefd * windspeed_SEFD_modification(ctx["windspeed"][site], wl["v0"], wl["w"])
    return sefd * np.exp(tau)


def antab_time(mjd):
    d = datetime.datetime(1970, 1, 1) + datetime.timedelta(days=float(mjd) - MJD_UNIX_ZERO)
    minutes = d.hour * 60 + d.minute + d.second / 60.0
    return f"{d.timetuple().tm_yday:3d} {int(minutes // 60):02d}:{minutes % 60:04.1f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--array", help="ngehtsim array whose sites and receiver settings apply")
    group.add_argument("--sites", help="ngehtsim site names, comma-separated, with ngehtsim's default receivers")
    parser.add_argument("--overrides", help="TOML file of ngehtsim per-station overrides")
    parser.add_argument("--stations", required=True, help="CODE=SITE pairs, comma-separated")
    parser.add_argument("--weather", default="exact")
    parser.add_argument("uvfits", nargs="+")
    args = parser.parse_args()
    stations = dict(pair.split("=") for pair in args.stations.split(","))
    receivers = dict(
        array=args.array,
        sites=args.sites.split(",") if args.sites else None,
        overrides=read_overrides(args.overrides),
    )

    freqs = None
    contexts = {}
    rows = {}
    for path in args.uvfits:
        per_station, f, bw, ra, dec = station_samples(path, stations)
        if freqs is None:
            freqs = f
        elif not np.array_equal(freqs, f):
            raise ValueError(f"{path}: its IF frequencies differ from the track's other files")
        for name, (t, el) in per_station.items():
            for scan in scans(t):
                tm = float(np.mean(t[scan]))
                day = ut_date(tm)
                values = []
                for fc in freqs:
                    key = (day, float(fc), ra, dec)
                    if key not in contexts:
                        contexts[key] = context(receivers, day, fc, bw, ra, dec, args.weather, stations.values())
                    values.append(float(np.mean(apriori_sefd(contexts[key], stations[name], el[scan]))))
                rows.setdefault(name, []).append((tm, values))

    nif = len(freqs)
    with open(args.out, "w") as io:
        io.write("! Simulated a-priori SEFDs (Jy) from ngehtsim's weather, receiver and dish tables.\n")
        io.write("! TSYS holds SEFD; DPFU = 1 and a flat gain curve make each value read as one.\n")
        for name in sorted(rows):
            io.write(f"GAIN {name} ELEV DPFU = 1.0, 1.0 POLY = 1.0 /\n")
        for name in sorted(rows):
            index = ", ".join(f"'L{i}|R{i}'" for i in range(1, nif + 1))
            io.write(f"TSYS {name} FT = 1.0 TIMEOFF = 0 INDEX = {index} /\n")
            for tm, values in sorted(rows[name]):
                io.write(f"{antab_time(tm)} " + " ".join(f"{x:.4f}" for x in values) + "\n")
            io.write("/\n")


if __name__ == "__main__":
    main()
