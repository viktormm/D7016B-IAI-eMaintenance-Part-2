"""Simulated monitoring of the track from ONE scan: degradation, detection, alarm, maintenance.

The track is "scanned" at regular intervals. Between runs it degrades by an assumed model; every run is measured by
the unchanged pipeline (track_curve.py); indicators are checked against limits; an alarm orders maintenance at that
place, which restores the geometry there; then the degradation starts again. Nothing is scripted to happen in a
given week: maintenance only happens when the measurements raise an alarm.

  the track       three defects that grow linearly within a maintenance cycle (as in track geometry degradation models,
                  see Soleimanmeigouni et al., 2018) and are reset by maintenance:
                    gauge widening, 100-120 m        the left rail moves outward
                    alignment defect, about 210 m    both rails move sideways
                    short-wave gauge irregularity, 227-248 m   the gauge varies with a 1 m wavelength
  a run           a different pseudo-random half of the point cloud (its own sampling noise, like a real run), with
                  the defects of that week added, built by the pipeline; optionally with a sideways positioning
                  error (--position-sigma), which is taken out again (median offset) before runs are compared
  indicators      against the first run: the gauge change and the sideways change of the centreline as a 5 m mean
                  along the curve, and the increase of the TQI (gauge) of every 25 m section
  limits          detection limits are calibrated: runs WITHOUT degradation are compared with the first run, and the
                  largest value they reach, times a safety factor, is the limit. Above it: ALERT.
                  intervention limits are assumed (INTERVENTION, TQI_ALARM_FACTOR). At or above them: ALARM.
  maintenance     an ALARM orders the matching work on the alarmed stretch (gauge correction, lining, fastening
                  renewal); it is done before the next run and resets the defects of that kind there.

The sequence shows the behaviour of the method, not of the track: the degradation and the intervention limits are
assumed, and the rates are chosen so that a whole cycle fits in the simulated period.

  python simulate_runs.py gps_pointcloud_00.laz --reference output/reference/track.json [--runs 13] [--interval-weeks 2]
  python simulate_runs.py gps_pointcloud_00.laz --reference output/reference/track.json --view-only

Writes the runs to output/runs/sim_wNN/ (NN = week; their meta.json says they are simulated and lists the findings),
a summary to output/simulation/ and the 3D view output/viewer/curve3d_simulation.html. All seeds are fixed, so the
output is the same every time.
"""
import argparse
import os
import webbrowser

import numpy as np

import track_curve as tc

OUT = os.path.join(tc.OUT, "simulation")
DEFECTS = (                       # the assumed degradation: where, how it moves the rails, growth per week (mm)
    {"name": "gauge widening", "kind": "widen", "lo": 100.0, "hi": 120.0, "rate": 0.25, "indicator": "gauge"},
    {"name": "alignment defect", "kind": "bump", "lo": 194.0, "hi": 226.0, "rate": 0.35, "indicator": "centre"},
    {"name": "short-wave gauge irregularity", "kind": "rough", "lo": 227.0, "hi": 248.0, "rate": 0.10, "indicator": "tqi"},
)
ROUGH_WAVELENGTH = 1.0            # wavelength of the short-wave gauge irregularity (m)
INTERVENTION = {"gauge": 3.0, "centre": 3.5}   # assumed intervention limits for the 5 m mean change (mm)
TQI_ALARM_FACTOR = 2.0            # assumed: TQI alarm at this multiple of its detection limit
WORK = {"gauge": "gauge correction", "centre": "lining", "tqi": "fastening renewal"}
LABEL = {"gauge": "gauge change", "centre": "alignment change", "tqi": "TQI (gauge) increase"}
MEAN_OVER = 5.0                   # the change along the curve is judged as a mean over this length (m)
CALIBRATION_RUNS = 4              # runs without degradation that set the detection limits
SAFETY = 1.2                      # detection limit = this times the largest value the calibration runs reach
EDGE = 5.0                        # the first and last metres of the reference are left out (m)
STEP = 0.25                       # grid of the comparison (m), as in track_curve.compare_curves


# ------------------------------------------------------------------ the assumed degradation


def _window(s, lo, hi, ramp=2.0):
    x = np.clip((s - lo) / ramp, 0, 1) * np.clip((hi - s) / ramp, 0, 1)
    return x * x * (3 - 2 * x)


def profile(defect, s):
    """Shape of a defect along the track; 1 is its full size."""
    lo, hi = defect["lo"], defect["hi"]
    if defect["kind"] == "widen":
        return _window(s, lo, hi)
    if defect["kind"] == "bump":
        return np.exp(-0.5 * ((s - 0.5 * (lo + hi)) / ((hi - lo) / 4)) ** 2)
    return np.sin(2 * np.pi * s / ROUGH_WAVELENGTH) * _window(s, lo, hi)


def shift_mm(s, d, sizes):
    """Sideways shift of every point: the alignment defect moves both rails, the other two the left rail only."""
    out = np.zeros_like(s)
    for defect, size in zip(DEFECTS, sizes):
        if size:
            out += size * profile(defect, s) * (1.0 if defect["kind"] == "bump" else (d > 0))
    return out


def expected_tqi_increase(size):
    """TQI of the short-wave irregularity alone over its section, which is about what it adds to the section's TQI."""
    rough = next(df for df in DEFECTS if df["kind"] == "rough")
    sec = int(rough["lo"] // tc.P.section)
    sc = np.arange(sec * tc.P.section + tc.P.slab / 2, (sec + 1) * tc.P.section, tc.P.slab)
    return float(tc._tqi(size / 1e3 * profile(rough, sc), tc.P.slab))


# ------------------------------------------------------------------ indicators, limits, findings


def moving_mean(v, n):
    """Mean over n samples; not given where less than 80 % of the window has a value (ends of the track, gaps)."""
    ok = np.isfinite(v)
    total = np.convolve(np.where(ok, v, 0.0), np.ones(n), mode="same")
    count = np.convolve(ok.astype(float), np.ones(n), mode="same")
    return np.where(count >= 0.8 * n, total / np.maximum(count, 1), np.nan)


def indicators(base, run, ref, registered):
    """What is judged in a run, all against the first run: change along the curve (mm) and TQI increase per section."""
    grid = np.arange(0, ref.length, STEP)
    gauge, centre = (np.full(len(grid), np.nan) for _ in range(2))
    diff = tc.compare_curves(base["curve"], run["curve"], step=STEP)
    i = np.round(diff["s"] / STEP).astype(int)
    gauge[i], centre[i] = diff["gauge"] * 1e3, diff["d_centre"] * 1e3
    if registered:                                        # a run's own sideways offset is taken out first
        centre -= np.nanmedian(centre)
    edge = (grid < EDGE) | (grid > ref.length - EDGE)
    gauge[edge], centre[edge] = np.nan, np.nan
    n = int(round(MEAN_OVER / STEP))
    rows = {int(r["section"]): r for r in tc.section_table(run)}
    first = {int(r["section"]): r for r in tc.section_table(base)}
    return {"grid": grid, "gauge": gauge, "centre": centre, "gauge_mean": moving_mean(gauge, n), "centre_mean": moving_mean(centre, n),
            "sections": rows, "dtqi": {k: rows[k]["tqi_gauge"] - first[k]["tqi_gauge"] for k in rows if k in first}}


def calibrate(base, quiet_runs, ref):
    """Detection limits: the largest value runs without degradation reach against the first run, times SAFETY."""
    ind = [indicators(base, r, ref, False) for r in quiet_runs]
    top = lambda key: SAFETY * max(float(np.nanmax(np.abs(i[key]))) for i in ind)
    return {"gauge": top("gauge_mean"), "centre": top("centre_mean"),
            "tqi": SAFETY * float(np.nanmax([abs(v) for i in ind for v in i["dtqi"].values()]))}


def findings(ind, limits):
    """Where an indicator is above its detection limit (ALERT) or at its intervention limit (ALARM)."""
    out, n = [], int(round(MEAN_OVER / STEP))
    for q in ("gauge", "centre"):
        m = ind[q + "_mean"]
        over = np.isfinite(m) & (np.abs(m) > limits[q])
        ends = np.flatnonzero(np.diff(np.concatenate([[0], over.astype(int), [0]]))).reshape(-1, 2).tolist()
        merged = []
        for i0, i1 in ends:                               # stretches closer together than MEAN_OVER are one finding
            if merged and i0 - merged[-1][1] < n:
                merged[-1][1] = i1
            else:
                merged.append([i0, i1])
        for i0, i1 in merged:
            peak = float(m[i0:i1][np.nanargmax(np.abs(m[i0:i1]))])
            alarm = abs(peak) >= INTERVENTION[q]
            out.append({"level": "ALARM" if alarm else "ALERT", "indicator": q, "from": float(ind["grid"][i0]), "to": float(ind["grid"][i1 - 1]),
                        "value": peak, "limit": INTERVENTION[q] if alarm else limits[q], "unit": " mm"})
    for sec, dt in sorted(ind["dtqi"].items()):
        if dt > limits["tqi"]:
            alarm = dt >= TQI_ALARM_FACTOR * limits["tqi"]
            out.append({"level": "ALARM" if alarm else "ALERT", "indicator": "tqi", "from": sec * tc.P.section, "to": (sec + 1) * tc.P.section,
                        "value": float(dt), "limit": (TQI_ALARM_FACTOR if alarm else 1) * limits["tqi"], "unit": ""})
    for f in out:
        f["text"] = f"{f['level']}: {LABEL[f['indicator']]} {f['value']:+.1f}{f['unit']} at {f['from']:.0f}–{f['to']:.0f} m"
    return out


# ------------------------------------------------------------------ the loop: degrade, scan, judge, maintain


def simulate(laz, ref, n_runs, interval, position_sigma):
    pts = tc.read_points(laz)
    s, d = ref.project(pts.E, pts.N)
    near = (np.abs(d) < tc.P.d_max + 0.2) & (s > -1) & (s < ref.length + 1)      # only these can end up in a run
    E, N, Z, s, d = pts.E[near], pts.N[near], pts.Z[near], s[near], d[near]
    nx, ny = ref.normal_at(s)
    blocks = np.arange(len(E)) // 1000

    def scan(run_id, seed, sizes=None, offset=0.0):
        shift = ((shift_mm(s, d, sizes) if sizes else 0.0) + offset) / 1e3
        half = (np.random.RandomState(seed).rand(blocks[-1] + 1) < 0.5)[blocks]
        return tc.build_run(tc.Points(E + shift * nx, N + shift * ny, Z, pts.name, pts.sha256), ref, tc.P, run_id, half)

    base = scan("sim_w00", 1000)
    limits = calibrate(base, [scan(f"quiet{i}", 3000 + i) for i in range(CALIBRATION_RUNS)], ref)
    print(f"detection limits from {CALIBRATION_RUNS} runs without degradation: gauge change {limits['gauge']:.2f} mm, "
          f"alignment change {limits['centre']:.2f} mm (5 m means), TQI increase {limits['tqi']:.1f}")
    print(f"intervention limits (assumed): gauge change {INTERVENTION['gauge']:.1f} mm, alignment change {INTERVENTION['centre']:.1f} mm, "
          f"TQI increase {TQI_ALARM_FACTOR * limits['tqi']:.1f}\n")

    since = [0] * len(DEFECTS)                            # week of the last maintenance of every defect
    history = []
    for k in range(n_runs):
        week = k * interval
        sizes = [df["rate"] * (week - since[j]) for j, df in enumerate(DEFECTS)]
        offset = float(np.random.RandomState(2000 + k).normal(0, position_sigma)) if position_sigma and k else 0.0
        run = base if k == 0 else scan(f"sim_w{week:02d}", 1000 + k, sizes, offset)
        ind = indicators(base, run, ref, position_sigma > 0)
        found = findings(ind, limits) if k else []
        for f in found:
            f["action"] = ""
            if f["level"] == "ALARM":                     # maintenance is ordered and done before the next run
                f["action"] = f"maintenance ordered: {WORK[f['indicator']]} at {f['from']:.0f}–{f['to']:.0f} m"
                for j, df in enumerate(DEFECTS):
                    if df["indicator"] == f["indicator"] and df["lo"] < f["to"] + MEAN_OVER / 2 and df["hi"] > f["from"] - MEAN_OVER / 2:
                        since[j] = week
        run["meta"]["simulation"] = {"note": "simulated run: a random half of one scan with an assumed degradation added", "week": week,
                                     "defects_mm": {df["name"]: round(size, 3) for df, size in zip(DEFECTS, sizes)},
                                     "position_error_mm": round(offset, 3), "seed": 1000 + k}
        run["meta"]["findings"] = [{"level": f["level"], "from": f["from"], "to": f["to"],
                                    "text": f["text"] + (" → " + f["action"] if f["action"] else "")} for f in found]
        tc.save_run(run, os.path.join(tc.OUT, "runs", run["meta"]["run_id"]))
        history.append({"week": week, "run": run["meta"]["run_id"], "sizes": sizes, "ind": ind, "found": found, "offset": offset})
        print(f"week {week:2d}  " + "  ".join(f"{df['name']} {size:4.2f} mm" for df, size in zip(DEFECTS, sizes)))
        for f in found:
            print(f"         {f['text']}" + (f"  ->  {f['action']}" if f["action"] else ""))
    return history, limits


# ------------------------------------------------------------------ summary


def at_defect(h, defect):
    """The indicator that should see this defect, measured at the defect, and the size that was put in."""
    ind, lo, hi = h["ind"], defect["lo"], defect["hi"]
    size = h["sizes"][DEFECTS.index(defect)]
    if defect["kind"] == "rough":
        return ind["dtqi"].get(int(lo // tc.P.section), np.nan), expected_tqi_increase(size)
    m = ind[defect["indicator"] + "_mean"]
    inside = (ind["grid"] >= lo + 3) & (ind["grid"] <= hi - 3) if defect["kind"] == "widen" else np.abs(ind["grid"] - 0.5 * (lo + hi)) < 2
    v = m[inside]
    return (float(v[np.nanargmax(np.abs(v))]) if np.isfinite(v).any() else np.nan), size


def write_csv(path, rows):
    cols = list(rows[0])
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) if isinstance(r[c], str) else f"{r[c]:.3f}" for c in cols) + "\n")


def summarise(history, limits):
    timeline, events, sections = [], [], []
    for h in history:
        row = {"run": h["run"], "week": h["week"], "position_error_mm": h["offset"]}
        for df in DEFECTS:
            measured, put_in = at_defect(h, df)
            key = df["name"].replace(" ", "_").replace("-", "_")
            row[key + "_size_mm"], row[key + "_indicator"], row[key + "_expected"] = h["sizes"][DEFECTS.index(df)], measured, put_in
        timeline.append(row)
        events += [{"week": h["week"], "run": h["run"], "level": f["level"], "indicator": LABEL[f["indicator"]], "from_m": f["from"],
                    "to_m": f["to"], "value": f["value"], "limit": f["limit"], "action": f["action"]} for f in h["found"]]
        sections += [dict(run=h["run"], week=h["week"], **r) for _, r in sorted(h["ind"]["sections"].items())]
    write_csv(os.path.join(OUT, "timeline.csv"), timeline)
    write_csv(os.path.join(OUT, "sections_over_time.csv"), sections)
    if events:
        write_csv(os.path.join(OUT, "events.csv"), events)
    with open(os.path.join(OUT, "limits.csv"), "w", newline="\n") as f:
        f.write("indicator,detection_limit,intervention_limit\n")
        f.write(f"gauge change (5 m mean) mm,{limits['gauge']:.3f},{INTERVENTION['gauge']:.3f}\n")
        f.write(f"alignment change (5 m mean) mm,{limits['centre']:.3f},{INTERVENTION['centre']:.3f}\n")
        f.write(f"TQI (gauge) increase,{limits['tqi']:.3f},{TQI_ALARM_FACTOR * limits['tqi']:.3f}\n")
    return timeline, events


def figures(history, limits, interval):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888", "axes.labelcolor": "#333", "xtick.color": "#555",
                         "ytick.color": "#555", "axes.spines.top": False, "axes.spines.right": False})
    weeks = np.array([h["week"] for h in history])
    grid = history[0]["ind"]["grid"]
    G, C = (np.array([h["ind"][q] for h in history]) for q in ("gauge", "centre"))
    ok = np.isfinite(G).any(0)
    extent = [grid[ok][0], grid[ok][-1], weeks[0] - interval / 2, weeks[-1] + interval / 2]
    cmap = LinearSegmentedColormap.from_list("change", ["#1f5fae", "#f1f3f3", "#d9822b"])     # narrower/right, none, wider/left
    cmap.set_bad("#ffffff")
    colour = {"gauge": "#D9822B", "centre": "#0E5A62", "tqi": "#6B4C9A"}
    alarm_red, ink = "#C62828", "#333"
    works = [(h["week"], f) for h in history for f in h["found"] if f["level"] == "ALARM"]

    fig, ax = plt.subplots(4, 1, figsize=(9, 11), gridspec_kw={"height_ratios": [1, 1, 1.1, 1.1]})
    for a, M, q, title in ((ax[0], G, "gauge", "Gauge change against the first run (mm)"),
                           (ax[1], C, "centre", "Sideways change of the centreline against the first run (mm)")):
        im = a.imshow(M[:, ok], aspect="auto", origin="lower", extent=extent, cmap=cmap, vmin=-5, vmax=5, interpolation="nearest")
        a.set_title(title, fontsize=10, loc="left"); a.set_ylabel("week"); a.set_yticks(weeks[::2])
        a.set_ylim(extent[2], extent[3] + interval)           # room for work ordered after the last run
        fig.colorbar(im, ax=a, pad=0.01, fraction=0.03)
        for week, f in works:                                 # where and when an alarm ordered maintenance
            if f["indicator"] == q or (q == "gauge" and f["indicator"] == "tqi"):
                a.plot([f["from"], f["to"]], [week + interval / 2] * 2, color=ink, lw=1.6, solid_capstyle="butt")
                near_end = f["to"] > extent[1] - 40           # no room on the right: the label goes on the left
                a.text(f["from"] - 2 if near_end else f["to"] + 2, week + interval / 2, WORK[f["indicator"]], fontsize=7.5, color=ink,
                       va="center", ha="right" if near_end else "left")
    ax[1].set_xlabel("chainage s (m)")

    def series(a, df, label):
        vals = [at_defect(h, df) for h in history]
        c, q = colour[df["indicator"]], df["indicator"]
        a.plot(weeks, [v[1] for v in vals], color=c, lw=1, ls=(0, (4, 3)))
        a.plot(weeks, [v[0] for v in vals], color=c, lw=1.6, marker="o", ms=4, label=label)
        hit = [(h["week"], v[0]) for h, v in zip(history, vals) if any(f["level"] == "ALARM" and f["indicator"] == q and
               f["from"] < df["hi"] + MEAN_OVER and f["to"] > df["lo"] - MEAN_OVER for f in h["found"])]
        if hit:
            a.plot(*zip(*hit), ls="none", marker="o", ms=10, mfc="none", mec=alarm_red, mew=1.6)
        return c

    a = ax[2]
    for df in DEFECTS[:2]:
        c = series(a, df, f"{df['name']}, {df['lo']:.0f}–{df['hi']:.0f} m")
        q = df["indicator"]
        a.axhline(limits[q], color=c, lw=0.8, ls=":"); a.axhline(INTERVENTION[q], color=c, lw=0.8)
        a.text(weeks[-1] + 0.4, limits[q], f"detection {limits[q]:.2f}", fontsize=7.5, color=c, va="bottom" if q == "gauge" else "top")
        a.text(weeks[-1] + 0.4, INTERVENTION[q], f"intervention {INTERVENTION[q]:.1f}", fontsize=7.5, color=c, va="bottom")
    a.set_title("Change along the curve at the two defects (5 m mean): measured (line), put in (dashed), alarm (red ring)", fontsize=10, loc="left")
    a.set_ylabel("mm"); a.legend(frameon=False, loc="upper left", fontsize=8)

    a = ax[3]
    rough = DEFECTS[2]
    sec = int(rough["lo"] // tc.P.section)
    for other in sorted(history[-1]["ind"]["dtqi"]):
        if other != sec:
            a.plot(weeks, [h["ind"]["dtqi"].get(other, np.nan) for h in history], color="#C8C8C8", lw=0.9)
    c = series(a, rough, f"{rough['name']}, section {sec * tc.P.section:.0f}–{(sec + 1) * tc.P.section:.0f} m")
    a.axhline(limits["tqi"], color=c, lw=0.8, ls=":"); a.axhline(TQI_ALARM_FACTOR * limits["tqi"], color=c, lw=0.8)
    a.text(weeks[-1] + 0.4, limits["tqi"], f"detection {limits['tqi']:.1f}", fontsize=7.5, color=c, va="center")
    a.text(weeks[-1] + 0.4, TQI_ALARM_FACTOR * limits["tqi"], f"intervention {TQI_ALARM_FACTOR * limits['tqi']:.1f}", fontsize=7.5, color=c, va="center")
    a.set_title("TQI (gauge) increase per section: measured (line), expected (dashed), alarm (red ring); grey: other sections", fontsize=10, loc="left")
    a.set_ylabel("TQI increase"); a.set_xlabel("week"); a.legend(frameon=False, loc="upper left", fontsize=8)
    for a in ax[2:]:
        a.set_xticks(weeks); a.set_xlim(weeks[0] - 0.5, weeks[-1] + 5); a.grid(axis="y", color="#eee", lw=0.6); a.set_axisbelow(True)
        for week in sorted({w for w, _ in works}):
            a.axvline(week + interval / 2, color=ink, lw=0.7, ls=(0, (4, 3)))
    plt.tight_layout(); plt.savefig(os.path.join(OUT, "sim_spacetime.png"), dpi=140); plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    hot = {int(df["lo"] // tc.P.section): colour[df["indicator"]] for df in (DEFECTS[0], DEFECTS[2])}
    for sec in sorted(history[0]["ind"]["sections"]):
        for a, key in zip(ax, ("gauge_dev_mean_mm", "tqi_gauge")):
            y = [h["ind"]["sections"][sec][key] if sec in h["ind"]["sections"] else np.nan for h in history]
            a.plot(weeks, y, color=hot.get(sec, "#C8C8C8"), lw=1.6 if sec in hot else 0.9, marker="o" if sec in hot else None, ms=4, zorder=3 if sec in hot else 2)
            if sec in hot:
                a.text(weeks[-1] + 0.5, y[-1], f"{sec * tc.P.section:.0f}–{(sec + 1) * tc.P.section:.0f} m", va="center", fontsize=8, color=hot[sec])
    for a, title, unit in zip(ax, ("Section mean gauge deviation over time", "Section TQI (gauge) over time"), ("mm", "TQI")):
        a.set_title(title, fontsize=10, loc="left"); a.set_xlabel("week"); a.set_ylabel(unit); a.set_xticks(weeks[::2])
        a.grid(axis="y", color="#eee", lw=0.6); a.set_axisbelow(True); a.set_xlim(weeks[0] - 0.5, weeks[-1] + 5)
    ax[0].text(0.02, 0.03, "grey: the other sections", transform=ax[0].transAxes, fontsize=8, color="#777")
    plt.tight_layout(); plt.savefig(os.path.join(OUT, "sim_sections.png"), dpi=140); plt.close()


def write_view(run_ids, laz, open_it):
    import view_curve as vc
    data, _ = vc.build_data(run_ids, 20000, laz)
    html = os.path.join(tc.OUT, "viewer", "curve3d_simulation.html")
    os.makedirs(os.path.dirname(html), exist_ok=True)
    with open(html, "w", encoding="utf-8", newline="\n") as f:
        f.write(vc.render(data))
    print(f"3D view in {html}")
    if open_it:
        webbrowser.open("file:///" + os.path.abspath(html).replace(os.sep, "/"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("laz"); ap.add_argument("--reference", required=True)
    ap.add_argument("--view-only", action="store_true", help="only rebuild the 3D view from the simulated runs already in output/runs/")
    ap.add_argument("--runs", type=int, default=13, help="number of runs, the first one is the baseline (default 13)")
    ap.add_argument("--interval-weeks", type=int, default=2, help="time between runs (default 2 weeks)")
    ap.add_argument("--position-sigma", type=float, default=0.0, help="sideways positioning error per run, standard deviation in mm (default 0)")
    ap.add_argument("--no-open", action="store_true", help="do not open the 3D view in the browser")
    a = ap.parse_args()
    if a.view_only:
        folder = os.path.join(tc.OUT, "runs")
        write_view(sorted(n for n in os.listdir(folder) if n.startswith("sim_w")), a.laz, not a.no_open)
        return
    os.makedirs(OUT, exist_ok=True)
    ref = tc.load_reference(a.reference)
    history, limits = simulate(a.laz, ref, a.runs, a.interval_weeks, a.position_sigma)
    timeline, events = summarise(history, limits)
    figures(history, limits, a.interval_weeks)
    alarms = [e for e in events if e["level"] == "ALARM"]
    print(f"\n{len(events)} findings in {len(history) - 1} runs, {len(alarms)} of them alarms that ordered maintenance:")
    for e in alarms:
        print(f"  week {e['week']:2d}: {e['indicator']} {e['value']:+.1f} at {e['from_m']:.0f}–{e['to_m']:.0f} m (limit {e['limit']:.1f})  ->  {e['action']}")
    print(f"summary in {OUT}")
    write_view([h["run"] for h in history], a.laz, not a.no_open)


if __name__ == "__main__":
    main()
