"""Interactive 3D view of the run curves, plus the same curves as a LAS file for CloudCompare.

  python view_curve.py run1 [run2 ...] [--points 60000] [--laz FILE] [--out DIR] [--no-open]

Writes to output/viewer/:
  curve3d_<runs>.html   one self-contained page (works offline): rotate, zoom, colour by gauge, cursor along the track
  <run>_curve.las       the Catmull-Rom curve of each rail and the centreline every 0.1 m, in the point cloud's own
                        coordinates, so it can be loaded into CloudCompare on top of the original .laz

The curves are the ones `track_curve.py sample` writes: recomputed from the run's control points (curve.csv), so the
view never shows a stale file. Nothing here changes a run or the pipeline, and the output is byte-identical when the
inputs are the same (the point cloud sample is thinned evenly, not at random).
"""
import argparse
import datetime
import json
import os
import sys
import webbrowser

import laspy
import numpy as np
from plotly.offline import get_plotlyjs

import catmull_rom as cr
import track_curve as tc

STEP = 0.1                        # spacing of the sampled curve (m), the same as `track_curve.py sample`
CROP_D = 1.5                      # point cloud context: within this sideways offset of the reference (m)
Z_BELOW, Z_ABOVE = 0.36, 0.15     # ... and within this height range around the rail heads (m)
GAUGE_HALF_RANGE = 4.0            # gauge colour scale: the median gauge of the runs +- this (mm)
PALETTE = ["#028090", "#D9822B", "#0B2027", "#6FBAC3"]
LAS_COLOURS = {"left": (217, 130, 43), "right": (2, 128, 144), "centre": (150, 160, 165)}
LAS_DATE = datetime.date(2026, 1, 1)   # fixed, so the LAS file is byte-identical every time


def _list(a, decimals):
    return [None if not np.isfinite(v) else round(float(v), decimals) for v in np.asarray(a, float)]


def find_reference(ref_id):
    folder = os.path.join(tc.OUT, "reference")
    for name in sorted(os.listdir(folder)) if os.path.isdir(folder) else []:
        path = os.path.join(folder, name)
        if name.endswith(".json") and json.load(open(path)).get("id") == ref_id:
            return tc.load_reference(path)
    sys.exit(f"no reference frame with id {ref_id} in {folder}")


def load_curves(run_id):
    """The run's control points and its Catmull-Rom curve sampled every STEP m (NaN where the curve has a gap)."""
    run = tc.load_run(os.path.join(tc.OUT, "runs", run_id))
    c = run["curve"]
    dense = tc.sample_curve(c, STEP)
    for side in ("left", "right"):
        dense["d_" + side] = cr.evaluate(c["s"], c["d_" + side], dense["s"], 3 * tc.P.ctrl)
    ok = np.flatnonzero(np.isfinite(dense["E_centre"]) & np.isfinite(dense["Z_centre"]))
    lo, hi = ok[0], ok[-1] + 1
    dense = {k: v[lo:hi] for k, v in dense.items()}          # trim the ends, keep gaps inside as gaps
    return run, c, dense


def json_run(run, c, dense, colour, origin):
    E0, N0 = origin
    g = dense
    gauge = (g["d_left"] - g["d_right"]) * 1e3
    cg = (c["d_left"] - c["d_right"]) * 1e3
    return {"id": run["meta"]["run_id"], "colour": colour,
            "findings": run["meta"].get("findings", []),      # alerts and alarms of a simulated run (simulate_runs.py), if any
            "meta": {"reference": run["meta"]["reference"], "code": run["meta"]["code"],
                     "input": run["meta"]["input"]["file"], "sha256": run["meta"]["input"]["sha256"][:12]},
            "dense": {"s": _list(g["s"], 3), "dL": _list(g["d_left"], 5), "dR": _list(g["d_right"], 5),
                      "zL": _list(g["Z_left"], 5), "zR": _list(g["Z_right"], 5),
                      "EL": _list(g["E_left"] - E0, 3), "NL": _list(g["N_left"] - N0, 3),
                      "ER": _list(g["E_right"] - E0, 3), "NR": _list(g["N_right"] - N0, 3), "gauge": _list(gauge, 2)},
            "ctrl": {"s": _list(c["s"], 3), "dL": _list(c["d_left"], 5), "dR": _list(c["d_right"], 5),
                     "zL": _list(c["z_left"], 5), "zR": _list(c["z_right"], 5),
                     "EL": _list(c["E_left"] - E0, 3), "NL": _list(c["N_left"] - N0, 3),
                     "ER": _list(c["E_right"] - E0, 3), "NR": _list(c["N_right"] - N0, 3), "gauge": _list(cg, 2)}}


def load_points(sha256, ref, laz, n_points, origin, s_range, z_rail):
    """An even, deterministic sample of the raw points around the rails of the first run, as context."""
    if not laz or not os.path.exists(laz):
        print(f"point cloud {laz} not found, showing the curves only")
        return None
    pts = tc.read_points(laz)
    if pts.sha256 != sha256:
        print("the point cloud is not the file this run was built from, showing the curves only")
        return None
    stride = max(1, len(pts.E) // (n_points * 8))
    E, N, Z = pts.E[::stride], pts.N[::stride], pts.Z[::stride]
    s, d = ref.project(E, N)
    keep = (np.abs(d) <= CROP_D) & (s >= s_range[0]) & (s <= s_range[1]) & (Z >= z_rail - Z_BELOW) & (Z <= z_rail + Z_ABOVE)
    idx = np.flatnonzero(keep)
    if len(idx) > n_points:
        idx = idx[np.linspace(0, len(idx) - 1, n_points).round().astype(int)]
    z = Z[idx]
    return {"s": _list(s[idx], 3), "d": _list(d[idx], 4), "z": _list(z, 4),
            "E": _list(E[idx] - origin[0], 3), "N": _list(N[idx] - origin[1], 3),
            "zrange": [round(float(np.percentile(z, 1)), 4), round(float(np.percentile(z, 99)), 4)],
            "file": os.path.basename(laz), "stride": int(stride)}


def export_las(dense, path):
    """Left rail, right rail and centreline as points every STEP m, with chainage and gauge as scalar fields.

    Each line keeps its own points: where one rail is missing, the other rail (and only that one) is still written.
    """
    gauge = (dense["d_left"] - dense["d_right"]) * 1e3
    parts = []
    for rail, side in enumerate(("left", "right", "centre")):
        E, N, Z = (dense[f"{a}_{side}"] for a in "ENZ")
        m = np.isfinite(E) & np.isfinite(N) & np.isfinite(Z)
        parts.append({"x": E[m], "y": N[m], "z": Z[m], "s": dense["s"][m], "gauge": gauge[m], "rail": np.full(m.sum(), rail, np.uint8),
                      "rgb": np.tile(np.array(LAS_COLOURS[side]) * 257, (m.sum(), 1)).astype(np.uint16)})
    cat = lambda k: np.concatenate([p[k] for p in parts])
    x, y, z = cat("x"), cat("y"), cat("z")
    hdr = laspy.LasHeader(point_format=2, version="1.2")
    hdr.creation_date = LAS_DATE
    hdr.scales = np.array([1e-4, 1e-4, 1e-4])
    hdr.offsets = np.floor([x.min(), y.min(), z.min()])
    hdr.add_extra_dim(laspy.ExtraBytesParams(name="s_m", type="f4", description="chainage along the reference"))
    hdr.add_extra_dim(laspy.ExtraBytesParams(name="gauge_mm", type="f4", description="gauge, mm"))
    hdr.add_extra_dim(laspy.ExtraBytesParams(name="rail", type="u1", description="0 left, 1 right, 2 centre"))
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = x, y, z
    las.s_m = cat("s").astype(np.float32)
    las.gauge_mm = cat("gauge").astype(np.float32)
    las.rail = cat("rail")
    las.red, las.green, las.blue = (cat("rgb")[:, k] for k in range(3))
    las.write(path)
    return len(x)


def build_data(run_ids, n_points, laz):
    loaded = [load_curves(r) for r in run_ids]
    for run, _, _ in loaded[1:]:
        if not tc.same_setup(loaded[0][0]["meta"], run["meta"]):
            sys.exit("the runs use different reference frames or settings, so they cannot be shown together")
    ref = find_reference(loaded[0][0]["meta"]["reference"])
    allE = np.concatenate([d[k] for _, _, d in loaded for k in ("E_left", "E_right")])
    allN = np.concatenate([d[k] for _, _, d in loaded for k in ("N_left", "N_right")])
    origin = (float(np.floor(np.nanmin(allE) / 1000) * 1000), float(np.floor(np.nanmin(allN) / 1000) * 1000))
    runs = [json_run(run, c, d, PALETTE[i % len(PALETTE)], origin) for i, (run, c, d) in enumerate(loaded)]
    s_all = np.concatenate([d["s"] for _, _, d in loaded])
    s_range = (float(np.nanmin(s_all)), float(np.nanmax(s_all)))
    zs = np.concatenate([d[k] for _, _, d in loaded for k in ("Z_left", "Z_right")])
    z_rail = float(np.nanmedian(zs))
    gauge = np.concatenate([(d["d_left"] - d["d_right"]) * 1e3 for _, _, d in loaded])
    g_med = float(np.nanmedian(gauge))
    points = None
    if n_points > 0:
        first = loaded[0][0]["meta"]["input"]
        points = load_points(first["sha256"], ref, laz or first["file"], n_points, origin, s_range, z_rail)
    d_rails = float(np.nanmax(np.abs(np.concatenate([d[k] for _, _, d in loaded for k in ("d_left", "d_right")])))) + 0.1
    span = lambda lo, hi, pad: [round(float(lo) - pad, 4), round(float(hi) + pad, 4)]
    zpad = lambda lo, hi: max(0.02 * (hi - lo), 0.002)
    z_rails = span(np.nanmin(zs), np.nanmax(zs), zpad(np.nanmin(zs), np.nanmax(zs)))
    ranges = {"s": span(s_range[0], s_range[1], 0.0),
              "d_rails": [-round(d_rails, 3), round(d_rails, 3)], "z_rails": z_rails,
              "E_rails": span(np.nanmin(allE) - origin[0], np.nanmax(allE) - origin[0], 1.0),
              "N_rails": span(np.nanmin(allN) - origin[1], np.nanmax(allN) - origin[1], 1.0)}
    # the axes also cover the point cloud sample when it is shown, so each has a second range
    ranges.update({"d": ranges["d_rails"], "z": z_rails, "E": ranges["E_rails"], "N": ranges["N_rails"]})
    if points:
        pz = np.array(points["z"])
        lo, hi = min(z_rails[0], pz.min()), max(z_rails[1], pz.max())
        ranges["d"] = [-CROP_D, CROP_D]
        ranges["z"] = span(lo, hi, zpad(lo, hi))
        ranges["E"] = span(min(np.min(points["E"]), ranges["E_rails"][0] + 1.0), max(np.max(points["E"]), ranges["E_rails"][1] - 1.0), 1.0)
        ranges["N"] = span(min(np.min(points["N"]), ranges["N_rails"][0] + 1.0), max(np.max(points["N"]), ranges["N_rails"][1] - 1.0), 1.0)
    data = {"step": STEP, "origin": {"E": origin[0], "N": origin[1]},
            "reference": {"id": ref.id, "length": round(ref.length, 2)},
            "gauge_range": [round(g_med - GAUGE_HALF_RANGE, 2), round(g_med + GAUGE_HALF_RANGE, 2)],
            "ranges": ranges, "runs": runs, "points": points}
    return data, loaded


def render(data):
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "view_curve_template.html"), encoding="utf-8") as f:
        page = f.read()
    page = page.replace("/*DATA*/null", json.dumps(data, separators=(",", ":"), allow_nan=False), 1)
    return page.replace("/*PLOTLY*/", get_plotlyjs(), 1)


def main():
    ap = argparse.ArgumentParser(description="Interactive 3D view of the run curves, plus a LAS file for CloudCompare.")
    ap.add_argument("runs", nargs="+", help="run ids in output/runs/, built against the same reference")
    ap.add_argument("--points", type=int, default=60000, help="raw points shown as context, 0 for none (default 60000)")
    ap.add_argument("--laz", help="point cloud for the context points (default: the file the first run was built from)")
    ap.add_argument("--out", default=os.path.join(tc.OUT, "viewer"), help="output folder (default output/viewer)")
    ap.add_argument("--no-open", action="store_true", help="do not open the page in the browser")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    data, loaded = build_data(a.runs, a.points, a.laz)
    html = os.path.join(a.out, "curve3d_" + "_".join(a.runs) + ".html")
    with open(html, "w", encoding="utf-8", newline="\n") as f:
        f.write(render(data))
    print(f"3D view written to {html} ({os.path.getsize(html) / 1e6:.1f} MB)")
    for run_id, (_, _, dense) in zip(a.runs, loaded):
        path = os.path.join(a.out, f"{run_id}_curve.las")
        n = export_las(dense, path)
        print(f"CloudCompare file written to {path}: {n} points, scalar fields s_m, gauge_mm, rail")
    if not a.no_open:
        webbrowser.open("file:///" + os.path.abspath(html).replace(os.sep, "/"))


if __name__ == "__main__":
    main()
