"""The overhead contact wire as a curve: its height above the rail and its stagger, every metre along the track.

  python contact_wire.py run1 [--laz FILE]

The scan contains the overhead line: the contact wire about 5.5 m above the rails and the messenger wire above it.
Both quantities that matter for the pantograph are measured against the track in the same cross-section, so the
sensor-relative height of the scan cancels out, exactly as for the gauge:
  height    the contact wire above the plane of the two rail heads (m)
  stagger   its sideways position relative to the track centreline (m, + to the left), the zigzag that spreads the
            wear over the pantograph strip
Per 1 m bin on the run's control-point grid the lowest dense cluster of points above the track is taken as the contact
wire (droppers and the messenger wire lie higher); height and stagger are its medians. Bins where the wire points
spread sideways by more than 0.15 m are flagged: there two wires run side by side (an overlap) or a fitting is in the
way. A centripetal Catmull-Rom curve through the 1 m points gives the wire every 0.1 m for the 3D view.

Writes output/runs/<run>/contact_wire.csv (per metre) and contact_wire_dense.csv (every 0.1 m, with E, N, Z).
"""
import argparse
import os

import numpy as np

import catmull_rom as cr
import track_curve as tc

HALF_WIDTH = 1.3          # points within this sideways distance of the reference are considered (m)
Z_RANGE = (4.3, 6.6)      # height above the rail head where the contact wire can be (m)
BIN = 1.0                 # along-track bin, the run's control-point spacing (m)
MIN_POINTS = 6            # a bin needs this many points on the wire
CLUSTER = 0.02            # height histogram step for the lowest-cluster search (m)
BAND = (0.05, 0.07)       # wire points: this far below / above the cluster's bin (m)
DOUBLE_SPREAD = 0.15      # sideways spread above which two wires (or a fitting) are present (m)


def rail_plane(run):
    """Height of the rail-head plane along the track (mean of the two heads), interpolated over gaps."""
    sl = run["slabs"]
    z = 0.5 * (sl["z_left"] + sl["z_right"])
    ok = np.isfinite(z)
    return lambda s: np.interp(s, sl["s"][ok], z[ok])


def centreline(run):
    cv = run["curve"]
    ok = np.isfinite(cv["d_centre"])
    return lambda s: np.interp(s, cv["s"][ok], cv["d_centre"][ok])


MAX_STAGGER = 0.9         # the contact wire is never further than this from the track centre (m); side wires and fittings are
MAX_JUMP = 0.30           # the wire does not change height by more than this from one metre to the next (m)
PRIOR_WEIGHT = 0.3        # how fast the tracked position follows the latest metre


def clusters(h, dd, dc):
    """Height clusters among the points of one metre: contiguous populated 2 cm bins, each with its medians and size."""
    edges = np.arange(Z_RANGE[0], Z_RANGE[1] + CLUSTER, CLUSTER)
    hist, _ = np.histogram(h, bins=edges)
    pop = hist >= 3
    out, i = [], 0
    while i < len(pop):
        if not pop[i]:
            i += 1; continue
        j = i
        while j + 1 < len(pop) and pop[j + 1]:
            j += 1
        lo, hi = edges[i] - BAND[0], edges[j + 1] + BAND[1]
        sel = (h >= lo) & (h < hi)
        if sel.sum() >= MIN_POINTS:
            out.append({"h": float(np.median(h[sel])), "d": float(np.median(dd[sel])), "n": int(sel.sum()), "spread": float(np.std(dd[sel])), "lo": lo, "hi": hi})
        i = j + 1
    return out


def measure(pts, ref, run):
    s, d = ref.project(pts.E, pts.N)
    keep = (np.abs(d) < HALF_WIDTH) & (s >= 0) & (s <= ref.length)
    s, d, Z = s[keep], d[keep], pts.Z[keep]
    zr = Z - rail_plane(run)(s)
    keep = (zr > Z_RANGE[0]) & (zr < Z_RANGE[1])
    s, d, zr = s[keep], d[keep], zr[keep]
    centres = run["curve"]["s"]                                   # the 1 m control-point grid
    dc = centreline(run)(centres)
    cands = []
    for i, c in enumerate(centres):
        m = (s >= c - BIN / 2) & (s < c + BIN / 2)
        cands.append(clusters(zr[m], d[m], dc[i]) if m.sum() >= MIN_POINTS else [])
    # where to start: the typical height of the lowest cluster that sits near the track centre
    starts = [min(cl["h"] for cl in cs if abs(cl["d"] - dc[i]) < MAX_STAGGER) for i, cs in enumerate(cands) if any(abs(cl["d"] - dc[i]) < MAX_STAGGER for cl in cs)]
    prior_h = float(np.median(starts)) if starts else 5.6
    out = {k: np.full(len(centres), np.nan) for k in ("height", "stagger", "d_wire", "spread", "messenger")}
    n = np.zeros(len(centres), int)
    double = np.zeros(len(centres), bool)
    prior_d = None
    for i, cs in enumerate(cands):
        ok = [cl for cl in cs if abs(cl["d"] - dc[i]) < MAX_STAGGER and abs(cl["h"] - prior_h) < MAX_JUMP]
        if not ok:
            continue
        cost = lambda cl: abs(cl["h"] - prior_h) + (0.5 * abs(cl["d"] - dc[i] - prior_d) if prior_d is not None else 0)
        best = min(ok, key=cost)
        out["height"][i], out["d_wire"][i], out["spread"][i], n[i] = best["h"], best["d"], best["spread"], best["n"]
        out["stagger"][i] = best["d"] - dc[i]
        double[i] = best["spread"] > DOUBLE_SPREAD or any(cl is not best and abs(cl["h"] - best["h"]) < 0.15 for cl in cs)
        higher = [cl for cl in cs if cl["h"] > best["hi"] + 0.3]
        if higher:
            out["messenger"][i] = min(cl["h"] for cl in higher)
        prior_h = (1 - PRIOR_WEIGHT) * prior_h + PRIOR_WEIGHT * best["h"]
        prior_d = out["stagger"][i] if prior_d is None else (1 - PRIOR_WEIGHT) * prior_d + PRIOR_WEIGHT * out["stagger"][i]
    out["double"] = double
    return centres, out, n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_id")
    ap.add_argument("--laz", help="the point cloud the run was built from (default: the file named in the run)")
    a = ap.parse_args()
    folder = os.path.join(tc.OUT, "runs", a.run_id)
    run = tc.load_run(folder)
    ref = tc.load_reference(os.path.join(tc.OUT, "reference", run["meta"]["reference_name"] + ".json")) if "reference_name" in run["meta"] \
        else tc.load_reference(os.path.join(tc.OUT, "reference", "track.json"))
    pts = tc.read_points(a.laz or run["meta"]["input"]["file"])
    s, w, n = measure(pts, ref, run)
    ok = np.isfinite(w["height"])
    double = ok & w["double"]
    cols = ["s", "height_m", "stagger_m", "d_wire_m", "spread_m", "n_points", "double", "messenger_height_m"]
    with open(os.path.join(folder, "contact_wire.csv"), "w", encoding="utf-8", newline="\n") as f:
        f.write(",".join(cols) + "\n")
        for i in range(len(s)):
            vals = [f"{s[i]:.3f}"] + [("" if not np.isfinite(v) else f"{v:.4f}") for v in (w["height"][i], w["stagger"][i], w["d_wire"][i], w["spread"][i])] + \
                   [str(n[i]), "1" if double[i] else "0", "" if not np.isfinite(w["messenger"][i]) else f"{w['messenger'][i]:.3f}"]
            f.write(",".join(vals) + "\n")
    # the wire as a Catmull-Rom curve every 0.1 m, in the world frame, for the 3D view
    grid = np.arange(np.ceil(s[0] / 0.1) * 0.1, s[-1], 0.1)
    height = cr.evaluate(s, w["height"], grid, 3 * BIN)
    stagger = cr.evaluate(s, w["stagger"], grid, 3 * BIN)
    d_wire = cr.evaluate(s, w["d_wire"], grid, 3 * BIN)
    E, N = ref.to_world(grid, np.nan_to_num(d_wire))
    Z = rail_plane(run)(grid) + height
    good = np.isfinite(height) & np.isfinite(d_wire)
    with open(os.path.join(folder, "contact_wire_dense.csv"), "w", encoding="utf-8", newline="\n") as f:
        f.write("s,height_m,stagger_m,E,N,Z\n")
        for i in np.where(good)[0]:
            f.write(f"{grid[i]:.2f},{height[i]:.4f},{stagger[i]:.4f},{E[i]:.4f},{N[i]:.4f},{Z[i]:.4f}\n")
    print(f"contact wire found in {ok.sum()} of {len(s)} metres ({100 * ok.mean():.0f} %); {double.sum()} metres with two wires or a fitting")
    print(f"height above the rail: {np.nanmin(w['height']):.3f} to {np.nanmax(w['height']):.3f} m (median {np.nanmedian(w['height']):.3f})")
    print(f"stagger: {np.nanmin(w['stagger']):+.3f} to {np.nanmax(w['stagger']):+.3f} m (median {np.nanmedian(w['stagger']):+.3f})")
    if np.isfinite(w["messenger"]).any():
        print(f"messenger wire seen in {np.isfinite(w['messenger']).sum()} metres, {np.nanmin(w['messenger']):.2f} to {np.nanmax(w['messenger']):.2f} m above the rail")
    print(f"written to {folder}: contact_wire.csv, contact_wire_dense.csv")


if __name__ == "__main__":
    main()
