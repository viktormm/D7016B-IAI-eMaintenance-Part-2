"""One feature vector per metre of track, from what the pipeline already measured, for the unsupervised track-state view.

  python track_features.py run1 [--laz FILE]

Per 1 m control point of the run (output/runs/<run>/features.csv):
  gauge_mm           gauge at the control point, as deviation from 1435 mm            (curve.csv)
  versine_mm         lateral alignment: the centreline's offset from the chord between the points 5 m before and after (curve.csv)
  above_ground_mm    rail head above the ballast between the rails, mean of both rails  (curve_attributes.csv, if present)
  above_sleeper_mm   rail head above the nearest sleeper top, mean of both rails         (curve_attributes.csv, if present)
  head_tilt          sideways tilt of the rail-head normal, mean of both rails           (curve_attributes.csv, if present)
  sleeper_spacing_m  spacing of the sleepers around this metre                            (sleepers.csv, if present)
  points_per_m       LiDAR points within 2.5 m of the track in this metre: point density (slabs.csv)
  intensity_left/right  median LiDAR intensity of the points on each rail head in this metre (the point cloud)
  wire_height_m, wire_stagger_m  the contact wire above the rail and its stagger          (contact_wire.csv, if present)
The clustering itself runs in the page (monitor.py), where the number of states can be changed.
"""
import argparse
import os
import warnings

import numpy as np

import track_curve as tc

warnings.filterwarnings("ignore", category=RuntimeWarning)       # means of empty metres are NaN by design


def per_metre_mean(s_from, v, centres, half=0.5):
    out = np.full(len(centres), np.nan)
    ok = np.isfinite(v)
    for i, c in enumerate(centres):
        m = ok & (s_from >= c - half) & (s_from < c + half)
        if m.sum():
            out[i] = np.mean(v[m])
    return out


def columns(path):
    if not os.path.exists(path):
        return None
    t = np.genfromtxt(path, delimiter=",", names=True)
    return {k: np.asarray(t[k], float) for k in t.dtype.names}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_id")
    ap.add_argument("--laz", help="the point cloud the run was built from (default: the file named in the run)")
    a = ap.parse_args()
    folder = os.path.join(tc.OUT, "runs", a.run_id)
    run = tc.load_run(folder)
    cv, sl = run["curve"], run["slabs"]
    s = cv["s"]
    F = {"s": s, "gauge_mm": (cv["gauge"] - tc.NOMINAL_GAUGE) * 1e3}
    # lateral alignment as a 10 m versine of the centreline
    dc = cv["d_centre"]
    ok = np.isfinite(dc)
    dci = np.interp(s, s[ok], dc[ok])
    vers = np.full(len(s), np.nan)
    for i in range(len(s)):
        if s[i] - 5 >= s[0] and s[i] + 5 <= s[-1]:
            vers[i] = dci[i] - 0.5 * (np.interp(s[i] - 5, s, dci) + np.interp(s[i] + 5, s, dci))
    F["versine_mm"] = vers * 1e3
    # point density per metre near the track
    F["points_per_m"] = np.array([sl["points"][(sl["s"] >= c - 0.5) & (sl["s"] < c + 0.5)].sum() for c in s], float)
    # per-point attributes, averaged per metre
    at = columns(os.path.join(folder, "curve_attributes.csv"))
    if at:
        for name, cols in (("above_ground_mm", ("above_ground_left", "above_ground_right")), ("above_sleeper_mm", ("above_sleeper_left", "above_sleeper_right"))):
            F[name] = 1e3 * np.nanmean([per_metre_mean(at["s"], at[c], s) for c in cols], axis=0)
        F["head_tilt"] = np.nanmean([per_metre_mean(at["s"], at[c], s) for c in ("n_d_left", "n_d_right")], axis=0)
    sp = columns(os.path.join(folder, "sleepers.csv"))
    if sp and len(np.atleast_1d(sp["s"])) > 2:
        ss = np.sort(np.atleast_1d(sp["s"]))
        mids, gaps = 0.5 * (ss[1:] + ss[:-1]), np.diff(ss)
        F["sleeper_spacing_m"] = np.array([np.mean(gaps[(mids >= c - 0.5) & (mids < c + 0.5)]) if ((mids >= c - 0.5) & (mids < c + 0.5)).any() else np.nan for c in s])
    cw = columns(os.path.join(folder, "contact_wire.csv"))
    if cw:
        F["wire_height_m"] = np.interp(s, cw["s"], cw["height_m"], left=np.nan, right=np.nan) if np.isfinite(cw["height_m"]).all() else \
            np.array([cw["height_m"][np.argmin(np.abs(cw["s"] - c))] if np.min(np.abs(cw["s"] - c)) < 0.01 else np.nan for c in s])
        F["wire_stagger_m"] = np.array([cw["stagger_m"][np.argmin(np.abs(cw["s"] - c))] if np.min(np.abs(cw["s"] - c)) < 0.01 else np.nan for c in s])
    # rail-head intensity per metre and rail, from the point cloud
    pts_path = a.laz or run["meta"]["input"]["file"]
    import laspy
    las = laspy.read(pts_path)
    ref = tc.load_reference(os.path.join(tc.OUT, "reference", "track.json"))
    E, N, Z, I = np.asarray(las.x), np.asarray(las.y), np.asarray(las.z), np.asarray(las.intensity, float)
    sp_, dp = ref.project(E, N)
    keep = (np.abs(dp) < 1.2) & (sp_ >= s[0] - 0.5) & (sp_ <= s[-1] + 0.5)
    sp_, dp, Z, I = sp_[keep], dp[keep], Z[keep], I[keep]
    for side in ("left", "right"):
        okf = np.isfinite(cv["d_" + side]) & np.isfinite(cv["z_" + side])
        d_face = np.interp(sp_, s[okf], cv["d_" + side][okf]); z_top = np.interp(sp_, s[okf], cv["z_" + side][okf])
        on_head = (np.abs(dp - d_face) < 0.08) & (np.abs(Z - z_top) < 0.02)
        F["intensity_" + side] = np.array([np.median(I[on_head & (sp_ >= c - 0.5) & (sp_ < c + 0.5)]) if (on_head & (sp_ >= c - 0.5) & (sp_ < c + 0.5)).sum() >= 10 else np.nan for c in s])
    cols = list(F)
    with open(os.path.join(folder, "features.csv"), "w", encoding="utf-8", newline="\n") as f:
        f.write(",".join(cols) + "\n")
        for i in range(len(s)):
            f.write(",".join("" if not np.isfinite(F[c][i]) else (f"{F[c][i]:.3f}" if c != "s" else f"{F[c][i]:.1f}") for c in cols) + "\n")
    complete = np.all([np.isfinite(F[c]) for c in cols if c != "s"], axis=0).sum()
    print(f"features written for {len(s)} metres ({complete} with every feature): " + ", ".join(c for c in cols if c != "s"))


if __name__ == "__main__":
    main()
