"""Per-point attributes of a run's 3D curve, and its short-wave height indicators.

  python curve_attributes.py run1 [--laz FILE] [--step 0.1]

For every point of the Catmull-Rom curve of each rail (every --step m), from the run and the point cloud it was built from:
  position      E, N, Z (as in curve_dense.csv), and s (chainage), d (sideways offset of the inner face)
  tangent       the direction of the curve: unit vector in the world frame (tx, ty, tz)
  head normal   the normal of the rail-head top surface, fitted (PCA plane) to the LiDAR points around the point;
                world frame (nx, ny, nz) and track frame (n_s along the track, n_d sideways, n_z up). n_d is the
                sideways tilt of the rail head: positive leans towards +d (the left). n_points is how many points
                the fit used; the normal is missing where there were too few.
  ground        ballast level between the rails next to the point (median height of the between-rail surface), and
                the nearest sleeper top (the periodic high points of that surface) with its distance along the track
  heights       the rail head above the ballast and above the nearest sleeper
Written to output/runs/<run>/curve_attributes.csv and, for CloudCompare, output/viewer/<run>_curve_attributes.las
(every attribute as a scalar field). The sleepers found go to output/runs/<run>/sleepers.csv. Once these files exist,
view_curve.py shows the attributes in the 3D view (colour modes, sleeper marks, values at the cursor).

Short-wave height indicators per 0.25 m slab (from slabs.csv): the mid-chord offset of each rail's height over a
3 m and a 10 m chord (surface), the height difference between the rails (cross-level) and its change over a 3 m
base (twist). Written to output/runs/<run>/shortwave.csv. The LiDAR height is relative to a sensor that rides on
the rails, so only the short-wave part of these is meaningful: a change the vehicle does not follow (a dip a metre
or two long) is seen; the slow part (cant, long-wave surface, the absolute cross-level) is cancelled. See the README.
"""
import argparse
import datetime
import os

import laspy
import numpy as np
from scipy.signal import find_peaks

import catmull_rom as cr
import track_curve as tc
import view_curve as vc

HEAD_HALF = 0.015         # half-width of the rail-head top used for the normal (m): the running band, not the crowned edges
NORMAL_ALONG = 0.15       # points within this distance along the track are used for the normal (m)
NORMAL_BAND = (0.012, 0.010)   # ... and from this far below to this far above the head top (m)
MIN_NORMAL_PTS = 20
GROUND_ZONE = (0.30, 0.55)     # between-rail surface: this far from the track centre (m), as in the pipeline
SLEEPER_ZONE = (0.20, 0.60)    # a wider zone for the sleeper tops
BELOW_HEAD = (0.08, 0.45)      # sleepers and ballast lie this far below the rail head (m): anything else between the rails is ignored
SLEEPER_BIN = 0.1              # bin along the track for the sleeper-top signal (m)
SLEEPER_MIN_SPACING = 0.45     # sleepers are at least this far apart (m)
SLEEPER_PROMINENCE = 0.015     # a sleeper top stands this much above the cribs (m)
LAS_DATE = datetime.date(2026, 1, 1)


def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(n > 0, n, 1)


def head_normals(s, d, z, s_pts, head_pts, z_pts):
    """PCA plane normal of the rail-head top around every curve point (track frame: along, sideways, up)."""
    out = np.full((len(s_pts), 3), np.nan)
    count = np.zeros(len(s_pts), int)
    lo = np.searchsorted(s, s_pts - NORMAL_ALONG)
    hi = np.searchsorted(s, s_pts + NORMAL_ALONG)
    for i, (a, b, h, zt) in enumerate(zip(lo, hi, head_pts, z_pts)):
        if not (np.isfinite(h) and np.isfinite(zt)) or b - a < MIN_NORMAL_PTS:
            continue
        ds, dd, dz = s[a:b], d[a:b], z[a:b]
        m = (np.abs(dd - h) <= HEAD_HALF) & (dz >= zt - NORMAL_BAND[0]) & (dz <= zt + NORMAL_BAND[1])
        if m.sum() < MIN_NORMAL_PTS:
            continue
        X = np.column_stack([ds[m], dd[m], dz[m]])
        X = X - X.mean(0)
        _, _, vt = np.linalg.svd(X, full_matrices=False)
        n = vt[-1]
        out[i] = n if n[2] >= 0 else -n
        count[i] = int(m.sum())
    return out, count


def ground_and_sleepers(s, d, z, centre_s, centre_d, centre_z, length):
    """Ballast level per 0.25 m along the track, and the sleeper tops found in the between-rail surface."""
    dc = np.interp(s, centre_s, centre_d)                 # track centre at every point's chainage
    zc = np.interp(s, centre_s, centre_z)                 # rail-head height there
    off = np.abs(d - dc)
    below = (z < zc - BELOW_HEAD[0]) & (z > zc - BELOW_HEAD[1])   # the between-rail surface lies below the rail head
    bins = np.arange(0, length + SLEEPER_BIN, SLEEPER_BIN)
    centres = 0.5 * (bins[:-1] + bins[1:])
    top = np.full(len(centres), np.nan)
    ground = np.full(len(centres), np.nan)
    zone = (off > SLEEPER_ZONE[0]) & (off < SLEEPER_ZONE[1]) & below
    gz = (off > GROUND_ZONE[0]) & (off < GROUND_ZONE[1]) & below
    idx_zone = np.flatnonzero(zone)
    idx_gz = np.flatnonzero(gz)
    b_zone = np.searchsorted(s[idx_zone], bins)
    b_gz = np.searchsorted(s[idx_gz], bins)
    for i in range(len(centres)):
        pz = z[idx_zone[b_zone[i]:b_zone[i + 1]]]
        if len(pz) >= 15:
            top[i] = np.percentile(pz, 90)
        pg = z[idx_gz[max(0, b_gz[max(0, i - 1)]):b_gz[min(len(bins) - 1, i + 2)]]]   # 0.3 m window, as the pipeline's slab
        if len(pg) >= 20:
            ground[i] = np.median(pg)
    ok = np.isfinite(top)
    filled = np.interp(centres, centres[ok], top[ok]) if ok.sum() > 2 else top
    peaks, props = find_peaks(filled, prominence=SLEEPER_PROMINENCE, distance=int(round(SLEEPER_MIN_SPACING / SLEEPER_BIN)))
    peaks = peaks[ok[peaks]]
    return centres, ground, centres[peaks], top[peaks]


def shortwave(run):
    """Surface (mid-chord offsets of each rail's height), cross-level and twist per slab, in mm."""
    sl = run["slabs"]
    s, zl, zr = sl["s"], sl["z_left"], sl["z_right"]
    half3, half10 = int(round(1.5 / tc.P.slab)), int(round(5.0 / tc.P.slab))
    cross = zl - zr
    twist = np.full(len(s), np.nan)
    twist[half3:-half3] = cross[2 * half3:] - cross[:-2 * half3]
    return {"s": s, "surface_3m_left": tc._midchord(zl, half3) * 1e3, "surface_3m_right": tc._midchord(zr, half3) * 1e3,
            "surface_10m_left": tc._midchord(zl, half10) * 1e3, "surface_10m_right": tc._midchord(zr, half10) * 1e3,
            "cross_level": cross * 1e3, "twist_3m": twist * 1e3}


def attributes(run, ref, pts, step):
    """All per-point attributes of both rails as a dict of arrays, plus the sleepers found."""
    c = run["curve"]
    dense = tc.sample_curve(c, step)
    grid = dense["s"]
    s_all, d_all = ref.project(pts.E, pts.N)
    keep = (np.abs(d_all) < tc.P.d_max) & (s_all >= 0) & (s_all <= ref.length)
    order = np.argsort(s_all[keep], kind="stable")
    s, d, z = s_all[keep][order], d_all[keep][order], pts.Z[keep][order]
    sl = run["slabs"]
    centre_ok = np.isfinite(sl["head_left"]) & np.isfinite(sl["head_right"]) & np.isfinite(sl["z_left"]) & np.isfinite(sl["z_right"])
    centre_d = 0.5 * (sl["head_left"] + sl["head_right"])
    centre_z = 0.5 * (sl["z_left"] + sl["z_right"])
    bin_s, ground, sleeper_s, sleeper_z = ground_and_sleepers(s, d, z, sl["s"][centre_ok], centre_d[centre_ok], centre_z[centre_ok], ref.length)
    tx_ref, ty_ref = np.interp(grid, ref.s, ref.tx), np.interp(grid, ref.s, ref.ty)
    nx_ref, ny_ref = np.interp(grid, ref.s, ref.nx), np.interp(grid, ref.s, ref.ny)
    out = {"s": grid}
    for side, rail in (("left", 0), ("right", 1)):
        head = np.interp(grid, sl["s"][np.isfinite(sl["head_" + side])], sl["head_" + side][np.isfinite(sl["head_" + side])])
        E, N, Z = dense["E_" + side], dense["N_" + side], dense["Z_" + side]
        dd = cr.evaluate(c["s"], c["d_" + side], grid, 3 * tc.P.ctrl)
        tangent = unit(np.column_stack([np.gradient(E, grid), np.gradient(N, grid), np.gradient(Z, grid)]))
        n_track, n_count = head_normals(s, d, z, grid, head, Z)
        n_world = np.column_stack([n_track[:, 0] * tx_ref + n_track[:, 1] * nx_ref, n_track[:, 0] * ty_ref + n_track[:, 1] * ny_ref, n_track[:, 2]])
        g = np.interp(grid, bin_s[np.isfinite(ground)], ground[np.isfinite(ground)]) if np.isfinite(ground).any() else np.full(len(grid), np.nan)
        j = np.clip(np.searchsorted(sleeper_s, grid), 1, len(sleeper_s) - 1) if len(sleeper_s) > 1 else np.zeros(len(grid), int)
        nearer = np.where(np.abs(grid - sleeper_s[j - 1]) < np.abs(grid - sleeper_s[j]), j - 1, j) if len(sleeper_s) > 1 else j
        sz = sleeper_z[nearer] if len(sleeper_s) else np.full(len(grid), np.nan)
        sd = (grid - sleeper_s[nearer]) if len(sleeper_s) else np.full(len(grid), np.nan)
        for name, v in (("E", E), ("N", N), ("Z", Z), ("d", dd), ("tx", tangent[:, 0]), ("ty", tangent[:, 1]), ("tz", tangent[:, 2]),
                        ("nx", n_world[:, 0]), ("ny", n_world[:, 1]), ("nz", n_world[:, 2]),
                        ("n_s", n_track[:, 0]), ("n_d", n_track[:, 1]), ("n_z", n_track[:, 2]), ("n_points", n_count.astype(float)),
                        ("ground_z", g), ("sleeper_z", sz), ("sleeper_dist", sd), ("above_ground", Z - g), ("above_sleeper", Z - sz)):
            out[f"{name}_{side}"] = v
    return out, sleeper_s, sleeper_z


def _decimals(col):
    if col.startswith("n_points"):
        return 0
    if col == "s" or col.startswith(("surface", "cross", "twist")):
        return 3
    return 4


def write_csv(path, cols, data):
    with open(path, "w", newline="\n") as f:
        f.write(",".join(cols) + "\n")
        for i in range(len(data[cols[0]])):
            f.write(",".join("" if not np.isfinite(data[c][i]) else f"{data[c][i]:.{_decimals(c)}f}" for c in cols) + "\n")


def write_las(path, attrs):
    parts = []
    for rail, side in enumerate(("left", "right")):
        m = np.isfinite(attrs["E_" + side]) & np.isfinite(attrs["Z_" + side])
        parts.append({k: attrs[f"{k}_{side}"][m] for k in ("E", "N", "Z", "d", "tx", "ty", "tz", "nx", "ny", "nz", "n_s", "n_d", "n_z",
                                                          "ground_z", "sleeper_z", "sleeper_dist", "above_ground", "above_sleeper")}
                     | {"s": attrs["s"][m], "rail": np.full(m.sum(), rail, np.uint8)})
    cat = lambda k: np.concatenate([p[k] for p in parts])
    hdr = laspy.LasHeader(point_format=2, version="1.2")
    hdr.creation_date = LAS_DATE
    hdr.scales = np.array([1e-4, 1e-4, 1e-4])
    hdr.offsets = np.floor([cat("E").min(), cat("N").min(), cat("Z").min()])
    fields = [("s_m", "f4"), ("d_m", "f4"), ("rail", "u1"), ("tx", "f4"), ("ty", "f4"), ("tz", "f4"), ("nx", "f4"), ("ny", "f4"), ("nz", "f4"),
              ("n_s", "f4"), ("n_d", "f4"), ("n_z", "f4"), ("ground_z", "f4"), ("sleeper_z", "f4"), ("sleeper_dist", "f4"),
              ("above_ground", "f4"), ("above_sleeper", "f4")]
    for name, typ in fields:
        hdr.add_extra_dim(laspy.ExtraBytesParams(name=name, type=typ))
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = cat("E"), cat("N"), cat("Z")
    src = {"s_m": "s", "d_m": "d"}
    for name, typ in fields:
        setattr(las, name, cat(src.get(name, name)).astype(np.float32 if typ == "f4" else np.uint8))
    las.write(path)
    return len(las.x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_id")
    ap.add_argument("--laz", help="the point cloud the run was built from (default: the file named in the run)")
    ap.add_argument("--step", type=float, default=0.1, help="spacing of the curve points (m, default 0.1)")
    a = ap.parse_args()
    folder = os.path.join(tc.OUT, "runs", a.run_id)
    run = tc.load_run(folder)
    ref = vc.find_reference(run["meta"]["reference"])
    pts = tc.read_points(a.laz or run["meta"]["input"]["file"])
    if pts.sha256 != run["meta"]["input"]["sha256"] and "simulation" not in run["meta"]:
        raise SystemExit("the point cloud is not the file this run was built from")

    attrs, sleeper_s, sleeper_z = attributes(run, ref, pts, a.step)
    cols = ["s"] + [f"{k}_{side}" for side in ("left", "right") for k in
                    ("E", "N", "Z", "d", "tx", "ty", "tz", "nx", "ny", "nz", "n_s", "n_d", "n_z", "n_points", "ground_z", "sleeper_z",
                     "sleeper_dist", "above_ground", "above_sleeper")]
    write_csv(os.path.join(folder, "curve_attributes.csv"), cols, attrs)
    write_csv(os.path.join(folder, "sleepers.csv"), ["s", "z"], {"s": sleeper_s, "z": sleeper_z})
    sw = shortwave(run)
    write_csv(os.path.join(folder, "shortwave.csv"), list(sw), sw)
    os.makedirs(os.path.join(tc.OUT, "viewer"), exist_ok=True)
    n = write_las(os.path.join(tc.OUT, "viewer", f"{a.run_id}_curve_attributes.las"), attrs)

    spacing = np.diff(sleeper_s)
    print(f"{a.run_id}: {len(attrs['s'])} curve points per rail every {a.step} m")
    print(f"sleepers found: {len(sleeper_s)}, spacing median {np.median(spacing):.2f} m (IQR {np.percentile(spacing, 25):.2f}–{np.percentile(spacing, 75):.2f})"
          if len(spacing) else "no sleepers found")
    for side in ("left", "right"):
        ok = np.isfinite(attrs[f"n_d_{side}"])
        print(f"{side:5s} rail: head normal fitted at {ok.mean():.0%} of the points; sideways tilt n_d median {np.nanmedian(attrs[f'n_d_{side}']):+.4f} "
              f"(1:20 inclination would be about {'-' if side == 'left' else '+'}0.050); "
              f"head above ballast {np.nanmedian(attrs[f'above_ground_{side}']) * 1e3:.0f} mm, above nearest sleeper "
              f"{np.nanmedian(attrs[f'above_sleeper_{side}']) * 1e3:.0f} mm (median)")
    for k in ("surface_3m_left", "surface_3m_right", "cross_level", "twist_3m"):
        v = sw[k][np.isfinite(sw[k])]
        print(f"{k:17s} std {v.std():.2f} mm, 95 % within ±{np.percentile(np.abs(v - np.median(v)), 95):.2f} mm of the median {np.median(v):+.2f} mm")
    print(f"written {folder}\\curve_attributes.csv, {folder}\\shortwave.csv and output\\viewer\\{a.run_id}_curve_attributes.las ({n} points)")


if __name__ == "__main__":
    main()
