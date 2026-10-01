"""Repeatable pipeline: LiDAR run (.laz) -> rail curves in a fixed reference frame -> comparison of runs.

Why it is repeatable
  * One reference frame (the smoothed track centreline of the first run) is defined once and stored. Every run,
    including the first, is projected onto it: chainage s along the reference, sideways offset d from it, and the
    height Z exactly as stored in the LiDAR file. Left is +d, so left/right do not depend on the travel direction.
  * Slabs, control points and sections sit on fixed grids anchored at s = 0, so runs line up point for point.
  * Nothing is random and no parameter is chosen from the data. Parameters, the input file hash, the reference id
    and a hash of the code are stored with every run, and `compare` refuses runs made with different settings.

Commands (run from the repo root):
  python track_curve.py reference RUN.laz --name track
  python track_curve.py build RUN.laz --run-id run1 --reference output/reference/track.json
  python track_curve.py sample run1 [--step 0.1]
  python track_curve.py sections run1 [--norm-run run1]
  python track_curve.py compare run1 run2 [--plot]
  python track_curve.py selftest RUN.laz --reference output/reference/track.json
"""
import argparse
import hashlib
import json
import os
import sys
from dataclasses import dataclass, asdict

import numpy as np
import laspy
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.spatial import cKDTree

import catmull_rom as cr

OUT = "output"
CODE_HASH = hashlib.sha256(open(__file__, "rb").read() + open(cr.__file__, "rb").read()).hexdigest()[:12]
NOMINAL_GAUGE = 1.435


@dataclass(frozen=True)
class Params:
    slab: float = 0.25            # slab length along the reference (m)
    d_max: float = 2.5            # keep points within this offset of the reference (m)
    win: float = 0.12             # half-width of the rail search window (m)
    head_min: float = 0.10        # rail-head candidates: height above the local ballast (m)
    head_max: float = 0.30
    top_band: float = 0.008       # head-top points lie within this of the top (m)
    min_pts: int = 8
    outlier: float = 0.02         # reject a slab this far from the running median (m)
    face_outlier: float = 0.01    # reject a face position this far from the running median (m)
    face_below: tuple = (0.020, 0.008)   # inner face sampled this far below the head top (m); gauge is defined at 14 mm
    face_win: float = 0.08
    min_face: int = 6
    ctrl: float = 1.0             # control point spacing (m)
    section: float = 25.0         # section length for the condition table (m)
    chord: float = 10.0           # chord for the alignment measure (m)
    ref_sigma: float = 25.0       # Gaussian smoothing of the reference centreline (m)
    ref_step: float = 0.25        # reference sample spacing (m)
    ref_margin: float = 2.0       # the reference extends this far beyond the first run's rails (m)


P = Params()


# ------------------------------------------------------------------ input


@dataclass
class Points:
    E: np.ndarray
    N: np.ndarray
    Z: np.ndarray
    name: str
    sha256: str


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_points(path):
    las = laspy.read(path)
    return Points(np.asarray(las.x, float), np.asarray(las.y, float), np.asarray(las.z, float),
                  os.path.basename(path), sha256_file(path))


# ------------------------------------------------------------------ rail extraction on a fixed slab grid


def _ballast_z(ds, zs, dc):
    m = (np.abs(ds - dc) > 0.3) & (np.abs(ds - dc) < 0.55)
    return np.median(zs[m]) if m.sum() > 20 else np.nan


def _find_top(ds, zs, d_exp, zb, p):
    m = (np.abs(ds - d_exp) < p.win) & (zs > zb + p.head_min) & (zs < zb + p.head_max)
    if m.sum() < p.min_pts:
        return None
    dd, zz = ds[m], zs[m]
    top = zz >= np.percentile(zz, 97) - p.top_band
    if top.sum() < p.min_pts // 2:
        return None
    return np.median(dd[top]), np.median(zz[top]), int(top.sum())


def _initial_rails(s, d, z, mid, p):
    m = np.abs(s - mid) < 5
    ds, zs = d[m], z[m]
    near = np.abs(ds) < 0.3
    dc0 = np.median(ds[near]) if near.any() else 0.0
    zb = np.median(zs[(np.abs(ds - dc0) > 0.3) & (np.abs(ds - dc0) < 0.55)])
    hi = (zs > zb + 0.12) & (zs < zb + p.head_max) & (np.abs(ds - dc0) < 1.5)
    h, e = np.histogram(ds[hi] - dc0, bins=np.arange(-1.5, 1.5001, 0.01))
    h = np.convolve(h, np.ones(3) / 3, mode="same")
    c = 0.5 * (e[:-1] + e[1:])
    return dc0 + c[c < 0][np.argmax(h[c < 0])], dc0 + c[c > 0][np.argmax(h[c > 0])]


def track_slabs(s, d, z, edges, p):
    """Rail-head top and inner face per slab. s must be sorted. Returns arrays with NaN where a rail is not found."""
    n = len(edges) - 1
    bounds = np.searchsorted(s, edges)
    covered = np.flatnonzero(np.diff(bounds) >= 200)
    if len(covered) == 0:
        raise RuntimeError("no slab has enough points: does the run cover the reference?")
    mid_i = covered[len(covered) // 2]
    dR0, dL0 = _initial_rails(s, d, z, 0.5 * (edges[mid_i] + edges[mid_i + 1]), p)
    out = {k: np.full(n, np.nan) for k in ("dR", "zR", "nR", "dL", "zL", "nL", "fR", "fL")}
    out["points"] = np.diff(bounds).astype(float)

    def walk(indices, dR, dL):
        for i in indices:
            ds, zs = d[bounds[i]:bounds[i + 1]], z[bounds[i]:bounds[i + 1]]
            zb = _ballast_z(ds, zs, 0.5 * (dR + dL))
            if np.isnan(zb):
                continue
            for r, dexp, pct in (("R", dR, 95), ("L", dL, 5)):
                top = _find_top(ds, zs, dexp, zb, p)
                if top is None:
                    continue
                out["d" + r][i], out["z" + r][i], out["n" + r][i] = top
                if r == "R":
                    dR = top[0]
                else:
                    dL = top[0]
                k = (np.abs(ds - top[0]) < p.face_win) & (zs < top[1] - p.face_below[1]) & (zs > top[1] - p.face_below[0])
                if k.sum() >= p.min_face:
                    out["f" + r][i] = np.percentile(ds[k], pct)   # inner face: the side towards the track centre

    walk(range(mid_i, n), dR0, dL0)
    walk(range(mid_i - 1, -1, -1), dR0, dL0)
    return out


def _mask_series(sc, v, thr):
    """v with NaN where it is further than thr from the running median of its neighbours (9 slabs)."""
    ok = np.isfinite(v)
    if ok.sum() < 3:
        return v
    fill = np.interp(sc, sc[ok], v[ok])
    out = v.copy()
    with np.errstate(invalid="ignore"):
        out[np.abs(v - median_filter(fill, size=9, mode="nearest")) > thr] = np.nan
    return out


def _reject_outliers(sc, d, z, p):
    ok = np.isfinite(d)
    if ok.sum() < 3:
        return ok
    for v in (d, z):
        fill = np.interp(sc, sc[ok], v[ok])
        with np.errstate(invalid="ignore"):
            ok &= ~(np.abs(v - median_filter(fill, size=9, mode="nearest")) > p.outlier)
    return ok


# ------------------------------------------------------------------ reference frame


class Reference:
    """The fixed frame every run is expressed in: s along the reference centreline, d sideways (left is +d)."""

    def __init__(self, data):
        self.data = data
        self.s, self.E, self.N = (np.asarray(data[k], float) for k in ("s", "E", "N"))
        tx, ty = np.gradient(self.E), np.gradient(self.N)
        norm = np.hypot(tx, ty)
        self.tx, self.ty = tx / norm, ty / norm
        self.nx, self.ny = -self.ty, self.tx
        self.length = float(self.s[-1])
        self.id = data["id"]
        self.tree = cKDTree(np.column_stack([self.E, self.N]))

    def project(self, E, N):
        _, i = self.tree.query(np.column_stack([E, N]))
        vx, vy = E - self.E[i], N - self.N[i]
        return self.s[i] + vx * self.tx[i] + vy * self.ty[i], vx * self.nx[i] + vy * self.ny[i]

    def normal_at(self, s):
        return np.interp(s, self.s, self.nx), np.interp(s, self.s, self.ny)

    def to_world(self, s, d):
        f = lambda a: np.interp(s, self.s, a)
        return f(self.E) + d * f(self.nx), f(self.N) + d * f(self.ny)


def load_reference(path):
    return Reference(json.load(open(path)))


def make_reference(pts, p, name):
    xy = np.column_stack([pts.E[::20], pts.N[::20]])
    c = xy.mean(0)
    w, v = np.linalg.eigh(np.cov((xy - c).T))
    u = v[:, np.argmax(w)]
    if u[0] < 0:
        u = -u
    n = np.array([-u[1], u[0]])
    x, y = pts.E - c[0], pts.N - c[1]
    s0, d0 = x * u[0] + y * u[1], x * n[0] + y * n[1]
    keep = np.abs(d0) < p.d_max
    order = np.argsort(s0[keep], kind="stable")
    s, d, z = s0[keep][order], d0[keep][order], pts.Z[keep][order]
    edges = np.arange(np.floor(s[0] / p.slab) * p.slab, s[-1] + p.slab, p.slab)
    sl = track_slabs(s, d, z, edges, p)
    sc = 0.5 * (edges[:-1] + edges[1:])
    dc = 0.5 * (sl["dL"] + sl["dR"])
    ok = np.isfinite(dc)
    dc = np.interp(sc, sc[ok], dc[ok])
    sigma = p.ref_sigma / p.slab
    pad = int(4 * sigma) + 1
    smooth = gaussian_filter1d(np.pad(dc, pad, mode="reflect", reflect_type="odd"), sigma, mode="nearest")[pad:-pad]
    k = int(round(5.0 / p.slab))
    slope0, slope1 = np.polyfit(sc[:k], smooth[:k], 1)[0], np.polyfit(sc[-k:], smooth[-k:], 1)[0]
    m = p.ref_margin
    se = np.concatenate([[sc[0] - m], sc, [sc[-1] + m]])
    de = np.concatenate([[smooth[0] - slope0 * m], smooth, [smooth[-1] + slope1 * m]])
    E, N = c[0] + se * u[0] + de * n[0], c[1] + se * u[1] + de * n[1]
    arc = np.concatenate([[0], np.cumsum(np.hypot(np.diff(E), np.diff(N)))])
    s_ref = np.arange(0, arc[-1], p.ref_step)
    E_ref, N_ref = np.round(np.interp(s_ref, arc, E), 4), np.round(np.interp(s_ref, arc, N), 4)
    rid = hashlib.sha256(np.round(np.column_stack([E_ref, N_ref]), 3).tobytes()).hexdigest()[:12]
    return {"name": name, "id": rid, "source": pts.name, "source_sha256": pts.sha256, "code": CODE_HASH,
            "params": asdict(p), "s": np.round(s_ref, 4).tolist(), "E": E_ref.tolist(), "N": N_ref.tolist()}


# ------------------------------------------------------------------ run -> curve


def _bin_value(sc, v, centre, width):
    m = (sc >= centre - width / 2) & (sc < centre + width / 2) & np.isfinite(v)
    k = int(m.sum())
    if k == 0:
        return np.nan
    if k >= 3:
        a, b = np.polyfit(sc[m] - centre, v[m], 1)
        return b
    return float(v[m].mean())


def build_run(pts, ref, p, run_id, mask=None):
    E, N, Z = (pts.E, pts.N, pts.Z) if mask is None else (pts.E[mask], pts.N[mask], pts.Z[mask])
    s, d = ref.project(E, N)
    keep = (np.abs(d) < p.d_max) & (s >= 0) & (s <= ref.length)
    s, d, z = s[keep], d[keep], Z[keep]
    order = np.argsort(s, kind="stable")
    s, d, z = s[order], d[order], z[order]
    edges = np.arange(0, ref.length + p.slab, p.slab)
    edges = edges[edges <= ref.length + 1e-9]
    sl = track_slabs(s, d, z, edges, p)
    sc = 0.5 * (edges[:-1] + edges[1:])
    okR = _reject_outliers(sc, sl["dR"], sl["zR"], p)
    okL = _reject_outliers(sc, sl["dL"], sl["zL"], p)
    faceR = _mask_series(sc, np.where(okR, sl["fR"], np.nan), p.face_outlier)
    faceL = _mask_series(sc, np.where(okL, sl["fL"], np.nan), p.face_outlier)
    both = np.isfinite(faceR) & np.isfinite(faceL) & okR & okL
    gauge = np.where(both, faceL - faceR, np.nan)
    slabs = {"s": sc, "d_right": faceR, "z_right": np.where(okR, sl["zR"], np.nan),
             "d_left": faceL, "z_left": np.where(okL, sl["zL"], np.nan), "gauge": gauge, "points": sl["points"],
             "head_right": np.where(okR, sl["dR"], np.nan), "head_left": np.where(okL, sl["dL"], np.nan),
             "n_head_right": sl["nR"], "n_head_left": sl["nL"]}

    cs = (np.arange(int(ref.length // p.ctrl)) + 0.5) * p.ctrl
    curve = {"s": cs}
    for name in ("d_left", "z_left", "d_right", "z_right", "gauge"):
        curve[name] = np.array([_bin_value(sc, slabs[name], c, p.ctrl) for c in cs])
    curve["n_slabs"] = np.array([float(np.sum(both & (sc >= c - p.ctrl / 2) & (sc < c + p.ctrl / 2))) for c in cs])
    curve["d_centre"] = 0.5 * (curve["d_left"] + curve["d_right"])
    curve["z_centre"] = 0.5 * (curve["z_left"] + curve["z_right"])
    for side in ("left", "right"):
        curve["E_" + side], curve["N_" + side] = ref.to_world(cs, np.nan_to_num(curve["d_" + side], nan=0.0))
        bad = ~np.isfinite(curve["d_" + side])
        curve["E_" + side][bad] = np.nan
        curve["N_" + side][bad] = np.nan

    gv = gauge[np.isfinite(gauge)]
    qa = {"slabs": int(len(sc)), "both_rails_found": round(float(np.mean(both)), 4),
          "rejected_outliers": int(np.isfinite(sl["dR"]).sum() - okR.sum() + np.isfinite(sl["dL"]).sum() - okL.sum()),
          "gauge_mean_mm": round(float(gv.mean() * 1e3), 2), "gauge_std_mm": round(float(gv.std() * 1e3), 2),
          "control_points": int(len(cs)), "control_points_with_both_rails": int(np.sum(np.isfinite(curve["d_centre"]))),
          "points_used": int(len(s))}
    meta = {"run_id": run_id, "input": {"file": pts.name, "sha256": pts.sha256, "points": int(len(pts.E) if mask is None else mask.sum())},
            "reference": ref.id, "code": CODE_HASH, "params": asdict(p), "qa": qa,
            "columns": {"s": "chainage along the reference (m)", "d_*": "sideways offset of the rail's inner (gauge) face from the reference, left is +d (m)",
                        "z_*": "height Z of the rail-head top as stored in the LiDAR file (m), relative to the sensor",
                        "gauge": "distance between the inner rail faces, d_left - d_right (m)"}}
    return {"slabs": slabs, "curve": curve, "meta": meta}


CURVE_COLS = ["s", "n_slabs", "d_left", "z_left", "d_right", "z_right", "gauge", "d_centre", "z_centre",
              "E_left", "N_left", "E_right", "N_right"]
SLAB_COLS = ["s", "points", "d_left", "z_left", "d_right", "z_right", "gauge", "head_left", "head_right", "n_head_left", "n_head_right"]


def _table_bytes(cols, data, fmt_by_col):
    header = ",".join(cols)
    lines = [header]
    for i in range(len(data[cols[0]])):
        lines.append(",".join(fmt_by_col.get(c, "%.5f") % data[c][i] for c in cols))
    return ("\n".join(lines) + "\n").encode()


def serialise(run):
    fm = {"n_slabs": "%.0f", "points": "%.0f", "n_head_left": "%.0f", "n_head_right": "%.0f", "E_left": "%.4f", "N_left": "%.4f",
          "E_right": "%.4f", "N_right": "%.4f"}
    return {"curve.csv": _table_bytes(CURVE_COLS, run["curve"], fm), "slabs.csv": _table_bytes(SLAB_COLS, run["slabs"], fm),
            "meta.json": (json.dumps(run["meta"], indent=1, sort_keys=True) + "\n").encode()}


def save_run(run, folder):
    os.makedirs(folder, exist_ok=True)
    for name, blob in serialise(run).items():
        with open(os.path.join(folder, name), "wb") as f:
            f.write(blob)


def load_run(folder):
    def table(name):
        t = np.genfromtxt(os.path.join(folder, name), delimiter=",", names=True)
        return {k: np.asarray(t[k], float) for k in t.dtype.names}
    return {"curve": table("curve.csv"), "slabs": table("slabs.csv"), "meta": json.load(open(os.path.join(folder, "meta.json")))}


# ------------------------------------------------------------------ sections and TQI


def _tqi(v, slab):
    v = v[np.isfinite(v)]
    if len(v) < 10:
        return np.nan
    return (np.sum(np.sqrt(slab ** 2 + np.diff(v) ** 2)) / (slab * (len(v) - 1)) - 1) * 1e6


def _midchord(v, half):
    out = np.full(len(v), np.nan)
    out[half:-half] = v[half:-half] - 0.5 * (v[:-2 * half] + v[2 * half:])
    return out


def section_table(run, p=P):
    sl = run["slabs"]
    sc = sl["s"]
    centre = np.where(np.isfinite(sl["d_left"]) & np.isfinite(sl["d_right"]), 0.5 * (sl["d_left"] + sl["d_right"]), np.nan)
    align = _midchord(centre, int(round(p.chord / p.slab / 2)))
    gdev = sl["gauge"] - NOMINAL_GAUGE
    rows = []
    for k in range(int(sc[-1] // p.section)):
        m = (sc >= k * p.section) & (sc < (k + 1) * p.section)
        if np.isfinite(gdev[m]).mean() < 0.5:
            continue
        g = gdev[m][np.isfinite(gdev[m])]
        rows.append({"section": k, "s_start": k * p.section, "s_end": (k + 1) * p.section, "gauge_dev_mean_mm": g.mean() * 1e3,
                     "gauge_dev_std_mm": g.std() * 1e3, "gauge_dev_min_mm": g.min() * 1e3, "gauge_dev_max_mm": g.max() * 1e3,
                     "align_std_mm": np.nanstd(align[m]) * 1e3, "tqi_gauge": _tqi(gdev[m], p.slab), "tqi_alignment": _tqi(align[m], p.slab)})
    return rows


def normalise(rows, norm_rows):
    ng = np.nanpercentile([r["tqi_gauge"] for r in norm_rows], 95)
    na = np.nanpercentile([r["tqi_alignment"] for r in norm_rows], 95)
    for r in rows:
        r["tqi_combined"] = 0.5 * (r["tqi_gauge"] / ng + r["tqi_alignment"] / na)
    return rows


# ------------------------------------------------------------------ the 3D curve itself


def sample_curve(curve, step=0.1, p=P):
    """The run's 3D curve (Catmull-Rom through the control points) sampled every `step` m: E, N, Z per rail and centreline."""
    grid = np.arange(np.ceil(curve["s"][0] / step) * step, curve["s"][-1], step)
    out = {"s": grid}
    for side in ("left", "right"):
        for axis, col in (("E", "E_" + side), ("N", "N_" + side), ("Z", "z_" + side)):
            out[f"{axis}_{side}"] = cr.evaluate(curve["s"], curve[col], grid, 3 * p.ctrl)
    for axis in "ENZ":
        out[f"{axis}_centre"] = 0.5 * (out[f"{axis}_left"] + out[f"{axis}_right"])
    return out


# ------------------------------------------------------------------ comparison of runs

QUANTITIES = ["d_left", "z_left", "d_right", "z_right", "gauge", "d_centre", "z_centre"]


def compare_curves(a, b, p=P, step=0.25):
    """Change from run a to run b on the common chainage: b - a, for every curve quantity."""
    lo, hi = max(a["s"][0], b["s"][0]), min(a["s"][-1], b["s"][-1])
    grid = np.arange(np.ceil(lo / step) * step, hi, step)
    diff = {"s": grid}
    for q in QUANTITIES:
        diff[q] = cr.evaluate(b["s"], b[q], grid, 3 * p.ctrl) - cr.evaluate(a["s"], a[q], grid, 3 * p.ctrl)
    return diff


def summarise(diff, floor=None):
    out = {}
    for q in QUANTITIES:
        v = diff[q][np.isfinite(diff[q])]
        if len(v) == 0:
            continue
        row = {"rms_mm": float(np.sqrt(np.mean(v ** 2)) * 1e3), "p95_abs_mm": float(np.percentile(np.abs(v), 95) * 1e3),
               "max_abs_mm": float(np.abs(v).max() * 1e3), "samples": int(len(v))}
        if floor and q in floor:
            row["fraction_above_noise_floor"] = float(np.mean(np.abs(v) * 1e3 > floor[q]["p95_abs_mm"]))
        out[q] = row
    return out


def same_setup(meta_a, meta_b):
    keys = ("slab", "ctrl", "d_max", "section")
    return meta_a["reference"] == meta_b["reference"] and all(meta_a["params"][k] == meta_b["params"][k] for k in keys)


# ------------------------------------------------------------------ command line


def _folder(run_id):
    return os.path.join(OUT, "runs", run_id)


def cmd_reference(a):
    pts = read_points(a.laz)
    ref = make_reference(pts, P, a.name)
    os.makedirs(os.path.join(OUT, "reference"), exist_ok=True)
    path = os.path.join(OUT, "reference", a.name + ".json")
    with open(path, "w") as f:
        json.dump(ref, f, sort_keys=True)
    r = Reference(ref)
    print(f"reference {ref['id']} written to {path}: {r.length:.1f} m, {len(r.s)} samples")


def cmd_build(a):
    ref = load_reference(a.reference)
    pts = read_points(a.laz)
    mask = None if not a.split_half else ((np.arange(len(pts.E)) // 1000) % 2 == (a.split_half == "b"))
    run = build_run(pts, ref, P, a.run_id, mask)
    save_run(run, _folder(a.run_id))
    print(json.dumps(run["meta"]["qa"], indent=1))
    print("written to", _folder(a.run_id))


def cmd_sections(a):
    run = load_run(_folder(a.run_id))
    rows = section_table(run)
    norm = section_table(load_run(_folder(a.norm_run or a.run_id)))
    rows = normalise(rows, norm)
    cols = list(rows[0])
    with open(os.path.join(_folder(a.run_id), "sections.csv"), "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(f"{r[c]:.3f}" for c in cols) + "\n")
    print("normalised by the 95th percentile of run", a.norm_run or a.run_id)
    print("sec  s_start  s_end | gauge dev mean  std   min   max | align std | TQI gauge  TQI align | combined")
    for r in rows:
        print(f"{r['section']:3d} {r['s_start']:8.0f} {r['s_end']:6.0f} | {r['gauge_dev_mean_mm']:13.2f} {r['gauge_dev_std_mm']:5.2f} "
              f"{r['gauge_dev_min_mm']:5.2f} {r['gauge_dev_max_mm']:5.2f} | {r['align_std_mm']:9.2f} | {r['tqi_gauge']:9.1f} "
              f"{r['tqi_alignment']:10.1f} | {r['tqi_combined']:7.2f}")


def cmd_sample(a):
    run = load_run(_folder(a.run_id))
    dense = sample_curve(run["curve"], a.step)
    keep = np.isfinite(dense["E_centre"]) & np.isfinite(dense["Z_centre"])
    dense = {k: v[keep] for k, v in dense.items()}
    cols = ["s"] + [f"{ax}_{k}" for k in ("left", "right", "centre") for ax in "ENZ"]
    path = os.path.join(_folder(a.run_id), "curve_dense.csv")
    with open(path, "wb") as f:
        f.write(_table_bytes(cols, dense, {c: "%.4f" for c in cols if c[0] in "EN"}))
    print(f"3D curve of {a.run_id} sampled every {a.step} m: {len(dense['s'])} points written to {path}")


def cmd_compare(a):
    ra, rb = load_run(_folder(a.run_a)), load_run(_folder(a.run_b))
    if not same_setup(ra["meta"], rb["meta"]):
        sys.exit("the runs use different reference frames or settings, so they cannot be compared")
    diff = compare_curves(ra["curve"], rb["curve"])
    floor_path = os.path.join(OUT, "noise_floor.json")
    floor = json.load(open(floor_path))["split_half"] if os.path.exists(floor_path) else None
    summary = summarise(diff, floor)
    os.makedirs(os.path.join(OUT, "comparisons"), exist_ok=True)
    name = f"{a.run_a}__{a.run_b}"
    cols = ["s"] + QUANTITIES
    with open(os.path.join(OUT, "comparisons", name + ".csv"), "wb") as f:
        f.write(_table_bytes(cols, diff, {}))
    print(f"change from {a.run_a} to {a.run_b} over {diff['s'][0]:.0f}-{diff['s'][-1]:.0f} m (b - a, mm):")
    for q, row in summary.items():
        extra = f" | above noise floor: {row['fraction_above_noise_floor']:.1%}" if "fraction_above_noise_floor" in row else ""
        print(f"  {q:9s} rms {row['rms_mm']:6.2f}  p95 {row['p95_abs_mm']:6.2f}  max {row['max_abs_mm']:6.2f}{extra}")
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
        for axi, q, label in zip(ax, ("d_centre", "z_centre", "gauge"), ("lateral offset", "height", "gauge")):
            axi.plot(diff["s"], diff[q] * 1e3, lw=.9)
            if floor and q in floor:
                axi.axhspan(-floor[q]["p95_abs_mm"], floor[q]["p95_abs_mm"], color="#D9EDEF", zorder=0)
            axi.set_ylabel(f"{label} (mm)")
        ax[0].set_title(f"change from {a.run_a} to {a.run_b} (shaded: 95 % noise floor)", fontsize=10)
        ax[-1].set_xlabel("chainage s (m)")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, "comparisons", name + ".png"), dpi=110)


def injection_test(pts, ref, block, run_a, floor):
    """Run A is one half of the scan. Run B is the other half with a known change added. Does the comparison find it?"""
    s, d = ref.project(pts.E, pts.N)
    nx, ny = ref.normal_at(s)
    lat_win, z_win, dlat, dz = (100.0, 120.0), (150.0, 170.0), 0.004, 0.005
    left = (d > 0) & (s >= lat_win[0]) & (s < lat_win[1])
    high = (s >= z_win[0]) & (s < z_win[1])
    E2, N2, Z2 = pts.E + np.where(left, dlat, 0.0) * nx, pts.N + np.where(left, dlat, 0.0) * ny, pts.Z + np.where(high, dz, 0.0)
    run_b = build_run(Points(E2, N2, Z2, pts.name, pts.sha256), ref, P, "injected", block == 1)
    diff = compare_curves(run_a["curve"], run_b["curve"])
    g = diff["s"]
    res = {}
    print("3. injection test: run A = one half of the scan, run B = the other half with a known change added:")
    for q, win, true in (("d_left", lat_win, dlat), ("gauge", lat_win, dlat), ("z_centre", z_win, dz)):
        inside = (g >= win[0] + 3) & (g <= win[1] - 3)
        outside = ((g < win[0] - 3) | (g > win[1] + 3)) & np.isfinite(diff[q])
        v_in, v_out = diff[q][inside], diff[q][outside]
        thr = floor[q]["p95_abs_mm"] / 1e3
        res[q] = {"injected_mm": true * 1e3, "recovered_mean_mm": round(float(np.nanmean(v_in) * 1e3), 3),
                  "detected_fraction": round(float(np.nanmean(np.abs(v_in) > thr)), 3),
                  "false_alarm_fraction": round(float(np.nanmean(np.abs(v_out) > thr)), 3)}
        print(f"   {q:9s} injected {true * 1e3:4.1f} mm, recovered {res[q]['recovered_mean_mm']:5.2f} mm, "
              f"detected on {res[q]['detected_fraction']:.0%} of the stretch, false alarms elsewhere {res[q]['false_alarm_fraction']:.0%}")
    return res


def cmd_selftest(a):
    ref = load_reference(a.reference)
    pts = read_points(a.laz)
    r1, r2 = build_run(pts, ref, P, "selftest"), build_run(pts, ref, P, "selftest")
    same = all(x == y for x, y in zip(serialise(r1).values(), serialise(r2).values()))
    print("1. same input and settings twice gives byte-identical output:", same)
    idx = np.arange(len(pts.E))
    block = (idx // 1000) % 2
    ha, hb = build_run(pts, ref, P, "half_a", block == 0), build_run(pts, ref, P, "half_b", block == 1)
    diff = compare_curves(ha["curve"], hb["curve"])
    floor = summarise(diff)
    print("2. two complementary halves of the scan (alternating blocks of 1000 points), full pipeline including projection:")
    for q, row in floor.items():
        print(f"   {q:9s} rms {row['rms_mm']:5.2f} mm   95 % of |change| below {row['p95_abs_mm']:5.2f} mm")
    inj = injection_test(pts, ref, block, ha, floor)
    with open(os.path.join(OUT, "noise_floor.json"), "w") as f:
        json.dump({"note": "change between two halves of the same scan: sampling noise of the pipeline, at half the point density; "
                           "positioning error between real runs comes on top",
                   "input": pts.name, "input_sha256": pts.sha256, "reference": ref.id, "code": CODE_HASH,
                   "split_half": {q: {"rms_mm": round(v["rms_mm"], 3), "p95_abs_mm": round(v["p95_abs_mm"], 3)} for q, v in floor.items()},
                   "injection_test": inj},
                  f, indent=1, sort_keys=True)
    if not same:
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("reference"); s.add_argument("laz"); s.add_argument("--name", default="track"); s.set_defaults(f=cmd_reference)
    s = sub.add_parser("build"); s.add_argument("laz"); s.add_argument("--run-id", required=True); s.add_argument("--reference", required=True)
    s.add_argument("--split-half", choices=["a", "b"], help="testing aid: use only half of the points (alternating blocks of 1000)"); s.set_defaults(f=cmd_build)
    s = sub.add_parser("sections"); s.add_argument("run_id"); s.add_argument("--norm-run"); s.set_defaults(f=cmd_sections)
    s = sub.add_parser("sample"); s.add_argument("run_id"); s.add_argument("--step", type=float, default=0.1); s.set_defaults(f=cmd_sample)
    s = sub.add_parser("compare"); s.add_argument("run_a"); s.add_argument("run_b"); s.add_argument("--plot", action="store_true"); s.set_defaults(f=cmd_compare)
    s = sub.add_parser("selftest"); s.add_argument("laz"); s.add_argument("--reference", required=True); s.set_defaults(f=cmd_selftest)
    a = ap.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
