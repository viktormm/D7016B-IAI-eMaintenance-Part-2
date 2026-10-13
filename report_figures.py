"""Figures for the report, from the outputs of run1: output/figures/r1..r4.png (column width, 300 dpi).

  python report_figures.py [run1]

  r1_gauge_limits.png   gauge per 1 m stretch along the track against Trafikverket's point-fault levels (class H4, the line's class)
  r2_gauge_100m.png     the 100 m sliding mean against the UH1 minimum of each speed class
  r3_contact_wire.png   contact wire height above the rail and stagger, against the TSI values
  r4_states.png         the track states from k-means on the per-metre features, with three of the features
Also prints the numbers the report quotes.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import monitor as mon
import track_curve as tc

RUN = sys.argv[1] if len(sys.argv) > 1 else "run1"
FOLDER = os.path.join(tc.OUT, "runs", RUN)
OUT = os.path.join(tc.OUT, "figures")
TEAL, ORANGE, PURPLE, GREY, RED, DARK = "#0f8b8d", "#e07a1f", "#6f5aa8", "#9fb0b5", "#c62828", "#172126"
STATE_COLOURS = [TEAL, ORANGE, PURPLE, "#8a7a66", "#2b6cb0"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#eef2f3", "grid.linewidth": 0.6})
W = 3.3   # column width in inches


def columns(path):
    t = np.genfromtxt(path, delimiter=",", names=True)
    return {k: np.asarray(t[k], float) for k in t.dtype.names}


def block_means(s, v, L):
    centres = np.arange(np.floor(s[0] / L) * L + L / 2, s[-1], L)
    out = np.full(len(centres), np.nan)
    for i, c in enumerate(centres):
        m = (s >= c - L / 2) & (s < c + L / 2) & np.isfinite(v)
        if m.sum() >= 0.5 * L / 0.25:
            out[i] = v[m].mean()
    return centres, out


def sliding_mean(s, v, L=100.0, step=1.0):
    centres = np.arange(np.ceil((s[0] + L / 2) / step) * step, s[-1] - L / 2, step)
    out = np.full(len(centres), np.nan)
    for i, c in enumerate(centres):
        m = (s >= c - L / 2) & (s < c + L / 2) & np.isfinite(v)
        if m.sum() >= 0.8 * L / 0.25:
            out[i] = v[m].mean()
    return centres, out


def stretches(mask, s, step, min_len=1.0):
    out, open_, last = [], None, -1e9
    gap = max(1, round(1.0 / step))
    for i, on in enumerate(mask):
        if not on:
            continue
        if open_ and i - last <= gap:
            open_[1] = i
        else:
            if open_:
                out.append(open_)
            open_ = [i, i]
        last = i
    if open_:
        out.append(open_)
    return [(s[a], s[b]) for a, b in out if (b - a + 1) * step >= min_len]


run = tc.load_run(FOLDER)
sl = run["slabs"]
s, g = sl["s"], sl["gauge"] * 1e3
nom = tc.NOMINAL_GAUGE * 1e3
s1, g1 = block_means(s, g, 1.0)
os.makedirs(OUT, exist_ok=True)

# ---- r1: gauge per 1 m against the H4 point-fault levels (the line's class: Södra stambanan, 180 km/h)
lv = mon.TRV_POINT["H4"]
fig, ax = plt.subplots(figsize=(W, 2.3))
ax.axhline(nom, color=GREY, lw=0.8, ls="--")
for name, col, ls in (("PLAN", "#f2b35c", ":"), ("UH1", "#e07a1f", ":"), ("UH2", RED, "--"), ("KRIT", "#7b1113", "-")):
    ax.axhline(nom + lv[name][0], color=col, lw=0.9, ls=ls, label=f"{name} {nom + lv[name][0]:.0f} mm")
ok = np.isfinite(g1)
ax.plot(s1[ok], g1[ok], color=TEAL, lw=0.9, label="gauge, 1 m mean")
beyond = ok & (g1 < nom + lv["UH2"][0])
ax.plot(s1[beyond], g1[beyond], ".", color=RED, ms=4, label="beyond UH2")
for a, b in stretches(ok & (g1 < nom + lv["UH1"][0]), s1, 1.0, 1.0):
    ax.axvspan(a - 0.5, b + 0.5, color="#c62828", alpha=0.08, lw=0)
ax.set_xlabel("chainage along the track (m)"); ax.set_ylabel("gauge (mm)"); ax.set_ylim(1423.0, 1437)
ax.legend(loc="lower left", ncol=2, frameon=True, framealpha=0.9, edgecolor="none", fontsize=6.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "r1_gauge_limits.png"), dpi=300); plt.close(fig)

# ---- r2: 100 m sliding mean against UH1min per class
sc, gm = sliding_mean(s, g)
fig, ax = plt.subplots(figsize=(W, 2.0))
ax.plot(sc, gm, color=TEAL, lw=1.2, label="gauge, 100 m sliding mean")
SPEED = {"H2": "80–120 km/h", "H3": "120–160 km/h", "H4": "160–200 km/h", "H5": "200–250 km/h"}
for cls, col in (("H2", "#f2b35c"), ("H3", "#e07a1f"), ("H4", RED), ("H5", "#7b1113")):
    u1 = mon.TRV_MEAN100[cls][0]
    ax.hlines(u1, sc[0], sc[-1], color=col, lw=0.9, ls="--")
    ax.text(sc[-1] + 3, u1, f"UH1 min {cls}, {SPEED[cls]}", color=col, fontsize=6.0, va="center")
ax.set_xlim(sc[0] - 5, sc[-1] + 85); ax.set_ylim(1429.5, 1435.2)
ax.set_xlabel("chainage of the window centre (m)"); ax.set_ylabel("gauge (mm)")
ax.legend(loc="upper left", frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "r2_gauge_100m.png"), dpi=300); plt.close(fig)
share = {cls: float(np.mean(gm[np.isfinite(gm)] < mon.TRV_MEAN100[cls][0])) for cls in mon.TRV_MEAN100}
print("100 m mean: %.1f to %.1f mm; share of windows below UH1min:" % (np.nanmin(gm), np.nanmax(gm)), {k: f"{100 * v:.0f} %" for k, v in share.items()})
for cls in ("H2", "H3", "H4", "H5"):
    L = mon.TRV_POINT[cls]
    print(f"point faults {cls}: 1 m stretches beyond PLAN {int((g1 < nom + L['PLAN'][0]).sum())}, UH1 {int((g1 < nom + L['UH1'][0]).sum())}, UH2 {int((g1 < nom + L['UH2'][0]).sum())}, KRIT {int((g1 < nom + L['KRIT'][0]).sum())} of {int(ok.sum())}")

# ---- r3: contact wire
cw = columns(os.path.join(FOLDER, "contact_wire.csv"))
fig, axes = plt.subplots(2, 1, figsize=(W, 3.2), sharex=True)
okw = np.isfinite(cw["height_m"])
dbl = okw & (cw["double"] > 0)
ax = axes[0]
ax.axhspan(5.00, 5.75, color=TEAL, alpha=0.07, lw=0); ax.axhline(5.75, color=RED, lw=0.9, ls="--"); ax.axhline(5.00, color=RED, lw=0.9, ls="--")
ax.plot(cw["s"][okw], cw["height_m"][okw], color=DARK, lw=1.0, label="contact wire height")
ax.plot(cw["s"][dbl], cw["height_m"][dbl], "x", color=GREY, ms=3.5, mew=0.8, label="two wires side by side")
ax.text(2, 5.76, "TSI nominal range 5.00–5.75 m", color=RED, fontsize=6.5, va="bottom")
ax.set_ylabel("height above rail (m)"); ax.set_ylim(5.4, 6.15); ax.legend(loc="upper left", frameon=False)
ax = axes[1]
ax.axhspan(-0.4, 0.4, color=TEAL, alpha=0.07, lw=0)
for y in (0.4, -0.4):
    ax.axhline(y, color=RED, lw=0.9, ls="--")
for y in (0.55, -0.55):
    ax.axhline(y, color="#7b1113", lw=0.9, ls=":")
ax.axhline(0, color=GREY, lw=0.8)
ax.plot(cw["s"][okw], cw["stagger_m"][okw], color=DARK, lw=1.0, label="stagger (+ left)")
ax.plot(cw["s"][dbl], cw["stagger_m"][dbl], "x", color=GREY, ms=3.5, mew=0.8)
ax.text(2, 0.41, "TSI ±0.40 m (1600 mm head), ±0.55 m (1950 mm head)", color=RED, fontsize=6.5, va="bottom")
ax.set_ylabel("stagger (m)"); ax.set_xlabel("chainage along the track (m)"); ax.set_ylim(-0.75, 0.8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "r3_contact_wire.png"), dpi=300); plt.close(fig)
print(f"wire: {int(okw.sum())} of {len(okw)} m, height {np.nanmin(cw['height_m']):.2f}-{np.nanmax(cw['height_m']):.2f} m (median {np.nanmedian(cw['height_m']):.2f}), "
      f"stagger {np.nanmin(cw['stagger_m']):+.2f}..{np.nanmax(cw['stagger_m']):+.2f}, above 5.75: {int((cw['height_m'] > 5.75).sum())} m, |stagger|>0.40: {int((np.abs(cw['stagger_m']) > 0.40).sum())} m, >0.55: {int((np.abs(cw['stagger_m']) > 0.55).sum())} m, double: {int(dbl.sum())} m")

# ---- r4: track states (k-means, k = 3, all features, standardised, k-means++ with fixed seeds)
ft = columns(os.path.join(FOLDER, "features.csv"))
keys = [k for k in ft if k != "s"]
X = np.column_stack([ft[k] for k in keys])
rows = np.where(np.isfinite(X).all(axis=1))[0]
Z = (X[rows] - X[rows].mean(0)) / X[rows].std(0)


def kmeans(Zm, k, seed):
    rng = np.random.default_rng(seed)
    centres = [Zm[rng.integers(len(Zm))]]
    while len(centres) < k:
        d = np.min([((Zm - c) ** 2).sum(1) for c in centres], axis=0)
        centres.append(Zm[rng.choice(len(Zm), p=d / d.sum())])
    centres = np.array(centres)
    for _ in range(100):
        lab = np.argmin(((Zm[:, None, :] - centres[None]) ** 2).sum(2), axis=1)
        new = np.array([Zm[lab == j].mean(0) if (lab == j).any() else centres[j] for j in range(k)])
        if np.allclose(new, centres):
            break
        centres = new
    inertia = ((Zm - centres[lab]) ** 2).sum()
    return lab, centres, inertia


best = min((kmeans(Z, 3, seed) for seed in range(1, 13)), key=lambda r: r[2])
lab, centres = best[0], best[1]
order = []
for l in lab:
    if l not in order:
        order.append(l)
lab = np.array([order.index(l) for l in lab])
state = np.full(len(ft["s"]), -1); state[rows] = lab
LAB = {"gauge_mm": "gauge\n(mm from 1435)", "intensity_left": "head intensity\n(left rail)", "wire_height_m": "wire height\n(m)", "above_ground_mm": "head above\nballast (mm)"}
show = ["gauge_mm", "above_ground_mm", "intensity_left", "wire_height_m"]
fig, axes = plt.subplots(len(show) + 1, 1, figsize=(W, 3.6), sharex=True, gridspec_kw={"height_ratios": [0.35] + [1] * len(show)})
ax = axes[0]
for j in range(3):
    m = state == j
    ax.bar(ft["s"][m], np.ones(m.sum()), width=1.0, color=STATE_COLOURS[j], label=f"state {j + 1}")
ax.set_yticks([]); ax.set_ylabel("state", rotation=0, ha="right", va="center"); ax.legend(loc="upper left", bbox_to_anchor=(0, 2.4), ncol=3, frameon=False); ax.grid(False)
for ax, key in zip(axes[1:], show):
    for j in range(3):
        m = state == j
        ax.plot(ft["s"][m], ft[key][m], ".", color=STATE_COLOURS[j], ms=2.5)
    ax.set_ylabel(LAB[key], fontsize=6.5)
axes[-1].set_xlabel("chainage along the track (m)")
fig.tight_layout(h_pad=0.4); fig.savefig(os.path.join(OUT, "r4_states.png"), dpi=300); plt.close(fig)


def ranges(j):
    ss = ft["s"][state == j]
    out, a, b = [], None, None
    for v in ss:
        if a is None:
            a = b = v
        elif v - b <= 1.01:
            b = v
        else:
            out.append((a, b)); a = b = v
    if a is not None:
        out.append((a, b))
    return [(round(a), round(b)) for a, b in out if b - a >= 2]


print("states (k=3):", {f"state {j + 1}": {"metres": int((state == j).sum()), "stretches": ranges(j)} for j in range(3)})
for key in ("gauge_mm", "intensity_left", "wire_height_m", "above_ground_mm", "sleeper_spacing_m"):
    print("  ", key, "mean per state:", [round(float(np.nanmean(ft[key][state == j])), 2) for j in range(3)])
print("written:", ", ".join(f"r{i}" for i in range(1, 5)), "to", OUT)
