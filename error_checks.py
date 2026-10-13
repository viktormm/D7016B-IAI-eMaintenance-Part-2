"""Two checks without ground truth: (1) how the Trafikverket verdict moves with a constant gauge offset of -2..+2 mm;
(2) how much the absolute gauge moves when the gauge face is sampled in different height bands below the rail top."""
import dataclasses, os, sys
import numpy as np

import track_curve as tc
import monitor as mon

run = tc.load_run(os.path.join(tc.OUT, "runs", "run1"))
sl = run["slabs"]; s, g = sl["s"], sl["gauge"] * 1e3
nom = 1435.0
cent = np.arange(0.5, 264, 1.0)
m1 = np.array([np.mean(g[(s >= c - 0.5) & (s < c + 0.5) & np.isfinite(g)]) if ((s >= c - 0.5) & (s < c + 0.5) & np.isfinite(g)).sum() >= 2 else np.nan for c in cent])
ok = np.isfinite(m1)

print("1) verdict against Trafikverket H4 (the line's class) / H2 point-fault levels, if every reading had a constant offset:")
print("   offset   H4: beyond PLAN  UH1  UH2  KRIT | H2: PLAN  UH1  UH2  KRIT | 100 m mean min..max")
for off in (-2, -1, 0, 1, 2):
    v = m1[ok] + off
    sm = []
    for c in np.arange(50, 215, 1.0):
        k = ok_ = (s >= c - 50) & (s < c + 50) & np.isfinite(g); sm.append(np.mean(g[k]) + off)
    L2, L3 = mon.TRV_POINT["H4"], mon.TRV_POINT["H2"]
    print("   %+d mm    %3d %3d %3d %3d           %3d %3d %3d %3d        %.1f..%.1f" % (off,
          (v < nom + L2["PLAN"][0]).sum(), (v < nom + L2["UH1"][0]).sum(), (v < nom + L2["UH2"][0]).sum(), (v < nom + L2["KRIT"][0]).sum(),
          (v < nom + L3["PLAN"][0]).sum(), (v < nom + L3["UH1"][0]).sum(), (v < nom + L3["UH2"][0]).sum(), (v < nom + L3["KRIT"][0]).sum(), min(sm), max(sm)))

print("\n2) the gauge face sampled in different bands below the rail top (pipeline default 8-20 mm; the standard says 14 mm):")
pts = tc.read_points("gps_pointcloud_00.laz")
ref = tc.load_reference(os.path.join(tc.OUT, "reference", "track.json"))
base = None
for band in ((0.020, 0.008), (0.016, 0.012), (0.018, 0.010), (0.024, 0.008)):
    p = dataclasses.replace(tc.P, face_below=band)
    r = tc.build_run(pts, ref, p, "bandtest")
    gg = r["slabs"]["gauge"] * 1e3
    mean = np.nanmean(gg)
    if base is None:
        base = gg
    diff = gg - base
    print("   band %2.0f-%2.0f mm below the top: mean gauge %.2f mm, slabs found %d, difference to the default: mean %+.2f mm, 95%% of slabs within %.2f mm" % (
        band[1] * 1e3, band[0] * 1e3, mean, np.isfinite(gg).sum(), np.nanmean(diff), np.nanpercentile(np.abs(diff), 95)))
