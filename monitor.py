"""Track gauge against reference limits:
output/viewer/monitor.html

  python monitor.py [run1] [--no-open]

One self-contained page for one run (the real scan by default):
  gauge                 every 25 cm slab and every 1 m control point, against Trafikverket's maintenance limits per
                        speed class (TRVINFRA-00013, tables K11.2 and K11.3)
The page has the series along the track against the chosen limits, the curve in 3D
and a table of the flagged stretches. All limit computation happens
in the page, so the controls react at once; the data are embedded and the page works offline.
"""
import argparse
import datetime
import json
import os
import re
import webbrowser

import numpy as np
from plotly.offline import get_plotlyjs

import track_curve as tc

NOMINAL_MM = tc.NOMINAL_GAUGE * 1e3
INCH = 25.4
# Trafikverket TRVINFRA-00013 v2.0 (2021-10-11), Spårläge. Speed classes of table K7.1; gauge point-fault limits of table K11.2 as
# (min, max) deviation from 1435 mm per level; 100 m sliding-mean minimums of table K11.3 (UH1min, UH2min; straight track and R > 2500 m).
TRV_CLASSES = [("H0", "V ≤ 40 km/h"), ("H1", "40 < V ≤ 80 km/h"), ("H2", "80 < V ≤ 120 km/h"), ("H3", "120 < V ≤ 160 km/h"), ("H4", "160 < V ≤ 200 km/h"), ("H5", "200 < V ≤ 250 km/h")]
TRV_POINT = {"H0": {"PLAN": (-4, 15), "UH1": (-4, 22), "UH2": (-5, 30), "KRIT": (-10, 35)},
             "H1": {"PLAN": (-4, 15), "UH1": (-4, 22), "UH2": (-5, 30), "KRIT": (-10, 35)},
             "H2": {"PLAN": (-4, 12), "UH1": (-4, 17), "UH2": (-5, 25), "KRIT": (-9, 33)},
             "H3": {"PLAN": (-3, 12), "UH1": (-4, 15), "UH2": (-5, 22), "KRIT": (-8, 33)},
             "H4": {"PLAN": (-3, 10), "UH1": (-4, 12), "UH2": (-5, 20), "KRIT": (-7, 28)},
             "H5": {"PLAN": (-3, 8), "UH1": (-3, 10), "UH2": (-4, 15), "KRIT": (-5, 28)}}
TRV_MEAN100 = {"H1": (1431, 1430), "H2": (1431, 1430), "H3": (1432, 1430), "H4": (1433, 1430), "H5": (1434, 1432)}
TRV_LEVELS = {"PLAN": "planning limit (alert limit): plan the work", "UH1": "maintenance limit 1: repair soon, before UH2 is reached",
              "UH2": "maintenance limit 2, a safety requirement: repair before the next inspection", "KRIT": "critical limit: act at once"}


def _list(a, dec):
    return [None if not np.isfinite(v) else round(float(v), dec) for v in np.asarray(a, float)]


def columns(path):
    if not os.path.exists(path):
        return None
    t = np.genfromtxt(path, delimiter=",", names=True)
    return {k: np.asarray(t[k], float) for k in t.dtype.names}


def collect(run_id):
    folder = os.path.join(tc.OUT, "runs", run_id)
    run = tc.load_run(folder)
    sl, cv, meta = run["slabs"], run["curve"], run["meta"]
    view = os.path.join(tc.OUT, "viewer", f"curve3d_{run_id}.html")
    D = {"generated": datetime.date.today().isoformat(),
         "run": {"id": run_id, "file": meta["input"]["file"], "points": meta["input"]["points"], "reference": meta["reference"], "code": meta["code"],
                 "length_m": round(float(sl["s"][-1] + tc.P.slab / 2), 1), "slabs": int(len(sl["s"])), "measured": int(np.isfinite(sl["gauge"]).sum())},
         "views": {"run": os.path.basename(view) if os.path.exists(view) else None}, "params": []}
    g = sl["gauge"] * 1e3
    D["params"].append({
        "key": "gauge", "label": "Gauge", "short": "Gauge", "unit": "mm", "decimals": 1, "nominal": round(NOMINAL_MM, 1),
        "what": "the distance between the inner faces of the two rails, nominally 1435 mm, read from the LiDAR points in every 25 cm slab",
        "base": {"name": "25 cm slab", "step": tc.P.slab, "s": _list(sl["s"], 3), "v": _list(g, 2)},
        "ctrl": {"name": "1 m control point of the curve", "step": tc.P.ctrl, "s": _list(cv["s"], 3), "v": _list(cv["gauge"] * 1e3, 2)},
        "limits": [{"group": "Trafikverket TRVINFRA-00013, point faults (table K11.2)", "name": f"Trafikverket {c}, {sp}, point faults", "series": None,
                    "note": "deviation of a single reading from 1435 mm, measured unfiltered 14 mm below the rail top; four levels, PLAN, UH1, UH2 and KRIT",
                    "levels": [{"name": lv, "min": round(NOMINAL_MM + lo, 1), "max": round(NOMINAL_MM + hi, 1), "meaning": TRV_LEVELS[lv]} for lv, (lo, hi) in TRV_POINT[c].items()]}
                   for c, sp in TRV_CLASSES] +
                  [{"group": "Trafikverket TRVINFRA-00013, 100 m sliding mean (table K11.3)", "name": f"Trafikverket {c}, {sp}, 100 m mean", "series": "m100",
                    "note": "the sliding mean of the gauge over 100 m, for straight track and curves with R > 2500 m; only minimums are set; choose the 100 m sliding mean as the points of the chart",
                    "levels": [{"name": "UH1min", "min": float(u1), "max": None, "meaning": TRV_LEVELS["UH1"]}, {"name": "UH2min", "min": float(u2), "max": None, "meaning": TRV_LEVELS["UH2"]}]}
                   for c, sp in TRV_CLASSES if c in TRV_MEAN100 for u1, u2 in [TRV_MEAN100[c]]],
        "default_limit": 4,
        "limits_text": "Trafikverket's requirements for track geometry, TRVINFRA-00013 version 2.0 (2021-10-11), tables K11.2 and K11.3: the deviation of a single gauge reading from 1435 mm, and the 100 m sliding mean, per speed class H0 to H5 and in four levels. PLAN is the planning limit, UH1 and UH2 the two maintenance limits (UH2 is a safety requirement: repair before the next inspection), and KRIT the critical limit where action is immediate. The scan lies on Södra stambanan south-west of Höör (km 566.2 to 566.4, 180 km/h), so class H4 (160 to 200 km/h) is preselected; it can be changed in the menu."})
    return D


def render(D):
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "monitor_template.html"), encoding="utf-8") as f:
        page = f.read()
    # the page holds two copies of the quantity panel (a tab shows one or two quantities); the second copy gets "_2" on every id
    m = re.search(r"<!--PANEL-->(.*?)<!--/PANEL-->", page, re.S)
    panel2 = re.sub(r'\b(id|for)="([A-Za-z0-9]+)"', r'\1="\2_2"', m.group(1))
    page = page[:m.end()] + "\n    " + panel2.replace('id="panel_2"', 'id="panel_2" style="display:none"') + page[m.end():]
    page = page.replace("/*DATA*/null", json.dumps(D, separators=(",", ":"), allow_nan=False), 1)
    return page.replace("/*PLOTLY*/", get_plotlyjs(), 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_id", nargs="?", default="run1")
    ap.add_argument("--no-open", action="store_true")
    a = ap.parse_args()
    D = collect(a.run_id)
    path = os.path.join(tc.OUT, "viewer", "monitor.html")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(render(D))
    print(f"monitor written to {path} ({os.path.getsize(path) / 1e6:.1f} MB): " + ", ".join(p["key"] for p in D["params"]))
    if not a.no_open:
        webbrowser.open("file:///" + os.path.abspath(path).replace(os.sep, "/"))


if __name__ == "__main__":
    main()
