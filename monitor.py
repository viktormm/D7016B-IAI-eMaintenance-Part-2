"""Track and overhead line against reference limits:
output/viewer/monitor.html

  python monitor.py [run1] [--no-open]

One self-contained page for one run (the real scan by default), with a Gauge tab and a Contact wire tab (height and stagger stacked):
  gauge                 every 25 cm slab and every 1 m control point, against Trafikverket's maintenance limits per
                        speed class (TRVINFRA-00013, tables K11.2 and K11.3)
  contact wire height   the wire above the rail-head plane every metre (contact_wire.py), against the TSI Energy values
                        (Regulation (EU) 1301/2014, 4.2.9.1) and Trafikverket's nominal height
  contact wire stagger  its sideways position relative to the track centre, against the TSI maximum lateral deviation
                        per pantograph length (4.2.9.2)
Every tab has the series along the track against the chosen limits, the curve in 3D
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
    cw = columns(os.path.join(folder, "contact_wire.csv"))
    if cw:
        D["params"].append({
            "key": "wire_height", "label": "Contact wire height above the rail", "short": "Wire height", "unit": "m", "decimals": 3, "nominal": 5.5,
            "what": "the height of the overhead contact wire above the plane of the two rail heads, measured every metre from the points of the lowest wire above the track (contact_wire.py); the pantograph must find the wire within a height window",
            "base": {"name": "1 m", "step": 1.0, "s": _list(cw["s"], 2), "v": _list(cw["height_m"], 4)},
            "flags": {"name": "two wires side by side (overlap) or a fitting", "s": _list(cw["s"][cw["double"] > 0], 2)},
            "limits": [{"group": "TSI Energy, Regulation (EU) 1301/2014", "name": "TSI Energy, nominal range", "note": "4.2.9.1: nominal height 5.00 to 5.75 m", "levels": [{"name": "nominal range", "min": 5.00, "max": 5.75, "meaning": "outside the nominal range of the TSI"}]},
                       {"group": "TSI Energy, Regulation (EU) 1301/2014", "name": "TSI Energy, maximum", "note": "4.2.9.1: at most 6.50 m with tolerances and uplift (6.20 m design); the minimum follows EN 50119 and is not set here", "levels": [{"name": "maximum", "min": None, "max": 6.50, "meaning": "above the TSI maximum"}]},
                       {"group": "Trafikverket", "name": "Trafikverket, nominal 5.50 m", "note": "Swedish nominal contact wire height 5.50 m; 5.65 m as the maximum of TRVINFRA-00153 is quoted from a search result and must be confirmed", "levels": [{"name": "maximum", "min": None, "max": 5.65, "meaning": "above the Swedish maximum (to be confirmed)"}]}],
            "default_limit": 0,
            "limits_text": "the height values of the European TSI for the energy subsystem (Regulation (EU) 1301/2014, point 4.2.9.1): a nominal height between 5.00 and 5.75 m and at most 6.50 m with tolerances and pantograph uplift. The Swedish nominal height is 5.50 m; the Swedish maximum is to be confirmed with the supervisor."})
        D["params"].append({
            "key": "wire_stagger", "label": "Contact wire stagger", "short": "Wire stagger", "unit": "m", "decimals": 3, "nominal": 0.0,
            "what": "the sideways position of the contact wire relative to the track centreline (positive to the left): the zigzag that spreads the wear over the pantograph strip; too far out and the wire leaves the strip",
            "base": {"name": "1 m", "step": 1.0, "s": _list(cw["s"], 2), "v": _list(cw["stagger_m"], 4)},
            "flags": {"name": "two wires side by side (overlap) or a fitting", "s": _list(cw["s"][cw["double"] > 0], 2)},
            "limits": [{"group": "TSI Energy, Regulation (EU) 1301/2014", "name": "TSI Energy, 1600 mm pantograph", "note": "table 4.2.9.2: at most 0.40 m from the track centre under cross wind", "levels": [{"name": "maximum", "min": -0.40, "max": 0.40, "meaning": "further from the centre than a 1600 mm head allows"}]},
                       {"group": "TSI Energy, Regulation (EU) 1301/2014", "name": "TSI Energy, 1950 mm pantograph", "note": "table 4.2.9.2: at most 0.55 m from the track centre under cross wind", "levels": [{"name": "maximum", "min": -0.55, "max": 0.55, "meaning": "further from the centre than a 1950 mm head allows"}]}],
            "default_limit": 0,
            "limits_text": "the maximum lateral deviation of the contact wire from the track centre in the European TSI for the energy subsystem (Regulation (EU) 1301/2014, table 4.2.9.2): 0.40 m for a 1600 mm pantograph head and 0.55 m for a 1950 mm head, under cross wind. The scan shows the wire at rest, so the wind allowance is on top. The value for Swedish 1800 mm heads is to be confirmed."})
        # a small 3D picture: the rails and the wire in the track frame, every metre
        okw = np.isfinite(cw["height_m"])
        zrail = 0.5 * (cv["z_left"] + cv["z_right"])
        D["wire3d"] = {"s": _list(cv["s"], 2), "dL": _list(cv["d_left"], 3), "dR": _list(cv["d_right"], 3), "zL": _list(cv["z_left"] - zrail, 3), "zR": _list(cv["z_right"] - zrail, 3),
                       "ws": _list(cw["s"][okw], 2), "wd": _list(cw["d_wire_m"][okw], 3), "wh": _list(cw["height_m"][okw], 3)}
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
