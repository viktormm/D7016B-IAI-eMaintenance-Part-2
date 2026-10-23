# Track Gauge from a LiDAR Point Cloud, Checked Against Trafikverket's Maintenance Limits

**Course:** D7016B – Industrial AI and eMaintenance, Part II — Luleå University of Technology
**Author:** Viktor Lindberg
**Supervisor:** Amit Patwardhan

## What this project is about

One LiDAR scan of a railway track is turned into a 3D curve of both rails. From that
curve the **gauge**, the distance between the inner faces of the rails, is read every
25 cm and every metre along the track and checked against Trafikverket's maintenance
limits per speed class (TRVINFRA-00013).

A last step, in `report_figures.py`, lets the metres
of track sort themselves into **states** with k-means on one feature vector per metre.

Everything is on one interactive page, `output/viewer/monitor.html`: the gauge along
the chainage against the limits, a table of flagged stretches, and the curve in 3D.

The project was scaled down to this after discussion with the supervisor (October
2026). Earlier parts, a simulated sequence of runs with assumed defects, planted
faults, a point-density study and a multi-step dashboard, were removed from the
repository; they remain in the git history before that date.

**Why the gauge.** The gauge is measured inside one cross-section of the track,
relative to the rails. It does not depend on how well two scans are positioned against
each other, and it does not depend on the absolute height. The height in this scan is
relative to the sensor, which rides on the rails, so height-based track parameters
(cross-level, cant, twist, long-wave surface) cancel out and cannot be read from it.
(An earlier version also read the overhead contact wire from the scan; it was removed
on 15 October 2026 to keep the project to the gauge, and remains in the git history.)

## From points to curve

`track_curve.py` turns a `.laz` scan into a run that can be compared with other runs
of the same track. What makes it repeatable:

- **One fixed reference frame.** The smoothed track centreline of the first run is
  stored once (`output/reference/track.json`). Every run, including the first, is
  projected onto it: chainage `s` along the reference, sideways offset `d` from it
  (left is `+d`), height `Z` as stored.
- **Fixed grids anchored at `s = 0`.** Slabs (25 cm), control points (1 m) and
  sections (25 m) sit on the same grid in every run.
- **No random steps and no data-dependent settings.** All settings live in one
  `Params` block and are stored with every run, together with the input file's
  SHA-256, the reference id and a hash of the code. `compare` refuses runs made with
  a different reference or grid.
- **Checked by a self-test** (`selftest`): rebuilding gives byte-identical files, and
  two halves of the same scan give the noise floor of the whole pipeline.

The curve of a run is:

1. Project the points onto the reference; keep those within 2.5 m of it.
2. In every 25 cm slab, find the rail-head top of each rail (highest points in a window that follows the rail) and the inner (gauge) face 8–20 mm below it. The gauge of the slab is the distance between the two faces.
3. Reject outliers against a running median.
4. Average into control points every 1 m on the fixed grid (a line fitted through the slabs in each bin, evaluated at the bin centre).
5. Join the control points with a centripetal Catmull-Rom curve. It passes through every control point, and a change in one place only moves the curve nearby.

### Running it

From the repo root, in this order:

```
python track_curve.py reference gps_pointcloud_00.laz --name track       # once, from the first run
python track_curve.py build gps_pointcloud_00.laz --run-id run1 --reference output/reference/track.json
python track_curve.py sample run1                                        # the 3D curve every 0.1 m
python track_curve.py sections run1                                      # 25 m sections and TQI
python track_curve.py selftest gps_pointcloud_00.laz --reference output/reference/track.json
python curve_attributes.py run1                                          # optional: ballast, sleepers, head normals (used as features)
python track_features.py run1                                            # one feature vector per metre
python view_curve.py run1                                                # the 3D view
python monitor.py run1                                                   # the page, opens in the browser
```

To add a second run and compare it with the first (same reference, same settings):

```
python track_curve.py build run2.laz --run-id run2 --reference output/reference/track.json
python track_curve.py compare run1 run2 --plot
```

Output per run in `output/runs/<run-id>/`:

| File | Content |
|---|---|
| `slabs.csv` | The per-slab measurements: `s`, the face position `d_left`/`d_right`, the head height `z_left`/`z_right`, `gauge` (m), and how many points were used |
| `curve.csv` | The curve as control points every 1 m: `s`, `d_left`, `z_left`, `d_right`, `z_right`, `gauge`, `d_centre`, `z_centre`, and `E`/`N` of both rails |
| `curve_dense.csv` | The 3D curve (E, N, Z for both rails and the centreline) sampled every 0.1 m |
| `sections.csv` | 25 m sections with gauge statistics, alignment and TQI |
| `meta.json` | Settings, input hash, reference id, code hash and quality figures |
| `features.csv` | One row per metre: gauge, 10 m alignment versine, points per metre, rail head above ballast and above sleeper, head tilt, sleeper spacing, rail-head intensity left and right |
| `curve_attributes.csv`, `sleepers.csv`, `shortwave.csv` | Optional per-point attributes from `curve_attributes.py` |

## The monitor page

```
python monitor.py run1          # writes output/viewer/monitor.html and opens it
```

One self-contained page (works offline) with one tab, **Gauge**: the limit menu, the
chart along the chainage, the curve in 3D and the table of flagged stretches.
Everything is computed in the page, so the controls react at once.

**On the page:**

- **Reference limits** in a menu, each with one or more levels; a point is marked by the most severe level it exceeds. Trafikverket's requirements for track geometry, TRVINFRA-00013 version 2.0 (2021-10-11): table K11.2, the deviation of a single reading from 1435 mm per speed class H0–H5 with the levels PLAN (planning limit), UH1 and UH2 (maintenance limits; UH2 is a safety requirement) and KRIT (critical, immediate action), and table K11.3, minimums for the 100 m sliding mean (choose the 100 m sliding mean as the points of the chart). H4 (160–200 km/h) is preselected, because the scan lies on Södra stambanan south-west of Höör (inferred from its coordinates, see Results); it can be changed in the menu. Points outside the limits are red; the stretches are shaded red.
- **Stretches, not points.** A fault is a stretch, so the gauge is averaged over 1 m by default (2 m and 5 m available, as are the raw slabs and the 1 m control points). Flagged stretches shorter than a chosen length are ignored. (The Western Electric control-chart rules of an earlier version are still in the page's code but switched off and hidden: they judged a stretch against the track's own mean, which is not how track maintenance is decided, and they were dropped from the report.)
- **The curve in 3D.** The page embeds the 3D viewer (`view_curve.py`).

**Track states** are not on the page any more (removed in October 2026 to keep the
page to the limit check). The k-means grouping of the per-metre feature vectors
(`features.csv`) lives in `report_figures.py`, which draws the report's Figure 5 and
prints the state profiles. k-means knows nothing about tracks: a state is a combination
of feature levels that occur together along this track. It runs with fixed seeds, so
it gives the same result every time.

## The 3D view

```
python view_curve.py run1            # writes output/viewer/curve3d_run1.html and opens it
```

One self-contained page: drag to rotate, scroll to zoom. It shows the Catmull-Rom
curve of each rail, the 1 m control points, the centreline and a cross-link every
5 m; **Colour the rails by gauge** paints the curve by the gauge; the **Window** zooms
to a stretch and the **Cursor** reads out gauge and offsets along the track.
The viewer also writes `output/viewer/<run>_curve.las`, the curve
as points in the scan's own coordinates, for CloudCompare on top of
`gps_pointcloud_00.laz`.

## Results (one run: `gps_pointcloud_00.laz`)

- 15.9 million points over 264 m of double track, with a gentle curve at the start. The LiDAR followed the left track, so only that track is used. The only per-point attribute with values is the intensity; return number, classification and the rest are empty.
- **Rails and gauge:** both gauge faces are found in 1035 of the 1057 slabs (98 %). Mean gauge 1432.3 mm, 2.7 mm below nominal, scatter 1.6 mm; tightest slab 1424.7 mm, widest 1436.1 mm. The whole track is tight: 97 % of the slabs lie below 1435 mm. Whether that is the track or a millimetre of bias in the gauge-face detection cannot be told from one scan; the differences along the track are unaffected by a constant offset.
- **A step at 155 m (about km 566.29):** the 5 m mean gauge falls from 1434.9 mm at 152 m to 1430.6 mm at 160 m. Before it the track averages 1433–1434 mm, after it 1431–1432 mm. The rail head also sits higher above the sleepers after it (184.6 mm over 130–154 m, 186.4 mm over 156–180 m, 188.3 mm over 184–258 m). The sleeper spacing has a 0.8 m gap at 155–156 m, but also at six other places, so that is not specific to the step; the left rail-head intensity rises locally after the step but is lower again after 184 m. Something physical sits there (a joint, a change of rail or fastening); nothing is mapped in OpenStreetMap at that point. It is the first thing to ask about.
- **Where the scan is (inferred, to be confirmed):** the file has no coordinate reference system, but its coordinates fit SWEREF 99 TM and place the scan on Södra stambanan between Stehag and Höör, 2 km south-west of Höör station, at km 566.44–566.18; the scan ends at Höör's distant signals Hö 22/Hö 52, which OpenStreetMap puts at km 566.180. OSM tags the stretch 180 km/h (Trafikverket lists the line at 200 km/h), so the speed class is H4 (160–200 km/h). The 155 m step is at about km 566.29 and the tight stretch 186–216 m at km 566.23–566.26.
- **Against Trafikverket's maintenance limits (TRVINFRA-00013, table K11.2, point faults):** the wide side is never reached. On the tight side, for the line's class H4, 103 of 260 one-metre stretches are beyond PLAN (below 1432 mm), 46 beyond UH1 (1431 mm) and 17 beyond UH2 (1430 mm), none beyond KRIT (1428 mm); the tight stretch 186–216 m holds most of the UH2 metres. For class H2 (80–120 km/h) PLAN and UH1 coincide at 1431 mm, so the counts are 46/46/17/0; for H5 they are 103/103/46/17. On the raw 25 cm slabs for H4 the counts are 420, 214, 88 and 6 (the KRIT slabs are at the scan edge and in the tight stretch). **Table K11.3, 100 m sliding mean:** the mean runs from 1431.1 to 1433.5 mm, within the H1–H2 limits, below UH1min for H3 (1432 mm, 55 of 165 windows), H4 (1433 mm, 94 windows) and H5 (1434 mm, every window), never below UH2min. So for its own class the stretch is a maintenance case: 17 metres past UH2 to repair and the gauge to be planned for correction over most of its length. A constant bias of a millimetre in the gauge-face detection would move these counts (H4, ±1 mm: UH2 between 3 and 46, KRIT 0 to 3), so the first check is against a track-recording-car run from Trafikverket's Optram system.
- **Noise floor of the pipeline** (two complementary halves of the scan compared): gauge 0.9 mm RMS, 95 % of the differences below 1.7 mm per slab; slab-to-slab noise 0.68 mm. A gauge change smaller than that cannot be told from the method's own scatter; positioning error between real scans comes on top.
- **Systematic error, bounded without a reference measurement** (`python error_checks.py`, rebuilds the run four times): sampling the gauge face 10–18 mm below the rail top raises the mean gauge by 0.37 mm; 12–16 mm by 0.84 mm (but finds only 450 of 1057 slabs); 8–24 mm lowers it by 0.37 mm. The band alone moves the mean by 0.4–0.8 mm, so the pipeline's own offset against the standard's band (the smallest distance anywhere from the rail top down to 14 mm below it, which lies higher on the rail than any band the scan can resolve) is of that order; the scanner's calibration is unknown. A constant offset of ±1 mm moves the class H4 count of 1 m stretches beyond UH2 between 46 and 3 (beyond PLAN between 172 and 46, beyond UH1 between 103 and 17); KRIT (1428 mm) is reached by 3 metres only if the reading is 1 mm too wide, and by 17 at 2 mm. The differences along the track are unaffected by an offset.
- **Track states** (`report_figures.py`, k = 3, all 9 features, 254 metres with every feature): state 1 is 5–63 m, state 2 mainly 64–155 m (78–132 m without a break) and state 3 mainly 156–258 m (183–258 m without a break); between 133 and 182 m states 2 and 3 alternate. The boundary between states 2 and 3 falls at the gauge step at 155 m. Mean per state: gauge −2.6 / −1.5 / −3.9 mm from 1435, rail head above ballast 225 / 227 / 242 mm, above sleeper 182 / 183 / 188 mm, left rail-head intensity 1469 / 1662 / 1622, points per metre 48 400 / 33 100 / 32 400. Only the first measured metre (5 m, at the scan's edge) scores above 2σ from its own state's centre (2.1); every other metre is within 1.5σ.

## How the solution was developed

The solution was developed with Claude Code, Anthropic's AI coding assistant, used as
a tool under the author's clear instructions: the author designed the pipeline, set
the scope, the research questions, the limits and the checks, instructed the tool step
by step and reviewed and tested each result; the tool assisted by producing and
revising the code, the figures and the page under those instructions.
The reason was speed under a scope that changed several times. Because that makes
the code easy to treat as a black box, correctness rests on tests rather than trust:
the half-scan noise floor (`track_curve.py selftest`), the band and offset checks
(`error_checks.py`), the automatic test of the monitor page, and every number in the
report traced back to the files in `output/`.

## Assumptions and limits

- One scan only; nothing about change over time can be shown yet. The pipeline is built for a second run, which would be processed with the same reference and compared point for point.
- The gauge is judged against Trafikverket's maintenance limits (TRVINFRA-00013). (The FRA safety limits of 49 CFR 213.53, the reference behind Sharma et al.'s "red-tag" defects, were used in an earlier version and dropped in October 2026: no reading came near them, and they are not how Swedish maintenance is decided.) The line and its speed class (Södra stambanan, H4) are inferred from the file's coordinates and OpenStreetMap, not stated in the data, so the Trafikverket tables are applied per class and H4 is only preselected.
- Trafikverket defines the gauge as the smallest distance between the rails measured anywhere from the rail top down to 14 mm below it, unfiltered (TRVINFRA-00013, terms and section 8.2.1); the pipeline samples the gauge face between 8 and 20 mm below the top, a little lower than that.
- Height is relative to the sensor, so of the track parameters only gauge and lateral alignment are measured in full.
- The reference covers the first run's extent plus 2 m at each end.

## Background: gauge, TQI and the source paper

Sharma et al. (2018) work with five track geometry parameters (gauge, cross-level,
cant, twist, surface), define the Track Quality Index per parameter and section as
`TQI = (Ls / L0 − 1) × 10^6` (the measured trace length against the section length),
and treat geometry measurements that violate the FRA thresholds as "red-tag" defects
that must be repaired before the next inspection. This project uses their definition
of gauge; the limits it checks against are Trafikverket's. `track_curve.py sections` still
computes the TQI per 25 m section (`sections.csv`), with 0.25 m samples in metres, so
the values are not comparable with the paper's foot-by-foot ones.

## Repository contents

| File | Purpose |
|---|---|
| `track_curve.py` | The pipeline: reference frame, per-run curve and gauge, sections and TQI, comparison of runs, self-test |
| `catmull_rom.py` | Centripetal Catmull-Rom curve math used by the pipeline |
| `track_features.py` | One feature vector per metre from the pipeline's outputs and the rail-head intensity |
| `error_checks.py` | The two checks of the systematic error without a reference measurement: the gauge-face band and the offset sensitivity of the Trafikverket verdict |
| `report_figures.py` | The report's figures from the outputs (`output/figures/r1`–`r3`) and the k-means track states |
| `monitor.py`, `monitor_template.html` | The page: the gauge against Trafikverket's limits, flagged stretches and the curve in 3D |
| `view_curve.py`, `view_curve_template.html` | The 3D viewer and the CloudCompare export |
| `curve_attributes.py` | Optional per-point attributes of the curve (normals, ballast, sleepers, short-wave height indicators); some of them are features |
| `gps_pointcloud_00.laz` | The scan: LiDAR point cloud of a two-track line (LAS 1.2, point format 0, no timestamps or coordinate system) |
| `output/reference/track.json` | The reference frame every run is projected onto |
| `output/runs/run1/` | The curve and measurements of run 1 (see the table above) |
| `output/noise_floor.json` | Noise floor and known-change check from `selftest` |
| `output/viewer/` | `monitor.html`, the 3D view `curve3d_run1.html` and the curve as LAS files for CloudCompare |
| `output/figures/` | The report's figures: `f2_rail_section.png` and `f6_catmull_rom.png` (method, from the first scope) and `r1`–`r3` from `report_figures.py` |
| `D7016B_Demo1_Viktor_Lindberg.pptx`, `D7016B_Demo2_Viktor_Lindberg.pptx` | Slides for Demonstrations 1 and 2 (the Demo 2 deck predates the scale-down) |
| `D7016B_Final_Viktor_Lindberg.pptx` | The final presentation, 12 slides: how the project changed, the point cloud, cleaning, method, the method's own error, Trafikverket's two checks and the results, two places to visit, k-means states, the next step (data over time, eMaintenance, predictive maintenance, Industrial AI); the speaker notes hold the manus |
| `D7016B_Final_manus_Viktor_Lindberg.docx` | The manus (speaker script) of the final presentation, one section per slide, same text as the speaker notes |
| `Industrial AI - Part II .docx` | Final report, based on the course report template. Local drafts of it (`Report rough draft.docx`) are ignored by git, see `.gitignore` |
| `PointCloud_Project_Time_Plan_Viktor - Rough_Final (2).xlsx` | The original time plan |
| `1-s2.0-S0968090X18302444-main.pdf` | Sharma et al. (2018), the source for the geometry parameters, the TQI and the maintenance models of the discussion |

## Tools

Python (needs `numpy`, `scipy`, `plotly` and `laspy[lazrs]`; the `lazrs` backend is
required to read `.laz` files, and `plotly` is used to embed plotly.js in the pages).
CloudCompare (`winget install CloudCompare.CloudCompare`) is used to inspect the point
cloud visually.

## References

Sharma, S., Cui, Y., He, Q., Mohammadi, R. and Li, Z. (2018). "Data-driven
optimization of railway maintenance for track geometry". Transportation Research
Part C, 90, pp. 34–58. https://doi.org/10.1016/j.trc.2018.02.019

Trafikverket (2021). Krav TRVINFRA-00013, Banöverbyggnad – Spårläge, version 2.0,
2021-10-11. Trafikverkets infrastrukturregelverk. Tables K7.1 (speed classes), K11.2
(gauge point faults) and K11.3 (gauge, 100 m sliding mean).
