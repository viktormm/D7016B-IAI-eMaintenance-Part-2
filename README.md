# Track Gauge and Overhead Contact Wire from a LiDAR Point Cloud, Checked Against Reference Limits

**Course:** D7016B – Industrial AI and eMaintenance, Part II — Luleå University of Technology
**Author:** Viktor Lindberg
**Supervisor:** Amit Patwardhan

## What this project is about

One LiDAR scan of a railway track is turned into a 3D curve of both rails. From that
curve, and from the same points, three things are read every metre along the track and
checked the same way:

1. the **gauge**, the distance between the inner faces of the rails, against the gauge
   limits of a track safety standard (the FRA limits per class of track, 49 CFR 213.53,
   which lie behind the "red-tag" geometry defects in Sharma et al., 2018);
2. the **height of the overhead contact wire** above the rails and
3. its **stagger**, the zigzag of the wire over the track centre, both against the values
   of the European TSI for the energy subsystem (Regulation (EU) 1301/2014).

Each is watched with **Western Electric control-chart rules**, which mark the stretches
where the value steps out of the track's own normal scatter. A last step lets the metres
of track sort themselves into **states** with k-means on one feature vector per metre,
with an anomaly score for the metres that fit no state well.

Everything is on one interactive page, `output/viewer/monitor.html`, with a tab per
quantity: the value along the chainage, the limits, the rule flags, a slider that reads
the value at any point, a table of flagged stretches, and the curve in 3D.

The project was scaled down to this after discussion with the supervisor (October
2026). Earlier parts, a simulated sequence of runs with assumed defects, planted
faults, a point-density study and a multi-step dashboard, were removed from the
repository; they remain in the git history before that date.

**Why these three.** Gauge, wire height and stagger are all measured inside one
cross-section of the track, relative to the rails. They do not depend on how well two
scans are positioned against each other, and they do not depend on the absolute height.
The height in this scan is relative to the sensor, which rides on the rails, so
height-based track parameters (cross-level, cant, twist, long-wave surface) cancel out
and cannot be read from it; the wire's height above the rails can, because both are in
the same frame.

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
python contact_wire.py run1                                              # the contact wire: height above the rail and stagger
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
| `contact_wire.csv`, `contact_wire_dense.csv` | The contact wire per metre (height above the rail, stagger, sideways position, points used, a flag for two wires side by side, the messenger wire's height where seen) and as a Catmull-Rom curve every 0.1 m with E, N, Z |
| `features.csv` | One row per metre: gauge, 10 m alignment versine, points per metre, rail head above ballast and above sleeper, head tilt, sleeper spacing, wire height and stagger, rail-head intensity left and right |
| `curve_attributes.csv`, `sleepers.csv`, `shortwave.csv` | Optional per-point attributes from `curve_attributes.py` |

## The contact wire

`contact_wire.py` finds the overhead contact wire in the same scan. Points within
1.3 m sideways of the reference and 4.3 to 6.6 m above the rail-head plane are binned
per metre; in every metre the height clusters are found, and the wire is tracked along
the track by taking, in each metre, the cluster nearest to where the wire was in the
previous metre (droppers, fittings and the messenger wire lie higher; wires more than
0.9 m from the track centre are not the contact wire). The height is the median of the
cluster's points above the rail plane; the stagger is the median sideways position
minus the track centreline from `curve.csv`. A metre is flagged when two wires run side
by side (an overlap) or a fitting sits in the cluster. A centripetal Catmull-Rom curve
through the per-metre points gives the wire every 0.1 m for the 3D view.

## The monitor page

```
python monitor.py run1          # writes output/viewer/monitor.html and opens it
```

One self-contained page (works offline), with a tab per quantity and a tab for the
track states. `monitor.html#wire_height`, `#wire_stagger` and `#states` open a tab
directly. Everything is computed in the page, so the controls react at once.

**On each quantity tab:**

- **Reference limits** in a menu, each with one or more levels; a point is marked by the most severe level it exceeds. Gauge: Trafikverket's requirements for track geometry, TRVINFRA-00013 version 2.0 (2021-10-11): table K11.2, the deviation of a single reading from 1435 mm per speed class H0–H5 with the levels PLAN (planning limit), UH1 and UH2 (maintenance limits; UH2 is a safety requirement) and KRIT (critical, immediate action), and table K11.3, minimums for the 100 m sliding mean (choose the 100 m sliding mean as the points of the chart). H4 (160–200 km/h) is preselected, because the scan lies on Södra stambanan south-west of Höör (inferred from its coordinates, see Results); it can be changed in the menu. The FRA limits per class of track (49 CFR 213.53) are kept for comparison; they are far wider. Wire height: the TSI nominal range 5.00–5.75 m, the TSI maximum 6.50 m, and the Swedish nominal 5.50 m with a maximum of 5.65 m that is quoted from a search result and must be confirmed. Stagger: the TSI maximum lateral deviation of 0.40 m (1600 mm pantograph head) or 0.55 m (1950 mm head), both for the wire under cross wind, so the wind allowance comes on top of what the scan shows at rest. Points outside the limits are red; the stretches are shaded red.
- **Stretches, not points.** A fault is a stretch, so the gauge is averaged over 1 m by default (2 m and 5 m available, as are the raw slabs and the 1 m control points); the wire is measured per metre and can be averaged over 2 or 5 m. Alarms shorter than a chosen length are ignored.
- **Western Electric (WECO) rules.** The series is a control chart with zones C (within 1σ of the centre line), B (1σ–2σ) and A (2σ–3σ). Centre line: the track's mean or the nominal value. σ: the scatter of the points, a robust estimate, or the point-to-point noise alone (the raw neighbour-to-neighbour scatter divided by √n for an n-point stretch). Rule 1: one point beyond zone A. Rule 2: two of three in zone A or beyond, same side. Rule 3: four of five in zone B or beyond, same side. Rule 4: eight in a row on one side. All four are on by default; the header shows what each gives; the markers are coloured by the most severe rule that fires. An alarm means "this stretch differs from the rest of this track", not that a limit is broken.
- **Slider and readout** for any chainage: the value, its deviation from nominal, its distance from the centre line in σ, the margin to the nearest limit, a status (OK, WECO alarm with the rules, outside the limits) and, for the wire, whether two wires are present there. Clicking the chart or a table row moves the cursor.
- **The curve in 3D.** The gauge tab embeds the 3D viewer (`view_curve.py`); the wire tabs show the rails and the wire together in the track frame.

**On the track-states tab:** every metre's feature vector (from `features.csv`) is
standardised and grouped by k-means into 2 to 5 states; the features used can be
switched on and off. The page shows the state of every metre along the track, each
feature coloured by state, a profile table (mean of each feature per state and where
the state lies), and the metres furthest from the centre of their own state with the
feature responsible. k-means knows nothing about tracks: a state is a combination of
feature levels that occur together along this track. The clustering runs with fixed
seeds, so it gives the same result every time.

## The 3D view

```
python view_curve.py run1            # writes output/viewer/curve3d_run1.html and opens it
```

One self-contained page: drag to rotate, scroll to zoom. It shows the Catmull-Rom
curve of each rail, the 1 m control points, the centreline and a cross-link every
5 m; **Colour the rails by gauge** paints the curve by the gauge; the **Window** zooms
to a stretch and the **Cursor** reads out gauge and offsets along the track.
**Contact wire** (when `contact_wire.py` has run) draws the wire along the top of the
box: its sideways zigzag is true and follows the sideways stretch, its height above
the rail is the colour, because the real 5.5 m cannot share a height axis that is
stretched a thousandfold for the rails. The wire tabs of the monitor page show it at
true height instead. The viewer also writes `output/viewer/<run>_curve.las`, the curve
as points in the scan's own coordinates, for CloudCompare on top of
`gps_pointcloud_00.laz`.

## Results (one run: `gps_pointcloud_00.laz`)

- 15.9 million points over 264 m of double track, with a gentle curve at the start. The LiDAR followed the left track, so only that track is used. The only per-point attribute with values is the intensity; return number, classification and the rest are empty.
- **Rails and gauge:** both gauge faces are found in 1035 of the 1057 slabs (98 %). Mean gauge 1432.3 mm, 2.7 mm below nominal, scatter 1.6 mm; tightest slab 1424.7 mm, widest 1436.1 mm. The whole track is tight: 97 % of the slabs lie below 1435 mm. Whether that is the track or a millimetre of bias in the gauge-face detection cannot be told from one scan; the differences along the track are unaffected by a constant offset.
- **A step at 155 m (about km 566.29):** the 5 m mean gauge falls from 1434.9 mm at 152 m to 1430.6 mm at 160 m. Before it the track averages 1433–1434 mm, after it 1431–1432 mm. The rail head also sits higher above the sleepers after it (184.6 mm over 130–154 m, 186.4 mm over 156–180 m, 188.3 mm over 184–258 m). The sleeper spacing has a 0.8 m gap at 155–156 m, but also at six other places, so that is not specific to the step; the left rail-head intensity rises locally after the step but is lower again after 184 m. Something physical sits there (a joint, a change of rail or fastening); nothing is mapped in OpenStreetMap at that point. It is the first thing to ask about.
- **Where the scan is (inferred, to be confirmed):** the file has no coordinate reference system, but its coordinates fit SWEREF 99 TM and place the scan on Södra stambanan between Stehag and Höör, 2 km south-west of Höör station, at km 566.44–566.18; the scan ends at Höör's distant signals Hö 22/Hö 52, which OpenStreetMap puts at km 566.180. OSM tags the stretch 180 km/h (Trafikverket lists the line at 200 km/h), so the speed class is H4 (160–200 km/h). The 155 m step is at about km 566.29 and the tight stretch 186–216 m at km 566.23–566.26.
- **Against Trafikverket's maintenance limits (TRVINFRA-00013, table K11.2, point faults):** the wide side is never reached. On the tight side, for the line's class H4, 103 of 260 one-metre stretches are beyond PLAN (below 1432 mm), 46 beyond UH1 (1431 mm) and 17 beyond UH2 (1430 mm), none beyond KRIT (1428 mm); the tight stretch 186–216 m holds most of the UH2 metres. For class H2 (80–120 km/h) PLAN and UH1 coincide at 1431 mm, so the counts are 46/46/17/0; for H5 they are 103/103/46/17. On the raw 25 cm slabs for H4 the counts are 420, 214, 88 and 6 (the KRIT slabs are at the scan edge and in the tight stretch). **Table K11.3, 100 m sliding mean:** the mean runs from 1431.1 to 1433.5 mm, within the H1–H2 limits, below UH1min for H3 (1432 mm, 55 of 165 windows), H4 (1433 mm, 94 windows) and H5 (1434 mm, every window), never below UH2min. So for its own class the stretch is a maintenance case: 17 metres past UH2 to repair and the gauge to be planned for correction over most of its length. A constant bias of a millimetre in the gauge-face detection would move these counts (H4, ±1 mm: UH2 between 3 and 46, KRIT 0 to 3), so the first check is against a track-recording-car run from Trafikverket's Optram system.
- **Against the FRA safety limits:** no slab lies outside any class's limits; the tightest point is 2.3 mm above the 1422.4 mm minimum.
- **Western Electric rules on 1 m stretches** (centre line = the track's mean, σ = 1.45 mm, all four rules, alarms shorter than 1 m ignored): 117 of 260 stretches fire in 12 alarmed stretches. Two are long: 82–132 m is wide (+1 to +3 mm above the mean) and 186–216 m is tight (down to 1428 mm). On raw slabs with rules 1 and 2 only, 28 slabs fire in 11 short stretches.
- **Noise floor of the pipeline** (two complementary halves of the scan compared): gauge 0.9 mm RMS, 95 % of the differences below 1.7 mm per slab; slab-to-slab noise 0.68 mm. A gauge change smaller than that cannot be told from the method's own scatter; positioning error between real scans comes on top.
- **Systematic error, bounded without a reference measurement** (`python error_checks.py`, rebuilds the run four times): sampling the gauge face 10–18 mm below the rail top raises the mean gauge by 0.37 mm; 12–16 mm by 0.84 mm (but finds only 450 of 1057 slabs); 8–24 mm lowers it by 0.37 mm. The band alone moves the mean by 0.4–0.8 mm, so the pipeline's own offset against the standard's band (the smallest distance anywhere from the rail top down to 14 mm below it, which lies higher on the rail than any band the scan can resolve) is of that order; the scanner's calibration is unknown. A constant offset of ±1 mm moves the class H4 count of 1 m stretches beyond UH2 between 46 and 3 (beyond PLAN between 172 and 46, beyond UH1 between 103 and 17); KRIT (1428 mm) is reached by 3 metres only if the reading is 1 mm too wide, and by 17 at 2 mm. The differences along the track are unaffected by an offset.
- **Contact wire:** found in 249 of 264 metres. Height above the rail 5.52–6.05 m, median 5.65 m; stagger −0.58 to +0.59 m. The height is above the TSI nominal upper bound of 5.75 m over 90–105 m and from 223 m to the end (up to 6.05 m), and the stagger is beyond ±0.40 m over 90–109 m and from 236 m to the end (48 metres; 12 beyond ±0.55 m, at 90–95 m and 258 m to the end). Two wires run side by side in 101 metres, mostly between 107 and 220 m. The excursions sit just before and just after that stretch, which is what an overlap section looks like (the outgoing wire rises and moves aside while the next takes over; the wire tracked at the end is probably the one going to its anchor), so they are most likely geometry of the overhead line rather than faults. The tracking does not tell which of the two wires is the running one inside the overlap. Nothing exceeds the TSI maximum of 6.50 m.
- **Track states** (k = 3, all 11 features): three contiguous stretches, 0–64 m, 69–170 m and 184–259 m, with the boundaries where the gauge, the rail-head intensity and the wire change together. The split near 155–170 m matches the gauge step. No metre scores above 2σ from its own state's centre.

## Assumptions and limits

- One scan only; nothing about change over time can be shown yet. The pipeline is built for a second run, which would be processed with the same reference and compared point for point.
- The gauge is judged against Trafikverket's maintenance limits (TRVINFRA-00013) and, for comparison, the FRA safety limits; the wire against the TSI values. The line and its speed class (Södra stambanan, H4) are inferred from the file's coordinates and OpenStreetMap, not stated in the data, so the Trafikverket tables are applied per class and H4 is only preselected. The Swedish wire maximum and the stagger limit for 1800 mm pantograph heads are to be confirmed.
- Trafikverket defines the gauge as the smallest distance between the rails measured anywhere from the rail top down to 14 mm below it, unfiltered (TRVINFRA-00013, terms and section 8.2.1); the pipeline samples the gauge face between 8 and 20 mm below the top, a little lower than that.
- Height is relative to the sensor, so of the track parameters only gauge and lateral alignment are measured in full; the wire is measured against the rails and is unaffected.
- The reference covers the first run's extent plus 2 m at each end.

## Background: gauge, TQI and the source paper

Sharma et al. (2018) work with five track geometry parameters (gauge, cross-level,
cant, twist, surface), define the Track Quality Index per parameter and section as
`TQI = (Ls / L0 − 1) × 10^6` (the measured trace length against the section length),
and treat geometry measurements that violate the FRA thresholds as "red-tag" defects
that must be repaired before the next inspection. This project uses their definition
of gauge and the FRA gauge limits they refer to. `track_curve.py sections` still
computes the TQI per 25 m section (`sections.csv`), with 0.25 m samples in metres, so
the values are not comparable with the paper's foot-by-foot ones.

## Repository contents

| File | Purpose |
|---|---|
| `track_curve.py` | The pipeline: reference frame, per-run curve and gauge, sections and TQI, comparison of runs, self-test |
| `catmull_rom.py` | Centripetal Catmull-Rom curve math used by the pipeline |
| `contact_wire.py` | The overhead contact wire from the same scan: height above the rail and stagger per metre, as a curve |
| `track_features.py` | One feature vector per metre from the pipeline's outputs and the rail-head intensity |
| `error_checks.py` | The two checks of the systematic error without a reference measurement: the gauge-face band and the offset sensitivity of the Trafikverket verdict |
| `report_figures.py` | The report's figures from the outputs (`output/figures/r1`–`r4`) |
| `monitor.py`, `monitor_template.html` | The page: gauge, wire height and stagger against limits with Western Electric rules, slider and readout, flagged stretches, 3D, and the track states |
| `view_curve.py`, `view_curve_template.html` | The 3D viewer and the CloudCompare export |
| `curve_attributes.py` | Optional per-point attributes of the curve (normals, ballast, sleepers, short-wave height indicators); some of them are features |
| `gps_pointcloud_00.laz` | The scan: LiDAR point cloud of a two-track line (LAS 1.2, point format 0, no timestamps or coordinate system) |
| `output/reference/track.json` | The reference frame every run is projected onto |
| `output/runs/run1/` | The curve and measurements of run 1 (see the table above) |
| `output/noise_floor.json` | Noise floor and known-change check from `selftest` |
| `output/viewer/` | `monitor.html`, the 3D view `curve3d_run1.html` and the curve as LAS files for CloudCompare |
| `output/figures/` | Figures made for the slides and the report (rail cross-section, curve through control points, gauge per slab and section) |
| `D7016B_Demo1_Viktor_Lindberg.pptx`, `D7016B_Demo2_Viktor_Lindberg.pptx` | Slides for Demonstrations 1 and 2 (the Demo 2 deck predates the scale-down) |
| `Industrial AI - Part II .docx` | Final report, based on the course report template. Local drafts of it (`Report rough draft.docx`) are ignored by git, see `.gitignore` |
| `PointCloud_Project_Time_Plan_Viktor - Rough_Final (2).xlsx` | The original time plan |
| `1-s2.0-S0968090X18302444-main.pdf` | Sharma et al. (2018), the source for the geometry parameters, the TQI and the FRA threshold reference |

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

Federal Railroad Administration. 49 CFR § 213.53 Gage. Track Safety Standards,
Subpart C, Track Geometry. https://www.ecfr.gov/current/title-49/section-213.53

European Commission (2014). Commission Regulation (EU) No 1301/2014 on the technical
specifications for interoperability relating to the 'energy' subsystem of the rail
system in the Union, points 4.2.9.1 (contact wire height) and 4.2.9.2 (lateral
deviation of the contact wire). https://eur-lex.europa.eu/eli/reg/2014/1301/oj

Western Electric Company (1956). Statistical Quality Control Handbook. Indianapolis.
