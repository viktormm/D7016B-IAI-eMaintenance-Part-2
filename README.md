# Following a Railway Track Over Time with 3D Point Cloud Curves

**Course:** D7016B – Industrial AI and eMaintenance, Part II — Luleå University of Technology
**Author:** Viktor Lindberg
**Supervisor:** Amit Patwardhan

## What this project is about

The project follows how a railway track changes over time. Every LiDAR run over
the same track is turned into a 3D curve, and the runs are compared to see where
faults appear and how they grow. Catmull-Rom splines are used as the curve
representation.

**Position and height.** A run's curve is `(s, d, z)`: chainage `s` along a fixed
reference line, sideways offset `d` from that line, and the height `Z` exactly as
stored in the LiDAR file. `Z` appears to be relative to the sensor, which rides on
the track (an inference from the data, not confirmed), so no absolute heights are
needed.

**Why the direction changed.** The original idea was to link the D7015B vibration
events (150 events from three switch measurement trials, classified normal/abnormal
with an SVM at 93 % accuracy) to the track geometry at the same positions, cluster
the abnormal events by geometry and predict the fault type from vibration alone.
The vibration data that would link events to the geometry is not available, so that
idea cannot be tested. After talking to Amit the project follows the track over time
instead. The original plan is kept in the Demo 1 slides (slides 1–14 of the deck).

## How a run becomes a repeatable 3D curve

`track_curve.py` turns a `.laz` run into a curve that can be compared with other
runs. What makes it repeatable:

- **One fixed reference frame.** The smoothed track centreline of the first run is
  stored once (`output/reference/track.json`). Every run, including the first, is
  projected onto it: `s` along the reference, `d` sideways from it (left is `+d`, so
  left and right do not depend on the travel direction), `Z` as stored.
- **Fixed grids anchored at `s = 0`.** Slabs (25 cm), control points (1 m) and
  sections (25 m) sit on the same grid in every run, so runs line up point for
  point. Nothing depends on where a particular run happens to start.
- **No random steps and no data-dependent settings.** All settings live in one
  `Params` block and are stored with every run, together with the input file's
  SHA-256, the reference id and a hash of the code. `compare` refuses runs made with
  a different reference or grid.
- **Checked by a self-test** (`selftest`): rebuilding gives byte-identical files, two
  halves of the same scan give the noise floor of the whole pipeline, and a known
  change added to one half is recovered (see the results below).

The curve of a run is:

1. Project the points onto the reference; keep those within 2.5 m of it.
2. In every 25 cm slab, find the rail-head top of each rail (highest points in a window that follows the rail) and the inner (gauge) face 8–20 mm below it.
3. Reject outliers against a running median.
4. Average into control points every 1 m on the fixed grid (a line fitted through the slabs in each bin, evaluated at the bin centre).
5. Join the control points with a centripetal Catmull-Rom curve. It passes through every control point, and a change in one place only moves the curve nearby.

Each rail's lateral position `d` is measured at the inner (gauge) face and its
height `z` at the head top. The face is a sharper edge than the head centre: it cut
the lateral noise from about 2.3 mm to 1.4 mm per rail.

### Running it

From the repo root:

```
python track_curve.py reference gps_pointcloud_00.laz --name track       # once, from the first run
python track_curve.py build gps_pointcloud_00.laz --run-id run1 --reference output/reference/track.json
python track_curve.py sample run1                                        # the 3D curve every 0.1 m
python track_curve.py sections run1                                      # 25 m sections and TQI
python track_curve.py selftest gps_pointcloud_00.laz --reference output/reference/track.json
```

To add a second run and compare it with the first (same reference, same settings):

```
python track_curve.py build run2.laz --run-id run2 --reference output/reference/track.json
python track_curve.py compare run1 run2 --plot
python track_curve.py sections run2 --norm-run run1       # TQI normalised with run 1, so runs stay comparable
```

Output per run in `output/runs/<run-id>/`:

| File | Content |
|---|---|
| `curve.csv` | The curve as control points every 1 m: `s`, `d_left`, `z_left`, `d_right`, `z_right`, `gauge`, `d_centre`, `z_centre`, and `E`/`N` of both rails |
| `curve_dense.csv` | The 3D curve (E, N, Z for both rails and the centreline) sampled every 0.1 m |
| `slabs.csv` | The per-slab measurements the curve is built from |
| `sections.csv` | 25 m sections with gauge statistics, alignment and TQI |
| `meta.json` | Settings, input hash, reference id, code hash and quality figures |

### Viewing the curve in 3D

```
python view_curve.py run1            # writes output/viewer/curve3d_run1.html and opens it in the browser
python view_curve.py run1 run2       # several runs of the same reference in one view
```

`curve3d_<runs>.html` is one self-contained page (it works offline, no server). Drag to rotate, scroll to zoom, right-drag to pan.

- **Curves:** the Catmull-Rom curve of each rail, the 1 m control points it passes through, the centreline (thick, midway between the rails; red for one run, in each run's own colour when several are shown) and a cross-link between the rails every 5 m. The curves are recomputed from the run's `curve.csv`, so the view is never out of date.
- **Track frame** (default): chainage `s`, sideways offset `d` and height `Z`. The rails are 1.43 m apart in a track hundreds of metres long, so sideways and height are stretched with two sliders (×1 is true scale). **World** shows `E`, `N` as in the LiDAR file at true horizontal scale.
- **Colour** the rails by gauge, run or chainage. Coloured by gauge (the default), the 3.0 mm step at 155 m chainage shows as a change from amber (wide) to blue (narrow).
- **Window** (from–to in metres) zooms to a stretch; the box keeps its shape. **Cursor** moves a marker along the track and reads out gauge, offsets and heights per run.
- **Point cloud sample** (off at first, tick the box to show it): 60 000 of the raw points around the rails as context (`--points N`, `0` for none, `--laz FILE` for another file). The sample is thinned evenly, not at random. The box keeps its shape when the cloud is switched on or off, so the height stretch changes with it (about ×1,160 without the cloud, about ×40 with it); next to the cloud's half metre of height the rails look flat.

`python view_curve.py` also writes `output/viewer/<run>_curve.las`: both rails and the centreline every 0.1 m in the point cloud's own coordinates, with the scalar fields `s_m`, `gauge_mm` and `rail` (0 left, 1 right, 2 centre); left is orange, right is teal, centre is grey. Where a rail is missing it has no points for that stretch, and the other rail is still written. To see the curve in the original point cloud, open `gps_pointcloud_00.laz` and this file together in CloudCompare (accept the same global shift for both when it asks) and increase the point size of the curve.

Both outputs are byte-identical when rebuilt from the same run. The viewer only reads the runs; it does not change them.

With three or more runs the view gets a **Runs over time** control: the first run stays on screen as the baseline and one later run is shown against it; the slider or Play steps through the runs. The changes are millimetres on a 1.4 m gauge, so the later run's sideways change against the baseline is drawn larger than it is (×100 at first, slider down to true size), and with **Colour → Change against the first run** its rails are coloured by gauge change and its centreline by sideways change. The cursor table always gives the measured values; its last row is the change between the two runs in mm. The sideways axis is the same for every run, so a zoomed view stays where it is while the runs play.

### Simulated sequence of runs

Only one real run exists, so `simulate_runs.py` shows what the pipeline would report if the track were scanned at regular intervals:

```
python simulate_runs.py gps_pointcloud_00.laz --reference output/reference/track.json     # 13 runs, every 2 weeks
python simulate_runs.py gps_pointcloud_00.laz --reference output/reference/track.json --view-only   # only rebuild the 3D view
```

It is a closed loop, and nothing is scripted to happen in a given week:

1. **The track degrades** by an assumed model: three defects that grow linearly within a maintenance cycle (as in track geometry degradation models, see Soleimanmeigouni et al., 2018) — a gauge widening at 100–120 m (0.25 mm a week), an alignment defect around 210 m (0.35 mm a week) and a short-wave gauge irregularity at 227–248 m (1 m wavelength, 0.10 mm a week).
2. **Every run is scanned and measured**: a different pseudo-random half of the point cloud (its own sampling noise, like a real run) with that week's defects added, built by the unchanged pipeline.
3. **Indicators are checked against limits**, all against the first run: the gauge change and the sideways change of the centreline as a 5 m mean along the curve, and the increase of the TQI (gauge) of every 25 m section. Detection limits are calibrated from four extra runs *without* degradation (the largest value they reach, times 1.2); above them a finding is an **alert**. Intervention limits are assumed (3.0 mm gauge, 3.5 mm alignment, twice the TQI detection limit); at or above them it is an **alarm**.
4. **An alarm orders maintenance** on the alarmed stretch (gauge correction, lining or fastening renewal). It is done before the next run and resets the defects of that kind there, after which they grow again.

`--position-sigma MM` adds a random sideways positioning error per run, which is taken out again (median offset) before the runs are compared.

The runs are written to `output/runs/sim_wNN/` (NN = week; their `meta.json` says they are simulated and lists their findings), the summary to `output/simulation/` (`events.csv`, `limits.csv`, `timeline.csv`, `sections_over_time.csv`, `sim_spacetime.png`, `sim_sections.png`) and the 3D view to `output/viewer/curve3d_simulation.html`, where the alerts and alarms of the run on screen are listed under the run slider and marked under the track; clicking one goes to that stretch.

With the default settings the alignment defect raises an alert in week 4 and alarms in weeks 10 and 20, the gauge widening an alert in week 4 and an alarm in week 14, and the short-wave irregularity a TQI alert in week 16 and an alarm in week 24. The curve finds the first two and not the third; the TQI finds the third and not the first two. Three of the 37 findings are false alerts (short and within 0.25 mm of the limit); none of the four alarms is.

The sequence shows the behaviour of the method, not of the track: the degradation and the intervention limits are assumed, the rates are chosen so that a whole cycle fits in 24 weeks, and all runs come from one acquisition, so the positioning error between real runs is not in it unless `--position-sigma` is used.

Soleimanmeigouni, I., Ahmadi, A. and Kumar, U. (2018). "Track geometry degradation and maintenance modelling: A review". Proceedings of the Institution of Mechanical Engineers, Part F: Journal of Rail and Rapid Transit, 232(1), pp. 73–102.

### Assumptions and limits

- The runs are of the same track and the LiDAR follows that track. A run on the other track is not handled.
- The reference covers the first run's extent plus 2 m at each end. A run that goes further is cut at the reference's ends; extending the reference means making a new one and rebuilding all runs.
- Runs are only as comparable as their coordinates. How well the coordinates of real runs line up is not known yet, since only one run exists.
- `Z` is relative to the sensor, so the height of the LiDAR track itself stays flat (see below).

## Results so far (one run: `gps_pointcloud_00.laz`)

- 15.9 million points over about 260 m of double track, with a gentle curve at the start. The LiDAR followed the left track (higher density, occlusion on both sides of the rails), so only that track is used.
- **Rails and gauge:** both rails and both gauge faces are found in 98 % of the 25 cm slabs (the rest are the 2 m margins of the reference and a few gaps). Gauge is 1432.3 mm on average with 1.6 mm scatter. There is a step of 3.0 mm in the gauge at 155 m chainage that is not explained.
- **Curve:** 261 control points per rail with both rails found.
- **Noise floor of the pipeline** (two complementary halves of the scan compared, at half the point density): lateral position 1.4–1.5 mm RMS per rail (95 % of the change below about 3 mm), height 0.3–0.4 mm RMS, gauge 0.9 mm RMS (95 % below 1.7 mm). Positioning error between real runs comes on top of this.
- **Known change recovered:** a 4 mm gauge widening added to one half of the scan is recovered as 3.9 mm and found across the whole stretch; a 5 mm height shift is recovered as 5.1 mm; false alarms are about 5 %, as expected from a 95 % threshold.
- **Condition table:** ten 25 m sections with gauge and alignment statistics and TQI (`output/runs/run1/sections.csv`). The curved start of the scan ranks worst.

## What can be tracked over time

| Parameter | Over time? | Note |
|---|---|---|
| Gauge | Yes | Measured inside one cross-section, so it does not depend on positioning between runs |
| Alignment | Yes, relative | Sideways offset from the shared reference frame; needs consistent positioning between runs |
| Height (`Z`) | Yes, relative | The LiDAR height as it is, relative to the sensor. The LiDAR track's own rails stay flat (within 3 mm over 260 m), so height change shows against the neighbouring track and the ballast, not as settlement of the LiDAR track itself |
| Cross-level, twist, surface | No | The rails keep the same height and roll relative to the sensor, so these are cancelled out; computed values are flat at the 1 mm resolution of the data |
| Cant | No | Needs cross-level and the designed values |

The neighbouring track's rails move about 28 cm relative to the sensor over 90 m
while the LiDAR track's own rails move 3 mm. That contrast is the evidence that `Z`
is relative to the LiDAR track.

## Background: the five parameters and TQI

Based on Sharma et al. (2018):

| # | Parameter   | What it measures                                   |
|---|-------------|-----------------------------------------------------|
| 1 | Gauge       | Distance between the two rails                       |
| 2 | Cross-level | How much one rail sits above the other               |
| 3 | Cant        | The cross-level a curve is designed to have           |
| 4 | Twist       | How much cross-level changes over a short distance    |
| 5 | Surface     | How even the top of the rail is                       |

The Track Quality Index per parameter and section is

```
TQI = (Ls / L0 − 1) × 10^6
```

where `L0` is the section's nominal length and `Ls` is the length of the actual
(wobbly) trace of that parameter along the section. A perfectly uniform section has
`TQI = 0`; every deviation makes the trace longer and raises the TQI.

In Sharma et al. (2018) the trace length is `Ls = Σ √(Δx² + Δy²)` over foot-by-foot
samples (Δx = 1 ft, Δy in ft) for 0.1-mile (528 ft ≈ 161 m) segments. Each parameter's
TQI is divided by its 95th percentile over the whole dataset and the parameters are
averaged. Absolute TQI values depend on the sampling interval and units, so only
compare TQI between sections computed the same way. Here TQI uses 0.25 m samples in
metres for gauge and alignment, so the values are not comparable with the paper's.
To keep runs comparable over time the normalisation is fixed from the first run
(`sections --norm-run run1`) instead of being redone for every run. The paper tracks
TQI over 13 inspection runs, which is the same idea this project applies to point
cloud runs.

## Open questions for Amit

- Are there more runs of this track, and how well do their coordinates line up?
- How repeatable are the positions between runs?
- Is the height relative to the sensor, as it looks, and do all runs use the same reference?

## Next steps

- Build a second run with the same reference and compare it with run 1.
- Check how well the coordinates of real runs line up, and use that with the noise floor to set change thresholds.
- Follow TQI per section over the runs.

## Repository contents

| File | Purpose |
|---|---|
| `track_curve.py` | The pipeline: reference frame, per-run curve, sections and TQI, comparison of runs, self-test |
| `catmull_rom.py` | Centripetal Catmull-Rom curve math used by the pipeline |
| `view_curve.py`, `view_curve_template.html` | The 3D viewer and the CloudCompare export (kept apart from the pipeline, so they do not change the code hash of the runs) |
| `simulate_runs.py` | A simulated sequence of runs from the one scan, with an assumed degradation, to show what the comparison would report over time |
| `gps_pointcloud_00.laz` | The first run: LiDAR point cloud of a two-track line (LAS 1.2, point format 0, no timestamps or coordinate system) |
| `output/reference/track.json` | The reference frame every run is projected onto |
| `output/runs/run1/` | The curve and measurements of run 1 (see the table above) |
| `output/noise_floor.json` | Noise floor and known-change test from `selftest` |
| `output/viewer/` | The 3D view (`curve3d_<runs>.html`) and the curve as a LAS file for CloudCompare (`<run>_curve.las`) |
| `output/figures/` | Figures made for the slides and the report. The scripts that drew them are not kept in the repo |
| `D7016B_Demo1_Viktor_Lindberg.pptx` | Slides: Demonstration 1 (slides 1–14, original idea) and the Demonstration 2 part (slides 15–23, follow the track over time) |
| `D7016B_Demo2_Viktor_Lindberg.pptx` | Slides for Demonstration 2 on their own, six slides: change of idea and the data, method, test on the real scan, simulated monitoring, where it stands |
| `Industrial AI - Part II .docx` | Final report, based on the course report template |
| `PointCloud_Project_Time_Plan_Viktor - Rough_Final.xlsx` | The original time plan, not updated for the new direction |
| `1-s2.0-S0968090X18302444-main.pdf` | Sharma et al. (2018), the reference paper for the geometry parameters and TQI |

## Tools

Python (needs `numpy`, `scipy`, `matplotlib`, `plotly` and `laspy[lazrs]`; the `lazrs`
backend is required to read `.laz` files, and `plotly` is only used to embed plotly.js in the 3D view). CloudCompare (`winget install CloudCompare.CloudCompare`)
is used to inspect the point cloud visually.

## Reference

Sharma, S., Cui, Y., He, Q., Mohammadi, R. and Li, Z. (2018). "Data-driven
optimization of railway maintenance for track geometry". Transportation
Research Part C, 90, pp. 34–58. https://doi.org/10.1016/j.trc.2018.02.019
(`1-s2.0-S0968090X18302444-main.pdf` in this repo) — defines the five track
geometry parameters and the Track Quality Index used in this project.
