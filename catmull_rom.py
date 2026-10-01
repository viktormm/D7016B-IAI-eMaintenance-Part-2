"""Centripetal Catmull-Rom curves through control points. Pure curve math, no I/O.

The curve passes through every control point, and moving one control point only changes the curve nearby.
The centripetal parameterisation (alpha = 0.5) avoids loops and overshoot when the points are unevenly spaced.
It is invariant to rotation and translation, so it gives the same curve in any frame.
"""
import numpy as np


def catmull_rom(P, per_segment=20, alpha=0.5):
    """Dense samples of the curve through the rows of P (n x k). Passes through every row of P."""
    P = np.asarray(P, float)
    ext = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])
    t = np.zeros(len(ext))
    t[1:] = np.cumsum(np.linalg.norm(np.diff(ext, axis=0), axis=1) ** alpha)
    out = []
    for i in range(1, len(ext) - 2):
        p0, p1, p2, p3 = ext[i - 1:i + 3]
        t0, t1, t2, t3 = t[i - 1:i + 3]
        tt = np.linspace(t1, t2, per_segment, endpoint=False)[:, None]
        a1 = (t1 - tt) / (t1 - t0) * p0 + (tt - t0) / (t1 - t0) * p1
        a2 = (t2 - tt) / (t2 - t1) * p1 + (tt - t1) / (t2 - t1) * p2
        a3 = (t3 - tt) / (t3 - t2) * p2 + (tt - t2) / (t3 - t2) * p3
        b1 = (t2 - tt) / (t2 - t0) * a1 + (tt - t0) / (t2 - t0) * a2
        b2 = (t3 - tt) / (t3 - t1) * a2 + (tt - t1) / (t3 - t1) * a3
        out.append((t2 - tt) / (t2 - t1) * b1 + (tt - t1) / (t2 - t1) * b2)
    out.append(P[-1:])
    return np.vstack(out)


def evaluate(s, values, s_grid, max_gap=None):
    """Curve through the points (s, values), evaluated at s_grid.

    NaN points are ignored. The result is NaN outside the range of the points, and wherever two neighbouring
    points are further apart than max_gap (so a missing stretch is not bridged silently).
    """
    s, values = np.asarray(s, float), np.asarray(values, float)
    ok = np.isfinite(values) & np.isfinite(s)
    s, values = s[ok], values[ok]
    out = np.full(len(s_grid), np.nan)
    if len(s) < 3:
        return out
    dense = catmull_rom(np.column_stack([s, values]))
    inside = (s_grid >= s[0]) & (s_grid <= s[-1])
    out[inside] = np.interp(s_grid[inside], dense[:, 0], dense[:, 1])
    if max_gap is not None:
        i = np.clip(np.searchsorted(s, s_grid) - 1, 0, len(s) - 2)
        out[(s[i + 1] - s[i]) > max_gap] = np.nan
    return out
