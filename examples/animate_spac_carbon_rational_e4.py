"""Animate E4 trajectories and every admissible positive frozen equilibrium."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import animation
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

from control_carbon.spac_carbon import YEAR_DAYS
from control_carbon.spac_carbon_rational import RationalSPAC
from analyze_spac_carbon_rational_e4_diagnostics import (
    EARLY_YEARS, ROOT, STABLE_COLOR, UNSTABLE_COLOR, VPD_HOLD_DAYS,
    VPD_PEAK, VPD_START, pulse_vpd, spectral_crossings,
)
from prepare_manuscript_phase_portraits import model_family


ACTUAL_COLOR = "#009E73"
NEUTRAL_COLOR = "#777777"


def frame_days():
    """Use the same nonuniform physical clock in both movies."""
    knots = [0., 10., 15., 10.+VPD_HOLD_DAYS, 20.+VPD_HOLD_DAYS,
             10*YEAR_DAYS, 11*YEAR_DAYS, 21*YEAR_DAYS, EARLY_YEARS*YEAR_DAYS]
    counts = [45, 60, 35, 45, 80, 30, 80, 45]
    return np.r_[np.concatenate([
        np.linspace(a, b, n, endpoint=False)
        for a, b, n in zip(knots[:-1], knots[1:], counts)
    ]), knots[-1]]


def interpolate_state(days, states, day):
    if not np.isfinite(day) or day < days[0] or day > days[-1]:
        raise ValueError("Cannot animate a state outside the integrated time interval")
    water = [np.interp(day, days, states[:, j]) for j in range(2)]
    carbon = [np.exp(np.interp(day, days, np.log(states[:, j]))) for j in range(2, 6)]
    return np.r_[water, carbon]


def frozen_equilibria(model, days, ramp_days):
    records = []
    for day in days:
        vpd = pulse_vpd(day, ramp_days)
        equilibria = []
        for state in model.equilibria(vpd):
            leading = float(np.linalg.eigvals(model.jacobian(state, vpd)).real.max())
            if np.max(np.abs(model.rhs(state, vpd))) > 1e-11:
                raise RuntimeError("Movie equilibrium fails the full six-state balance")
            equilibria.append({
                "state": state.tolist(), "leading_real_per_day": leading,
                "stability": "stable" if leading < 0 else "unstable",
            })
        records.append({"day": float(day), "vpd_mpa": vpd,
                        "positive_equilibria": equilibria})
    return records


def colored_curve(ax, coordinates, stable, *, three_d=False, alpha=1., lw=2.):
    segments = np.stack([coordinates[:-1], coordinates[1:]], axis=1)
    colors = np.where(stable[:-1] & stable[1:], STABLE_COLOR, UNSTABLE_COLOR)
    cls = Line3DCollection if three_d else LineCollection
    collection = cls(segments, colors=colors, linewidths=lw, alpha=alpha)
    if three_d:
        ax.add_collection3d(collection)
    else:
        ax.add_collection(collection)
    return collection


def build_movie(model, data, records, days, ramp_days, label, output, fps):
    states, integration_days = data["state"], data["day"]
    early = integration_days <= days[-1]
    path_states, path_days = states[early], integration_days[early]
    forcing_days = np.unique(np.r_[days, np.linspace(0, days[-1], 2500)])
    forcing = np.array([pulse_vpd(day, ramp_days) for day in forcing_days])
    curve_vpd = np.linspace(VPD_START, VPD_PEAK, 401)
    curves = [model.equilibria(vpd) for vpd in curve_vpd]
    if any(len(row) != 1 for row in curves):
        raise RuntimeError("This E4 baseline branch reference expects exactly one positive root")
    curve_states = np.array([row[0] for row in curves])
    curve_stable = np.array([
        np.linalg.eigvals(model.jacobian(state, vpd)).real.max() < 0
        for state, vpd in zip(curve_states, curve_vpd)
    ])
    _, _, crossings = spectral_crossings(model)
    relevant_crossings = [x for x in crossings if VPD_START <= x["vpd_mpa"] <= VPD_PEAK]
    if relevant_crossings:
        curve_vpd = np.sort(np.r_[curve_vpd, [x["vpd_mpa"] for x in relevant_crossings]])
        curve_states = np.array([model.equilibria(v)[0] for v in curve_vpd])
        curve_stable = np.array([
            np.linalg.eigvals(model.jacobian(y, v)).real.max() < 0
            for y, v in zip(curve_states, curve_vpd)
        ])
    fig = plt.figure(figsize=(12, 7.6), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, width_ratios=(1., 1.25))
    time_axes = [fig.add_subplot(grid[row, 0]) for row in range(3)]
    ax3 = fig.add_subplot(grid[:2, 1], projection="3d")
    phase = fig.add_subplot(grid[2, 1])
    time_text = fig.suptitle("", fontsize=12)
    time_axes[0].plot(forcing_days/YEAR_DAYS, forcing*1000, color=".25", lw=1.3)
    time_axes[0].set_ylabel("VPD (kPa)")
    time_axes[0].set_ylim(2.40, 3.05)
    cursors = [ax.axvline(0., color=".45", lw=.8) for ax in time_axes]
    actual_time_lines = [time_axes[j].plot([], [], color=ACTUAL_COLOR, lw=1.7)[0]
                         for j in (1, 2)]
    for ax in time_axes:
        ax.set(xlim=(0, EARLY_YEARS), xlabel="Time (years)")
        ax.grid(alpha=.18)
    time_axes[1].set(ylabel=r"$\psi_l$ (MPa)", ylim=(-4.1, -.5))
    time_axes[2].set(ylabel=r"$C_L$ (kg C m$^{-2}$)", ylim=(-.004, .078))
    max_roots = max(len(rec["positive_equilibria"]) for rec in records)
    for root_index in range(max_roots):
        leaf = np.array([
            rec["positive_equilibria"][root_index]["state"][2]
            if root_index < len(rec["positive_equilibria"]) else np.nan for rec in records
        ])
        psi = np.array([
            rec["positive_equilibria"][root_index]["state"][1]
            if root_index < len(rec["positive_equilibria"]) else np.nan for rec in records
        ])
        stable = np.array([
            root_index < len(rec["positive_equilibria"])
            and rec["positive_equilibria"][root_index]["stability"] == "stable"
            for rec in records
        ])
        for ax, values in zip(time_axes[1:], [psi, leaf]):
            colored_curve(ax, np.column_stack([days/YEAR_DAYS, values]), stable,
                          alpha=.65, lw=1.4)
    actual3, = ax3.plot([], [], [], color=ACTUAL_COLOR, lw=1.9)
    actual2, = phase.plot([], [], color=ACTUAL_COLOR, lw=1.7)
    dot3, = ax3.plot([], [], [], "o", color=ACTUAL_COLOR, markersize=5)
    dot2, = phase.plot([], [], "o", color=ACTUAL_COLOR, markersize=5)
    colored_curve(ax3, np.column_stack([curve_vpd*1000, curve_states[:, 1],
                                       curve_states[:, 2]]), curve_stable,
                  three_d=True, alpha=.5)
    colored_curve(phase, curve_states[:, [1, 2]], curve_stable, alpha=.5)
    eq3, eq2 = [], []
    for _ in range(max_roots):
        eq3.append(ax3.plot([], [], [], marker="*", ls="None", markersize=11)[0])
        eq2.append(phase.plot([], [], marker="*", ls="None", markersize=10)[0])
    carbon_neutral_psi = 4*(model.opening(0.)-1.)
    for ax in (phase,):
        ax.plot([-4.1, carbon_neutral_psi], [0, 0], color=NEUTRAL_COLOR, lw=4, alpha=.6)
        ax.plot([carbon_neutral_psi, 0], [0, 0], color=UNSTABLE_COLOR, lw=4, alpha=.6)
    zero3, = ax3.plot([], [], [], color=NEUTRAL_COLOR, lw=3)
    zero3_unstable, = ax3.plot([], [], [], color=UNSTABLE_COLOR, lw=3)
    ax3.set(xlabel="VPD (kPa)", ylabel=r"$\psi_l$ (MPa)",
            zlabel=r"$C_L$ (kg C m$^{-2}$)", xlim=(2.4, 3.05),
            ylim=(-4.1, 0.), zlim=(-.004, .078))
    ax3.view_init(elev=24, azim=-58)
    phase.set(xlabel=r"$\psi_l$ (MPa)", ylabel=r"$C_L$ (kg C m$^{-2}$)",
              xlim=(-4.1, 0.), ylim=(-.004, .078))
    phase.grid(alpha=.18)
    handles = [
        Line2D([], [], color=ACTUAL_COLOR, label="Actual path"),
        Line2D([], [], color=STABLE_COLOR, marker="*", label="Stable positive equilibrium"),
        Line2D([], [], color=UNSTABLE_COLOR, marker="*", label="Unstable equilibrium"),
        Line2D([], [], color=NEUTRAL_COLOR, lw=3,
               label="Zero-C continuum (non-hyperbolic; not a stable point)"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=8)
    zoom = time_axes[0].inset_axes([.12, .47, .4, .43])
    zoom_days = np.linspace(0, 10, 201)
    zoom.plot(zoom_days, [1000*pulse_vpd(t, ramp_days) for t in zoom_days], color=".25")
    zoom.set(xlim=(0, 10), ylim=(2.4, 3.05), xticks=(0, 10), yticks=(2.44, 3.),
             xlabel="First 10 days")
    zoom.tick_params(labelsize=7)
    zoom.xaxis.label.set_size(7)
    zoom_cursor = zoom.axvline(0., color=ACTUAL_COLOR, lw=1.)
    stopped_note = phase.text(.02, .96, "", va="top", transform=phase.transAxes, fontsize=8)

    def update(index):
        day = days[index]
        year = day/YEAR_DAYS
        vpd = pulse_vpd(day, ramp_days)*1000
        for cursor in cursors:
            cursor.set_xdata([year, year])
        zoom_cursor.set_visible(day <= 10.)
        zoom_cursor.set_xdata([day, day])
        inside = day <= integration_days[-1]
        last_day = min(day, float(integration_days[-1]))
        end_state = interpolate_state(integration_days, states, last_day)
        mask = path_days < last_day
        shown_days = np.r_[path_days[mask], last_day]
        shown_states = np.vstack([path_states[mask], end_state])
        shown_vpd = np.array([1000*pulse_vpd(t, ramp_days) for t in shown_days])
        actual3.set_data_3d(shown_vpd, shown_states[:, 1], shown_states[:, 2])
        actual2.set_data(shown_states[:, 1], shown_states[:, 2])
        for line, column in zip(actual_time_lines, (1, 2)):
            line.set_data(shown_days/YEAR_DAYS, shown_states[:, column])
        dot3.set_visible(inside)
        dot2.set_visible(inside)
        dot3.set_data_3d([vpd], [end_state[1]], [end_state[2]])
        dot2.set_data([end_state[1]], [end_state[2]])
        stop_message = "" if inside else (
            f"Integration stopped at {integration_days[-1]/YEAR_DAYS:.2f} years;\n"
            "no trajectory is extrapolated.")
        stopped_note.set_text(stop_message)
        positive = records[index]["positive_equilibria"]
        for root_index, (marker3, marker2) in enumerate(zip(eq3, eq2)):
            exists = root_index < len(positive)
            marker3.set_visible(exists)
            marker2.set_visible(exists)
            if exists:
                eq = positive[root_index]
                psi, leaf = eq["state"][1:3]
                color = STABLE_COLOR if eq["stability"] == "stable" else UNSTABLE_COLOR
                marker3.set_color(color)
                marker2.set_color(color)
                marker3.set_data_3d([vpd], [psi], [leaf])
                marker2.set_data([psi], [leaf])
        zero3.set_data_3d([vpd, vpd], [-4.1, carbon_neutral_psi], [0, 0])
        zero3_unstable.set_data_3d([vpd, vpd], [carbon_neutral_psi, 0.], [0, 0])
        time_text.set_text(
            f"E4, {label} | t = {day:.2f} days ({year:.4f} years), "
            f"VPD = {vpd:.3f} kPa\n"
            "Shared nonuniform playback clock; all plotted time axes are linear")

    stem = "e4_fast_10day" if ramp_days == 10. else "e4_slow_10year"
    movie_path = output/f"{stem}.mp4"
    writer = animation.FFMpegWriter(fps=fps, codec="libx264", bitrate=2200,
                                    extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    with writer.saving(fig, str(movie_path), dpi=100):
        for index in range(len(days)):
            update(index)
            writer.grab_frame()
    for name, index in [("start", 0), ("peak", int(np.flatnonzero(days == 10.)[0])
                       if ramp_days == 10. else int(np.flatnonzero(days == 10*YEAR_DAYS)[0])),
                        ("end", len(days)-1)]:
        update(index)
        fig.savefig(output/f"{stem}_{name}.png", dpi=120)
    plt.close(fig)
    return {"file": movie_path.name, "frames": len(days), "fps": fps,
            "duration_seconds": len(days)/fps, "size_pixels": [1200, 760],
            "integrated_until_day": float(integration_days[-1]),
            "frame_equilibria": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    import imageio_ffmpeg
    matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    args.output.mkdir(parents=True, exist_ok=False)
    source_provenance = json.loads((args.input/"provenance.json").read_text())
    for name, expected in source_provenance["output_sha256"].items():
        if hashlib.sha256((args.input/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Diagnostic artifact hash mismatch: {name}")
    for name, expected in source_provenance["source_sha256"].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Diagnostic source hash mismatch: {name}")
    model = RationalSPAC(model_family()["E4"])
    days = frame_days()
    results = {
        "scope": "Rational E4; same physical frame times, axes and initial condition",
        "physical_window_years": [0., EARLY_YEARS],
        "playback": "Nonuniform time intervals; identical frame clock in both movies",
        "state_interpolation": "Linear water coordinates, log-linear carbon; display only",
        "zero_equilibria": {
            "conducting": "(psi_l, psi_l, 0, 0, 0, 0), -3 < psi_l <= 0",
            "nonconducting": "(psi_x, psi_l, 0, 0, 0, 0), psi_l <= -3, psi_x <= 0",
            "projection": "All zero-carbon equilibria in -4.1 <= psi_l <= 0 shown as a line",
            "red": "Positive transverse leaf-carbon eigenvalue for psi_l > -20/9",
            "grey": "Non-hyperbolic continuum; no claim of asymptotic stability of a point",
        },
        "movies": {},
    }
    for label, ramp in [("fast_10day_ramp", 10.), ("slow_10year_ramp", 10*YEAR_DAYS)]:
        with np.load(args.input/f"{label}.npz") as archive:
            data = {key: archive[key] for key in archive.files}
        eq = frozen_equilibria(model, days, ramp)
        results["movies"][label] = build_movie(
            model, data, eq, days, ramp,
            "10-day ramp" if ramp == 10. else "10-year ramp", args.output, 15)
    sources = [Path(__file__), ROOT/"pyproject.toml"]
    results["source_sha256"] = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
    }
    results["input_provenance_sha256"] = hashlib.sha256(
        (args.input/"provenance.json").read_bytes()).hexdigest()
    results["created_at_utc"] = datetime.now(timezone.utc).isoformat()
    (args.output/"results.json").write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")
    html = """<!doctype html>
<html lang="en"><meta charset="utf-8"><title>E4 VPD trajectory movies</title>
<style>body{font:16px sans-serif;max-width:1500px;margin:2em auto;padding:0 1em}
.movies{display:flex;gap:1em}section{width:50%}video{width:100%}button{margin:.5em}</style>
<h1>E4: VPD, leaf water potential, leaf carbon, and frozen equilibria</h1>
<p>Both movies use the same nonuniform physical clock and linear time axes (0–26 years).
Blue: stable positive equilibrium; red: unstable equilibrium; green: actual path.
Grey: non-hyperbolic zero-carbon continuum, not an isolated stable point.</p>
<button id="play">Play both</button><button id="pause">Pause both</button>
<button id="reset">Restart both</button><div class="movies">
<section><h2>10-day ramp</h2><video controls preload="metadata"
poster="e4_fast_10day_start.png" src="e4_fast_10day.mp4"></video></section>
<section><h2>10-year ramp</h2><video controls preload="metadata"
poster="e4_slow_10year_start.png" src="e4_slow_10year.mp4"></video></section></div>
<p>The fast path stops at its numerical leaf boundary; it is not extrapolated.
The slow movie shows 26 years, not the full 400-year post-pulse recovery.</p>
<script>
const videos=[...document.querySelectorAll('video')];
document.querySelector('#play').onclick=()=>Promise.all(videos.map(v=>v.play())).catch(e=>alert(e));
document.querySelector('#pause').onclick=()=>videos.forEach(v=>v.pause());
document.querySelector('#reset').onclick=()=>videos.forEach(v=>{v.pause();v.currentTime=0;});
</script></html>
"""
    (args.output/"index.html").write_text(html, encoding="utf-8")
    outputs = sorted(p for p in args.output.iterdir() if p.is_file())
    provenance = {
        "source_sha256": results["source_sha256"],
        "input_provenance_sha256": results["input_provenance_sha256"],
        "output_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},
    }
    (args.output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps({key: {k: v for k, v in value.items() if k != "frame_equilibria"}
                      for key, value in results["movies"].items()}, indent=2))


if __name__ == "__main__":
    main()
