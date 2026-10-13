"""Three-dimensional E3 portraits from verified B1/B2 outputs; no new physiology."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.integrate import solve_ivp

from control_carbon.manzoni_hydraulics import HydraulicParameters
from control_carbon.spac_carbon import YEAR_DAYS, CarbonParameters
from control_carbon.spac_carbon_two_way import StructuralFeedback, equilibria_at_vpd
from control_carbon.spac_carbon_bifurcation import core_rhs, core_jacobian
from control_carbon.spac_carbon_mpa import PressureSPAC

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--b1", type=Path, default=ROOT/"artifacts/spac_carbon/stage_b1_final_20261009")
    parser.add_argument("--b2", type=Path, default=ROOT/"artifacts/spac_carbon/stage_b2_verified_20261012")
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon/phase3d_20261012")
    parser.add_argument("--atmospheric-pressure-mpa", type=float, default=.101325,
                        help="Explicit synthetic pressure assumption for archived mol/mol input conversion")
    args = parser.parse_args()
    b1 = json.loads((args.b1/"results.json").read_text())
    b2 = json.loads((args.b2/"results.json").read_text())
    if b1["config"] != b2["config"]:
        raise ValueError("B1/B2 input parameters must match")
    cfg = b2["config"]
    h = HydraulicParameters(**cfg["hydraulic"])
    k = cfg["carbon"]
    c = CarbonParameters(k["assimilation_max_per_year"]/YEAR_DAYS, tuple(k["allocation"]),
                         tuple(v/YEAR_DAYS for v in k["losses_per_year"]), tuple(k["litter_fractions"]),
                         k["specific_leaf_area"], k["carbon_dry_fraction"], k["extinction"])
    reference = b2["reference_stocks"]
    feedback = StructuralFeedback(reference[0], reference[1], cfg["leaf_floor"], cfg["minimum_leaf"],
                                  leaf_supply=1., leaf_demand=1.)
    root = cfg["root_potential"]
    pressure_model = PressureSPAC(h, c, feedback, root, args.atmospheric_pressure_mpa)
    local_summary = next(r for r in rows(args.b2/"periodic_summary.csv")
                         if r["variant"] == "E3" and float(r["amplitude"]) == .002)
    local_d = float(local_summary["D"])
    local_vpd_mpa = pressure_model.pressure(local_d)
    local_diagnostic = pressure_model.air_diagnostic(local_vpd_mpa, 25.)
    final_vpd_mpa = pressure_model.pressure(.03)
    equilibrium = equilibria_at_vpd(local_d, cfg["curve_bounds"], (), root, h, c, feedback)[0][:3]
    local = np.array([[float(r[key]) for key in ["xylem_MPa", "leaf_MPa", "L"]]
                      for r in rows(args.b2/"periodic_states.csv")
                      if r["variant"] == "E3" and float(r["amplitude"]) == .002])
    if len(local) < 100 or np.max(abs(local[-1]-local[0])) > 1e-8:
        raise RuntimeError("saved local cycle is missing or fails periodic closure")
    eigen, vectors = np.linalg.eig(core_jacobian(equilibrium, root, local_d, h, c, feedback))
    stable = vectors[:, np.argmin(eigen.real)].real
    near = []
    for sign in [-1, 1]:
        initial = equilibrium+sign*.001*stable
        result = solve_ivp(
            lambda t, y: core_rhs(y, root, local_d, h, c, feedback),
            (0., 5.), initial, method="Radau", rtol=1e-10, atol=1e-12,
            jac=lambda t, y: core_jacobian(y, root, local_d, h, c, feedback),
            t_eval=np.linspace(0., 5., 201),
        )
        if not result.success:
            raise RuntimeError(result.message)
        if np.linalg.norm(result.y[:, -1]-equilibrium) >= .01*np.linalg.norm(initial-equilibrium):
            raise RuntimeError("fast-direction perturbation did not contract over the plotted window")
        near.append(result.y.T)
    global_eq = equilibria_at_vpd(.03, cfg["curve_bounds"], (), root, h, c, feedback)[0][:3]
    trajectories = {}
    for duration in [10., 10*YEAR_DAYS]:
        data = [r for r in rows(args.b1/"trajectories.csv")
                if r["variant"] == "E3" and float(r["ramp_days"]) == duration]
        time = np.array([float(r["day"]) for r in data])
        state = np.array([[float(r[key]) for key in ["xylem_MPa", "leaf_MPa", "leaf_carbon"]] for r in data])
        if len(state) < 2:
            raise RuntimeError("saved B1 trajectory missing")
        trajectories[duration] = (time, state)
    args.output.mkdir(parents=True, exist_ok=False)
    scale = np.array([1e3, 1e3, 1e6])
    cycle_plot = (local-equilibrium)*scale

    def local_panel(ax, elevation=25, azimuth=-65):
        ax.plot(*cycle_plot.T, color="crimson", lw=2, label="Repelling periodic orbit")
        for number, index in enumerate([0, 100, 200, 300], start=1):
            ax.scatter(*cycle_plot[index], color="crimson", s=18)
            ax.text(*cycle_plot[index], f" {number}", color="crimson", fontsize=10)
        for i, path in enumerate(near):
            z = (path-equilibrium)*scale
            ax.plot(*z.T, color="seagreen", lw=1.8, label="Fast stable-direction transients (5 days)" if i == 0 else None)
            ax.scatter(*z[0], color="seagreen", marker="o", s=25)
        ax.scatter(0., 0., 0., color="black", marker="*", s=110, label="Stable equilibrium")
        ax.set(xlabel=r"$\Delta\psi_x$ (1e-3 MPa)", ylabel=r"$\Delta\psi_l$ (1e-3 MPa)",
               zlabel=r"$\Delta L$ (1e-6 kg C/m2)")
        ax.view_init(elev=elevation, azim=azimuth)
        ax.set_box_aspect((1, 1, 1))
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            axis.set_major_locator(MaxNLocator(5))

    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(121, projection="3d")
    for duration, color, label in [(10., "darkorange", "Fast ramp: 10 days"),
                                   (10*YEAR_DAYS, "royalblue", "Slow ramp: 10 years")]:
        time, state = trajectories[duration]
        ramp = time <= duration
        hold = time >= duration
        ax.plot(*state[ramp].T, ls="--", color=color, lw=1.7)
        ax.plot(*state[hold].T, color=color, lw=2.3, label=label)
        ax.scatter(*state[-1], color=color, marker="x" if duration == 10 else "o", s=65)
        ax.scatter(*state[0], color=color, marker="o", s=25)
    ax.scatter(*global_eq, color="black", marker="*", s=140, label=f"Stable equilibrium at VPD={final_vpd_mpa:.8f} MPa")
    ax.set(xlabel="xylem potential (MPa)", ylabel="leaf potential (MPa)", zlabel="Leaf carbon L (kg C/m2)",
           title=f"E3: ramp and final hold\nDashed: varying VPD; solid: frozen VPD={final_vpd_mpa:.8f} MPa")
    ax.view_init(elev=25, azim=-65)
    ax.set_box_aspect((1, 1, 1))
    ax.legend(fontsize=8, loc="upper left")
    ax = fig.add_subplot(122, projection="3d")
    local_panel(ax)
    ax.set_title(f"E3: FORMAL Hopf zoom, VPD={local_vpd_mpa:.8f} MPa\n"
                 f"Period T={float(local_summary['period_days']):.3f} days; red 1 -> 2 -> 3 -> 4 = time order")
    ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("Three active states: xylem potential, leaf potential, leaf carbon | SYNTHETIC\n"
                 f"Panels use DIFFERENT VPD | p_atm={args.atmospheric_pressure_mpa} MPa | orange endpoint is domain exit, NOT mortality")
    fig.subplots_adjust(left=.02, right=.95, bottom=.12, top=.82, wspace=.16)
    for suffix in ["png", "pdf"]:
        fig.savefig(args.output/f"phase3d.{suffix}", dpi=180)
    plt.close(fig)
    fig = plt.figure(figsize=(14, 11))
    for i, (elevation, azimuth) in enumerate([(25, -65), (25, 25), (60, -65), (10, -120)], start=1):
        ax = fig.add_subplot(2, 2, i, projection="3d")
        local_panel(ax, elevation, azimuth)
        ax.set_title(f"Same E3 orbit, view {i}")
    fig.suptitle(f"FORMAL E3 orbit | VPD={local_vpd_mpa:.8f} MPa | p_atm={args.atmospheric_pressure_mpa} MPa\n"
                 f"Red: repelling orbit; star: stable equilibrium | implied RH at 25 C = {local_diagnostic['implied_relative_humidity']:.3f}")
    fig.subplots_adjust(left=.03, right=.94, bottom=.06, top=.92, wspace=.14, hspace=.25)
    fig.savefig(args.output/"phase3d_views.png", dpi=160)
    plt.close(fig)
    with (args.output/"near_equilibrium.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["initial_sign", "day", "xylem_MPa", "leaf_MPa", "L"])
        for sign, path in zip([-1, 1], near):
            for t, y in zip(np.linspace(0., 5., 201), path):
                writer.writerow([sign, t, *y])
    inputs = [Path(__file__), args.b1/"results.json", args.b1/"trajectories.csv",
              args.b2/"results.json", args.b2/"periodic_summary.csv", args.b2/"periodic_states.csv",
              ROOT/"src/control_carbon/spac_carbon_bifurcation.py",
              ROOT/"src/control_carbon/spac_carbon_two_way.py",
              ROOT/"src/control_carbon/spac_carbon_bifurcation.py",
              ROOT/"src/control_carbon/spac_carbon_mpa.py",
              ROOT/"src/control_carbon/visitc_hydrology.py"]
    metadata = {"model": "unchanged synthetic E3", "archived_input_vpd_unit": "mol/mol",
                "vpd_output_unit": "MPa", "atmospheric_pressure_mpa": args.atmospheric_pressure_mpa,
                "global_final_vpd_mpa": final_vpd_mpa, "local_vpd_mpa": local_vpd_mpa,
                "local_air_diagnostic_at_25c": local_diagnostic,
                "local_equilibrium": equilibrium.tolist(),
                "local_cycle_radial_floquet": float(local_summary["radial_multiplier"]),
                "local_cycle_period_days": float(local_summary["period_days"]),
                "local_axes": "deviation from equilibrium; potentials scaled by 1e3, leaf carbon by 1e6",
                "near_trajectory_scope": "5-day fast-direction contraction, NOT global basin/no-R proof",
                "sha256": {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (args.output/"metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
