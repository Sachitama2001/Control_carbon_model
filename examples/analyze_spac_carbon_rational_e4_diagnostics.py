"""Generate focused E4 nullcline, capacity and bounded-VPD-path diagnostics."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from control_carbon.spac_carbon import YEAR_DAYS, carbon_matrix
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_rational import RationalSPAC
from prepare_manuscript_phase_portraits import model_family


ROOT = Path(__file__).resolve().parents[1]
POOL_NAMES = ("Leaf", "Wood", "Root", "Soil")
VPD_START = 0.00244
VPD_PEAK = 0.003
VPD_HOLD_DAYS = YEAR_DAYS
TAIL_DAYS = 400*YEAR_DAYS
LEAF_EXIT = 1e-6


def pulse_vpd(day, ramp_days):
    up_end = ramp_days
    hold_end = up_end+VPD_HOLD_DAYS
    down_end = hold_end+ramp_days
    if day < 0:
        return VPD_START
    if day < up_end:
        s = day/ramp_days
        return VPD_START+(VPD_PEAK-VPD_START)*s*s*(3-2*s)
    if day < hold_end:
        return VPD_PEAK
    if day < down_end:
        s = (day-hold_end)/ramp_days
        return VPD_PEAK+(VPD_START-VPD_PEAK)*s*s*(3-2*s)
    return VPD_START


def integrate_pulse(model, initial, ramp_years, scales):
    ramp_days = ramp_years*YEAR_DAYS
    pulse_end = 2*ramp_days+VPD_HOLD_DAYS
    end = pulse_end+TAIL_DAYS
    samples = np.unique(np.r_[
        np.linspace(0., ramp_days, 500),
        np.linspace(ramp_days, ramp_days+VPD_HOLD_DAYS, 400),
        np.linspace(ramp_days+VPD_HOLD_DAYS, pulse_end, 500),
        np.linspace(pulse_end, end, 3000),
    ])
    samples = samples[(samples >= 0) & (samples <= end)]
    def rhs(day, z):
        stocks = np.exp(z[2:])
        if not np.all(np.isfinite(stocks)) or np.any(stocks == 0):
            raise FloatingPointError("Log-carbon coordinates left floating-point range")
        state = np.r_[z[:2], stocks]
        values = model.rhs(state, pulse_vpd(day, ramp_days))
        return np.r_[values[:2], values[2:]/stocks]
    def jacobian(day, z):
        stocks = np.exp(z[2:])
        state = np.r_[z[:2], stocks]
        forcing = pulse_vpd(day, ramp_days)
        state_jacobian = model.jacobian(state, forcing)
        transform = np.diag(np.r_[1., 1., stocks])
        transformed = state_jacobian@transform
        transformed[2:] /= stocks[:, None]
        transformed[2:, 2:] -= np.diag(model.rhs(state, forcing)[2:]/stocks)
        return transformed
    def leaf_exit(day, z):
        return z[2]-np.log(LEAF_EXIT)
    leaf_exit.terminal = True
    leaf_exit.direction = -1
    run = solve_ivp(rhs, (0., end), np.r_[initial[:2], np.log(initial[2:])],
                    method="Radau", rtol=1e-9, atol=1e-11, t_eval=samples,
                    dense_output=True, jac=jacobian, events=leaf_exit)
    if not run.success:
        raise RuntimeError(f"E4 VPD-path integration failed: {run.message}")
    exit_day = float(run.t_events[0][0]) if run.status else None
    if exit_day is not None:
        days = np.unique(np.r_[run.t, exit_day])
        z = run.sol(days)
    else:
        days, z = run.t, run.y
    states = np.column_stack([z[:2].T, np.exp(z[2:].T)])
    forcing = np.array([pulse_vpd(t, ramp_days) for t in days])
    capacity = np.array([model.capacity(state) for state in states])
    potential = capacity-states[:, 2:]
    carbon_rhs = np.array([model.rhs(state, v)[2:] for state, v in zip(states, forcing)])
    identity = (carbon_rhs-(states[:, 2:]-capacity)@carbon_matrix(model.base.carbon).T)
    if np.max(np.abs(identity)) > 1e-12:
        raise RuntimeError("Original frozen-input capacity identity failed")
    return {
        "day": days, "state": states, "vpd_mpa": forcing,
        "capacity": capacity, "potential": potential,
        "status": "leaf_domain_exit" if exit_day is not None else "400_year_hold_completed",
        "exit_day": exit_day, "pulse_end_day": pulse_end,
        "nfev": int(run.nfev), "max_capacity_identity_residual": float(np.max(np.abs(identity))),
        "max_leaf_rate_abs": float(np.max(np.abs(carbon_rhs[:, 0]))),
    }


def spectral_crossings(model):
    leaves = np.unique(np.r_[np.geomspace(1e-7, .25, 400), np.linspace(.001, .25, 401)])
    abscissa = np.array([
        max(np.linalg.eigvals(model.jacobian(model.curve(leaf)[0], model.curve(leaf)[1])).real)
        for leaf in leaves
    ])
    found = []
    for left, right, a, b in zip(leaves[:-1], leaves[1:], abscissa[:-1], abscissa[1:]):
        if a*b >= 0:
            continue
        def leading_real(leaf):
            state, vpd, _ = model.curve(leaf)
            return float(max(np.linalg.eigvals(model.jacobian(state, vpd)).real))
        leaf = brentq(leading_real, left, right, xtol=1e-15)
        state, vpd, _ = model.curve(leaf)
        eigenvalues = np.linalg.eigvals(model.jacobian(state, vpd))
        neutral = eigenvalues[np.argmax(eigenvalues.real)]
        found.append({
            "leaf": float(leaf), "vpd_mpa": float(vpd), "state": state.tolist(),
            "eigenvalues_per_day": [[x.real, x.imag] for x in eigenvalues],
            "frequency_per_day": float(abs(neutral.imag)),
            "type": "fold" if abs(neutral.imag) < 1e-8 else "numerical complex-pair crossing",
        })
    return leaves, abscissa, found


def nullcline_slices(model, states):
    psi_leaf = np.linspace(-2.75, -.55, 240)
    leaf = np.linspace(.001, .16, 260)
    xx, ll = np.meshgrid(psi_leaf, leaf)
    records = {}
    for label, vpd in [("stable_2.44kPa", VPD_START), ("unstable_3.00kPa", VPD_PEAK)]:
        eq = states[label]
        psi_x, wood = eq[0], eq[3]
        state = np.tile(eq, (xx.size, 1))
        state[:, 1], state[:, 2] = xx.ravel(), ll.ravel()
        values = np.array([model.rhs(y, vpd) for y in state])
        psi_l_rate = values[:, 1].reshape(xx.shape)
        leaf_rate = values[:, 2].reshape(xx.shape)
        actual = model.rhs(eq, vpd)
        if max(abs(actual[1]), abs(actual[2])) > 1e-11:
            raise RuntimeError("Full equilibrium does not lie on the selected nullclines")
        records[label] = {
            "psi_leaf_mpa": psi_leaf, "leaf_carbon": leaf,
            "leaf_water_rate": psi_l_rate, "leaf_carbon_rate": leaf_rate,
            "psi_xylem_frozen_mpa": float(psi_x),
            "wood_carbon_frozen": float(wood),
            "equilibrium": eq.tolist(), "vpd_mpa": vpd,
            "slice": "psi_x, C_W, C_R and C_S frozen at the stated full equilibrium",
        }
    return records


def save_figure(fig, output, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(output/f"{stem}.{ext}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    base = model_family()["E4"]
    model = RationalSPAC(base)
    low = model.equilibria(VPD_START)[0]
    high = model.equilibria(VPD_PEAK)[0]
    states = {"stable_2.44kPa": low, "unstable_3.00kPa": high}
    scales = np.asarray(reference_scales(base).state)

    # The two displayed nullclines are slices through the six-dimensional field,
    # not a closed two-dimensional reduction.
    slices = nullcline_slices(model, states)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), constrained_layout=True)
    for ax, (name, data) in zip(axes, slices.items()):
        x, y = data["psi_leaf_mpa"], data["leaf_carbon"]
        xx, ll = np.meshgrid(x, y)
        a = data["leaf_water_rate"]
        b = data["leaf_carbon_rate"]*1000
        ax.contourf(xx, ll, b, levels=[-1e6, 0, 1e6],
                    colors=["#f2e2d5", "#dce9f0"], alpha=.8)
        water = ax.contour(xx, ll, a, levels=[0], colors=["#0072B2"], linewidths=2)
        carbon = ax.contour(xx, ll, b, levels=[0], colors=["#D55E00"], linewidths=2,
                            linestyles="--")
        eq = data["equilibrium"]
        ax.scatter(eq[1], eq[2], marker="*", s=140, color="black", zorder=5)
        ax.set(title=f"E4, VPD = {data['vpd_mpa']*1000:.2f} kPa",
               xlabel=r"Leaf water potential, $\psi_l$ (MPa)",
               ylabel=r"Leaf carbon, $C_L$ (kg C m$^{-2}$)")
        ax.text(.02, .02,
                rf"Fixed: $\psi_x={data['psi_xylem_frozen_mpa']:.3f}$ MPa, "
                rf"$C_W={data['wood_carbon_frozen']:.3f}$ kg C m$^{{-2}}$",
                transform=ax.transAxes, fontsize=8)
        ax.grid(alpha=.16)
        ax.legend([Line2D([], [], color="#0072B2", lw=2),
                   Line2D([], [], color="#D55E00", lw=2, ls="--"),
                   Line2D([], [], marker="*", color="black", ls="None", markersize=9),
                   Patch(facecolor="#f2e2d5"),
                   Patch(facecolor="#dce9f0")],
                  [r"$\dot{\psi}_l=0$", r"$\dot{C}_L=0$", "full equilibrium",
                   r"$\dot{C}_L<0$", r"$\dot{C}_L>0$"],
                  loc="best", fontsize=9)
    save_figure(fig, args.output, "e4_leaf_nullclines_slices")

    # A second slice shows both water nullclines in the (psi_x, psi_l) plane.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
    psi_x = np.linspace(-2.2, -.15, 260)
    psi_leaf = np.linspace(-3.8, -.55, 280)
    xx, ll = np.meshgrid(psi_x, psi_leaf)
    water_records = {}
    for ax, (name, eq) in zip(axes, states.items()):
        vpd = VPD_START if name == "stable_2.44kPa" else VPD_PEAK
        all_states = np.tile(eq, (xx.size, 1))
        all_states[:, 0], all_states[:, 1] = xx.ravel(), ll.ravel()
        water_rhs = np.array([model.rhs(y, vpd)[:2] for y in all_states])
        xrate, leafrate = water_rhs.T.reshape(2, *xx.shape)
        c1 = ax.contour(xx, ll, xrate, levels=[0], colors="#0072B2", linewidths=2)
        c2 = ax.contour(xx, ll, leafrate, levels=[0], colors="#D55E00",
                        linewidths=2, linestyles="--")
        ax.scatter(eq[0], eq[1], marker="*", color="black", s=140, zorder=5)
        ax.set(title=f"E4 water slice, VPD = {vpd*1000:.2f} kPa",
               xlabel=r"Xylem water potential, $\psi_x$ (MPa)",
               ylabel=r"Leaf water potential, $\psi_l$ (MPa)")
        ax.legend([Line2D([], [], color="#0072B2", lw=2),
                   Line2D([], [], color="#D55E00", lw=2, ls="--"),
                   Line2D([], [], marker="*", color="black", ls="None", markersize=9)],
                  [r"$\dot{\psi}_x=0$", r"$\dot{\psi}_l=0$", "full equilibrium"],
                  loc="best", fontsize=9)
        ax.grid(alpha=.16)
        water_records[name] = {"vpd_mpa": vpd, "equilibrium": eq.tolist(),
                               "frozen_carbon": eq[2:].tolist(),
                               "slice": "C_L,C_W,C_R,C_S frozen at the stated full equilibrium"}
    save_figure(fig, args.output, "e4_water_nullclines_slices")

    leaves, abscissa, crossings = spectral_crossings(model)
    vpd_curve = np.array([model.curve(x)[1] for x in leaves])
    stable = abscissa < 0
    fig, ax = plt.subplots(figsize=(8.2, 5.3), constrained_layout=True)
    ax.plot(vpd_curve[stable]*1000, leaves[stable], color="#0072B2", lw=2.3,
            label="Numerically stable equilibrium branch")
    ax.plot(vpd_curve[~stable]*1000, leaves[~stable], color="#D55E00", lw=2.3,
            label="Numerically unstable equilibrium branch")
    for item in crossings:
        ax.scatter(item["vpd_mpa"]*1000, item["leaf"], marker="D", color="black", zorder=4)
        ax.annotate(f"{item['vpd_mpa']*1000:.3f} kPa",
                    (item["vpd_mpa"]*1000, item["leaf"]),
                    xytext=(7, 7), textcoords="offset points", fontsize=9)
    ax.axvspan(VPD_START*1000, VPD_PEAK*1000, color="#999999", alpha=.16,
               label="VPD range of the pulse experiments")
    ax.set(xlabel="VPD (kPa)", ylabel=r"Leaf carbon at equilibrium, $C_L^*$ (kg C m$^{-2}$)",
           xlim=(0, 5.1), yscale="log")
    ax.grid(alpha=.18)
    ax.legend(fontsize=8)
    save_figure(fig, args.output, "e4_equilibrium_branch_and_stability")

    paths = {}
    for label, ramp_years in [("fast_10day_ramp", 10/YEAR_DAYS),
                              ("slow_10year_ramp", 10.)]:
        paths[label] = integrate_pulse(model, low, ramp_years, scales)

    # For each dynamic trajectory compute the frozen full equilibrium at that
    # instant by solving the monotone rational E4 equilibrium curve.
    for data in paths.values():
        ystar, capacity_residual = [], 0.
        for vpd in data["vpd_mpa"]:
            leaf = brentq(lambda c_leaf: model.curve(c_leaf)[1]-vpd,
                          1e-12, model.wet_leaf, xtol=1e-14)
            state = model.curve(leaf)[0]
            ystar.append(state)
            capacity_residual = max(capacity_residual,
                                    float(max(abs(model.capacity(state)-state[2:]))))
        data["frozen_equilibrium"] = np.asarray(ystar)
        data["frozen_equilibrium_capacity"] = data["frozen_equilibrium"][:, 2:].copy()
        data["max_frozen_capacity_residual"] = capacity_residual
        if capacity_residual > 1e-12:
            raise RuntimeError("Instantaneous frozen equilibrium does not satisfy K=C")

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    for col, (label, data) in enumerate(paths.items()):
        time_years = data["day"]/YEAR_DAYS
        ramp_name = "10-day ramp" if label.startswith("fast") else "10-year ramp"
        axes[0, col].plot(time_years, data["vpd_mpa"]*1000, color="#444444", lw=1.8)
        axes[0, col].axhline(3.167489, color=".55", ls=":",
                             label="25°C saturation VPD")
        axes[0, col].set(title=f"{ramp_name}; 1-year peak hold",
                         ylabel="VPD (kPa)")
        axes[1, col].plot(time_years, data["state"][:, 0],
                          label=r"Actual $\psi_x$", lw=1.8)
        axes[1, col].plot(time_years, data["frozen_equilibrium"][:, 0], ls="--",
                          label=r"Frozen-equilibrium $\psi_x^*$", lw=1.5)
        axes[1, col].plot(time_years, data["state"][:, 1], alpha=.75, lw=1.2,
                          label=r"Actual $\psi_l$")
        axes[1, col].set(xlabel="Time (years)", ylabel="Water potential (MPa)")
        early_window = data["pulse_end_day"]/YEAR_DAYS+5
        for row in range(2):
            axes[row, col].set_xlim(0, early_window)
            axes[row, col].grid(alpha=.18)
        axes[0, col].legend(fontsize=8)
        axes[1, col].legend(fontsize=8)
    save_figure(fig, args.output, "e4_vpd_pulses_and_xylem_trajectories")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8),
                             constrained_layout=True)
    for col, (label, data) in enumerate(paths.items()):
        time_years = data["day"]/YEAR_DAYS
        axes[0, col].plot(time_years, data["vpd_mpa"]*1000, color=".25", lw=1.4)
        axes[0, col].set(title=("10-day ramp" if col == 0 else "10-year ramp"),
                         ylabel="VPD (kPa)")
        axes[1, col].plot(data["state"][:, 1], data["state"][:, 2],
                          color="#0072B2", lw=1.7)
        points = axes[1, col].scatter(data["state"][:, 1], data["state"][:, 2],
                                      c=data["vpd_mpa"]*1000, s=5, cmap="viridis",
                                      norm=Normalize(VPD_START*1000, VPD_PEAK*1000))
        axes[1, col].scatter(low[1], low[2], marker="*", color="black", s=100,
                             label="Initial stable equilibrium")
        axes[1, col].set(xlabel=r"Leaf water potential, $\psi_l$ (MPa)",
                         ylabel=r"Leaf carbon, $C_L$ (kg C m$^{-2}$)")
        axes[1, col].legend(fontsize=8)
        if data["status"] == "leaf_domain_exit":
            axes[1, col].scatter(data["state"][-1, 1], data["state"][-1, 2],
                                 marker="v", color="black", s=40,
                                 label="Numerical leaf-domain exit")
        axes[0, col].grid(alpha=.18)
        axes[1, col].grid(alpha=.18)
        axes[0, col].set_xlim(0, data["pulse_end_day"]/YEAR_DAYS+5)
        axes[1, col].set_xlim(-4.1, -.5)
    fig.colorbar(points, ax=axes[1, :], label="VPD (kPa)", shrink=.86)
    save_figure(fig, args.output, "e4_leaf_carbon_phase_trajectories")

    fig, axes = plt.subplots(4, 2, figsize=(12, 12),
                             constrained_layout=True)
    for col, (label, data) in enumerate(paths.items()):
        time_years = data["day"]/YEAR_DAYS
        for pool, name in enumerate(POOL_NAMES):
            ax = axes[pool, col]
            ax.plot(time_years, data["state"][:, pool+2], color="#0072B2",
                    label=rf"$C_{{{name[0]}}}$")
            ax.plot(time_years, data["capacity"][:, pool], color="#D55E00", ls="--",
                    label=rf"$K_{{{name[0]}}}$")
            ax.set_yscale("symlog", linthresh=1e-10)
            ax.set_ylabel(rf"{name} pool ($\mathrm{{kg\,C\,m^{{-2}}}}$)")
            ax.set_xlim(0, data["pulse_end_day"]/YEAR_DAYS+5)
            ax.grid(alpha=.18)
            ax.legend(fontsize=8)
            if pool == 0:
                ax.set_title("10-day ramp" if col == 0 else "10-year ramp")
            if pool == len(POOL_NAMES)-1:
                ax.set_xlabel("Time (years)")
    save_figure(fig, args.output, "e4_classical_capacity_and_carbon_pools")

    fig, axes = plt.subplots(4, 2, figsize=(10, 12), constrained_layout=True)
    for col, (label, data) in enumerate(paths.items()):
        time_years = data["day"]/YEAR_DAYS
        for pool, name in enumerate(POOL_NAMES):
            ax = axes[pool, col]
            ax.plot(time_years, data["state"][:, pool+2], color="#0072B2",
                    label=rf"$C_{{{name[0]}}}$")
            ax.plot(time_years, data["capacity"][:, pool], color="#D55E00", ls="--",
                    label=rf"$K_{{{name[0]}}}$")
            ax.set_yscale("log")
            ax.set_ylabel(rf"{name} ($\mathrm{{kg\,C\,m^{{-2}}}}$)")
            ax.grid(alpha=.18)
            ax.legend(fontsize=8)
            if pool == 0:
                ax.set_title("10-day ramp" if col == 0 else "10-year ramp")
            if pool == len(POOL_NAMES)-1:
                ax.set_xlabel("Time (years)")
    save_figure(fig, args.output, "e4_capacity_and_stocks_full_horizon")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
    for ax, (label, data) in zip(axes, paths.items()):
        scale = scales[2:]
        time_years = data["day"]/YEAR_DAYS
        for pool, name in enumerate(POOL_NAMES):
            ax.plot(time_years, data["potential"][:, pool]/scale[pool],
                    label=rf"$(K_{{{name[0]}}}-C_{{{name[0]}}})/s_{{{name[0]}}}$")
        ax.axhline(0, color=".4", lw=.8)
        ax.set(title="10-day ramp" if label.startswith("fast") else "10-year ramp",
               xlabel="Time (years)", ylabel="Scaled capacity-minus-stock")
        ax.grid(alpha=.18)
        ax.legend(fontsize=8)
        ax.set_xlim(0, 25)
    save_figure(fig, args.output, "e4_capacity_minus_stock")

    fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    for col, (label, data) in enumerate(paths.items()):
        time = data["day"]/YEAR_DAYS
        points = None
        for pool, name in enumerate(POOL_NAMES):
            ax = axes[pool//2, pool%2]
            scale = scales[pool+2]
            ax.plot(data["state"][:, pool+2]/scale, data["capacity"][:, pool]/scale,
                    color=".72", lw=1)
            points = ax.scatter(data["state"][:, pool+2]/scale,
                                data["capacity"][:, pool]/scale,
                                c=time, s=5, cmap="viridis")
            upper = max(np.max(data["state"][:, pool+2]/scale),
                        np.max(data["capacity"][:, pool]/scale))*1.04
            ax.plot([0, upper], [0, upper], color="#D55E00", ls="--", lw=1.2)
            ax.set(title=f"{name} pool; " +
                   ("10-day ramp" if col == 0 else "10-year ramp"),
                   xlabel=rf"$C_{{{name[0]}}}/s_{{{name[0]}}}$",
                   ylabel=rf"$K_{{{name[0]}}}/s_{{{name[0]}}}$")
            ax.grid(alpha=.18)
    fig.colorbar(points, ax=axes, label="Time (years)", shrink=.84)
    save_figure(fig, args.output, "e4_capacity_vs_current_carbon")

    fig = plt.figure(figsize=(12, 5.5), constrained_layout=True)
    for index, (label, data) in enumerate(paths.items(), start=1):
        ax = fig.add_subplot(1, 2, index, projection="3d")
        trajectory = data["state"]
        colors = data["vpd_mpa"]*1000
        points = ax.scatter(trajectory[:, 1], trajectory[:, 2], trajectory[:, 3],
                            c=colors, s=5, cmap="viridis",
                            norm=Normalize(VPD_START*1000, VPD_PEAK*1000))
        ax.plot(trajectory[:, 1], trajectory[:, 2], trajectory[:, 3],
                color=".3", lw=.5, alpha=.6)
        ax.scatter(low[1], low[2], low[3], marker="*", s=100, color="black")
        ax.set(title="10-day ramp" if index == 1 else "10-year ramp",
               xlabel=r"$\psi_l$ (MPa)", ylabel=r"$C_L$ (kg C m$^{-2}$)",
               zlabel=r"$C_W$ (kg C m$^{-2}$)")
    fig.colorbar(points, ax=fig.axes, label="VPD (kPa)", shrink=.72)
    save_figure(fig, args.output, "e4_three_dimensional_phase_trajectories")

    records = {
        "scope": "Rational-canopy E4 only; classical frozen-input capacity; bounded smooth VPD pulse",
        "capacity_definition": "K=-M^{-1}u, u=(b_L,b_W,b_R,0)^T A; M constant",
        "units": {"time": "day; plots also show years", "water_potential": "MPa",
                  "carbon": "kg C m^-2", "VPD": "MPa; plots show kPa"},
        "forcing": {"initial_mpa": VPD_START, "peak_mpa": VPD_PEAK,
                    "hold_days": VPD_HOLD_DAYS, "smoothstep": "s^2(3-2s)",
                    "diagnostic_temperature_c": 25., "saturation_vpd_kpa": 3.167489},
        "initial_full_equilibrium": low.tolist(),
        "unstable_peak_full_equilibrium": high.tolist(),
        "initial_equilibrium_eigenvalues_per_day": [
            [x.real, x.imag] for x in np.linalg.eigvals(model.jacobian(low, VPD_START))],
        "peak_equilibrium_eigenvalues_per_day": [
            [x.real, x.imag] for x in np.linalg.eigvals(model.jacobian(high, VPD_PEAK))],
        "numerical_spectral_crossings": crossings,
        "nullcline_slices": {
            "leaf_plane": {name: {k: v for k, v in data.items()
                                  if k not in ["psi_leaf_mpa", "leaf_carbon",
                                               "leaf_water_rate", "leaf_carbon_rate"]}
                            for name, data in slices.items()},
            "water_plane": water_records,
        },
        "paths": {},
        "figure_interpretation": {
            "capacity": "K is computed from the original frozen-input definition at the instantaneous state; it is not a direct function of VPD alone",
            "nullclines": "Each nullcline figure is an explicitly frozen slice through the six-dimensional system, not a two-dimensional reduced model",
            "vpd_paths": "Smoothstep VPD pulse from 2.44 to 3.00 kPa and back, with a one-year peak hold",
            "trajectory_limit": "Leaf-domain exit at C_L=1e-6 is a numerical boundary, not mortality or proof of an absorbing zero state",
        },
    }
    for label, data in paths.items():
        arrays = {k: data[k] for k in ["day", "state", "vpd_mpa", "capacity",
                                       "potential", "frozen_equilibrium",
                                       "frozen_equilibrium_capacity"]}
        filename = f"{label}.npz"
        np.savez_compressed(args.output/filename, **arrays)
        records["paths"][label] = {
            "ramp_years": .02737850787132101 if label.startswith("fast") else 10.,
            "status": data["status"], "exit_day": data["exit_day"],
            "pulse_end_day": data["pulse_end_day"],
            "last_day": float(data["day"][-1]), "last_state": data["state"][-1].tolist(),
            "max_leaf_carbon": float(max(data["state"][:, 2])),
            "max_capacity_by_pool": data["capacity"].max(axis=0).tolist(),
            "max_abs_capacity_identity_residual": data["max_capacity_identity_residual"],
            "max_frozen_equilibrium_capacity_residual": data["max_frozen_capacity_residual"],
            "nfev": data["nfev"], "file": filename,
        }
        peak_index = int(np.argmin(abs(data["day"]-data["pulse_end_day"]/2)))
        pulse_end_index = int(np.argmin(abs(data["day"]-data["pulse_end_day"])))
        records["paths"][label]["state_at_peak_vpd"] = data["state"][peak_index].tolist()
        records["paths"][label]["capacity_at_peak_vpd"] = data["capacity"][peak_index].tolist()
        records["paths"][label]["state_at_pulse_end"] = data["state"][pulse_end_index].tolist()
        records["paths"][label]["capacity_at_pulse_end"] = data["capacity"][pulse_end_index].tolist()
        records["paths"][label]["final_capacity"] = data["capacity"][-1].tolist()
        records["paths"][label]["final_scaled_capacity_minus_stock"] = (
            data["potential"][-1]/scales[2:]).tolist()
        if data["max_capacity_identity_residual"] > 1e-12:
            raise RuntimeError("Capacity identity residual exceeded tolerance")
    if paths["fast_10day_ramp"]["status"] != "leaf_domain_exit":
        raise RuntimeError("Fast pulse did not reach the declared numerical leaf boundary")
    if paths["slow_10year_ramp"]["status"] != "400_year_hold_completed":
        raise RuntimeError("Slow pulse did not complete the declared 400-year hold")
    sources = [Path(__file__), ROOT/"examples/prepare_manuscript_phase_portraits.py",
               ROOT/"src/control_carbon/spac_carbon_rational.py",
               ROOT/"src/control_carbon/spac_carbon_capacity.py",
               ROOT/"src/control_carbon/spac_carbon.py",
               ROOT/"configs/spac_carbon_stage_b_synthetic.json",
               ROOT/"configs/spac_carbon_pressure_review.json"]
    records["source_sha256"] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in sources}
    records["python"] = {"numpy": np.__version__,
                         "scipy": __import__("scipy").__version__,
                         "matplotlib": matplotlib.__version__}
    records["created_at_utc"] = datetime.now(timezone.utc).isoformat()
    results_path = args.output/"results.json"
    results_path.write_text(json.dumps(records, indent=2, allow_nan=False)+"\n")
    outputs = sorted(path for path in args.output.iterdir()
                     if path.is_file() and path.name != "provenance.json")
    provenance = {
        "source_sha256": records["source_sha256"],
        "output_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in outputs},
        "scope": records["scope"],
    }
    (args.output/"provenance.json").write_text(
        json.dumps(provenance, indent=2, allow_nan=False)+"\n")
    print(json.dumps({
        "output": str(args.output),
        "crossings_kpa": [x["vpd_mpa"]*1000 for x in crossings],
        "paths": {k: {x: v[x] for x in ["status", "last_day", "last_state",
                                        "max_leaf_carbon", "max_abs_capacity_identity_residual"]}
                  for k, v in records["paths"].items()},
    }, indent=2))


if __name__ == "__main__":
    main()
