"""Analytically tractable deterministic feedback benchmarks, not VISIT laws.

The hypothetical carbon influx is P*x**2/(H**2+x**2), with P=m*g*H.
Loss is m*x. Vary g alone for a fold, or H alone at fixed g for moving
basins. Neither effective driver is calibrated to temperature or rainfall.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .nonlinear import ContinuousNonlinear, NonlinearTrajectory, simulate_continuous


FeedbackDriver = Literal["gain", "half_saturation"]


@dataclass(frozen=True)
class SaturatingFeedback:
    """A hypothetical one-pool carbon model with a state-feedback switch.

    Stocks and H use normalized carbon units; m has inverse model-time units.
    The baseline removes only the stock-dependent input multiplier, retaining
    P=m*g*H and the same loss. Its unique frozen equilibrium is g*H.
    """

    mortality: float = 1.0
    gain: float = 3.0
    half_saturation: float = 1.0

    def __post_init__(self) -> None:
        if (not np.all(np.isfinite([self.mortality, self.gain, self.half_saturation]))
                or self.mortality <= 0 or self.gain < 0 or self.half_saturation <= 0):
            raise ValueError("mortality and half_saturation must be positive; gain nonnegative")

    def parameters(self, driver: FeedbackDriver, value: float) -> tuple[float, float]:
        if driver not in ("gain", "half_saturation"):
            raise ValueError("driver must be gain or half_saturation")
        if not np.isfinite(value) or value < 0 or (driver == "half_saturation" and value == 0):
            raise ValueError("driver is outside its physical parameter domain")
        gain = value if driver == "gain" else self.gain
        half_saturation = value if driver == "half_saturation" else self.half_saturation
        return gain, half_saturation

    def system(self, driver: FeedbackDriver, *, feedback: bool = True) -> ContinuousNonlinear:
        self.parameters(driver, getattr(self, driver, np.nan))

        def rhs(state: np.ndarray, value: float) -> np.ndarray:
            gain, half_saturation = self.parameters(driver, value)
            production = self.mortality * gain * half_saturation
            if feedback:
                scaled = state / half_saturation
                production = production * scaled**2 / (1 + scaled**2)
            return production - self.mortality * state

        def jacobian(state: np.ndarray, value: float) -> np.ndarray:
            gain, half_saturation = self.parameters(driver, value)
            derivative = -self.mortality
            if feedback:
                scaled = state[0] / half_saturation
                derivative += 2 * self.mortality * gain * scaled / (1 + scaled**2)**2
            return np.array([[derivative]])

        return ContinuousNonlinear(
            rhs, jacobian, ("carbon",), ("normalized carbon",), driver,
            "1" if driver == "gain" else "normalized carbon", "model time",
            "model-agnostic derived analysis; hypothetical saturating feedback"
            if feedback else "model-agnostic derived analysis; affine baseline",
        )

    def equilibria(self, driver: FeedbackDriver, value: float, *, feedback: bool = True) -> np.ndarray:
        """All nonnegative frozen roots, ascending; the fold is returned once."""
        gain, half_saturation = self.parameters(driver, value)
        if not feedback:
            return np.array([gain * half_saturation])
        if gain < 2:
            return np.array([0.0])
        if gain == 2:
            return np.array([0.0, half_saturation])
        discriminant = np.sqrt(gain**2 - 4)
        return np.array([
            0.0, 2 * half_saturation / (gain + discriminant),
            half_saturation * (gain + discriminant) / 2,
        ])


@dataclass(frozen=True)
class SmoothRamp:
    """Finite C1 ramp with exact endpoint holds and duration 1/rate.

    rate is inverse time, not driver units/time. Maximum absolute driver
    velocity is 1.5*abs(end-start)*rate. This distinction matters in rate scans.
    """

    start: float
    end: float
    rate: float
    start_time: float = 0.0

    def __post_init__(self) -> None:
        if not np.all(np.isfinite([self.start, self.end, self.rate, self.start_time])) or self.rate <= 0:
            raise ValueError("ramp values must be finite and rate positive")

    @property
    def end_time(self) -> float:
        return self.start_time + 1 / self.rate

    @property
    def max_driver_speed(self) -> float:
        return 1.5 * abs(self.end - self.start) * self.rate

    def __call__(self, time: float) -> float:
        phase = float(np.clip((time - self.start_time) * self.rate, 0, 1))
        return self.start + (self.end - self.start) * phase**2 * (3 - 2 * phase)


@dataclass(frozen=True)
class FeedbackRampResult:
    model: SaturatingFeedback
    protocol: SmoothRamp
    ramp: NonlinearTrajectory
    hold: NonlinearTrajectory
    final_equilibria: np.ndarray
    final_basin: str
    distance_to_attractor: float
    final_residual: float


def run_feedback_ramp(
    model: SaturatingFeedback,
    driver: FeedbackDriver,
    protocol: SmoothRamp,
    *,
    hold_time: float = 40.0,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    steps_per_ramp: int = 100,
) -> FeedbackRampResult:
    """Start at the high frozen equilibrium, ramp, then hold the final driver.

    Classification uses the known scalar autonomous basin after forcing stops,
    not transient tracking error or a carbon-flux sign. Nonhyperbolic final
    forcing is left unresolved. No generic B/R label is inferred automatically.
    Increasing steps_per_ramp refines both the ramp and relaxation step caps.
    """
    if not np.isfinite(hold_time) or hold_time <= 0:
        raise ValueError("hold_time must be finite and positive")
    if not isinstance(steps_per_ramp, int) or steps_per_ramp < 2:
        raise ValueError("steps_per_ramp must be an integer >= 2")
    initial_roots = model.equilibria(driver, protocol.start)
    final_roots = model.equilibria(driver, protocol.end)
    if initial_roots.size != 3:
        raise ValueError("initial forcing must have a hyperbolic high equilibrium")
    system = model.system(driver)
    relaxation_step = 25 / (steps_per_ramp * model.mortality)
    ramp = simulate_continuous(
        system, [initial_roots[-1]],
        np.linspace(protocol.start_time, protocol.end_time, steps_per_ramp + 1), protocol,
        rtol=rtol, atol=atol,
        max_step=min(1 / (protocol.rate * steps_per_ramp), relaxation_step),
    )
    hold = simulate_continuous(
        system, ramp.states[-1],
        np.linspace(protocol.end_time, protocol.end_time + hold_time, 201), protocol,
        rtol=rtol, atol=atol, max_step=relaxation_step,
    )
    final_stock = hold.states[-1, 0]
    boundary_tolerance = 10 * (atol + rtol * max(1.0, final_roots[-1]))
    if final_stock < -boundary_tolerance:
        raise RuntimeError("integration left the nonnegative carbon domain")
    final_basin = "unresolved"
    target = np.nan
    if final_roots.size == 1:
        final_basin, target = "low", 0.0
    elif final_roots.size == 3:
        if final_stock < final_roots[1] - boundary_tolerance:
            final_basin, target = "low", 0.0
        elif final_stock > final_roots[1] + boundary_tolerance:
            final_basin, target = "high", final_roots[-1]
    return FeedbackRampResult(
        model, protocol, ramp, hold, final_roots, final_basin,
        float(abs(final_stock - target)),
        float(abs(system.evaluate([final_stock], protocol.end)[0])),
    )