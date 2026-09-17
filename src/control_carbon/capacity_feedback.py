"""Capacity feedback examples for nonlinear matrix carbon systems.

These are model-agnostic counterexamples, not VISIT equations or forest fits.
They separate three objects that coincide only in the exogenous linear case:
the pointwise Luo-style capacity, a self-consistent frozen equilibrium, and an
attracting trajectory.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class CubicCapacityFeedback:
    """Positive scalar carbon example with two nonzero stable equilibria.

    ``dx/dt = relaxation * (capacity(x, driver) - x)`` is exactly the scalar
    matrix form with M=-relaxation and nonnegative input
    ``mu=relaxation*capacity`` on the documented domain [0, 12].  The cubic
    is chosen for transparent geometry, not calibrated ecology.
    """

    rate: float = 0.005
    relaxation: float = 1.0
    low: float = 1.0
    threshold_initial: float = 3.0
    threshold_shift: float = 4.0
    high_initial: float = 5.0
    high_shift: float = 5.0

    def __post_init__(self):
        values = (
            self.rate, self.relaxation, self.low, self.threshold_initial,
            self.high_initial,
        )
        if not all(np.isfinite(values)) or min(values) <= 0:
            raise ValueError("rates and initial positive branches must be finite and positive")
        if not all(np.isfinite((self.threshold_shift, self.high_shift))):
            raise ValueError("branch shifts must be finite")
        for driver in (0.0, 1.0):
            low, threshold, high = self.equilibria(driver)
            if not low < threshold < high:
                raise ValueError("branches must remain ordered on driver in [0, 1]")

    def equilibria(self, driver: float) -> tuple[float, float, float]:
        driver = float(driver)
        return (
            self.low,
            self.threshold_initial + self.threshold_shift * driver,
            self.high_initial + self.high_shift * driver,
        )

    def capacity(self, stock, driver):
        """Pointwise capacity map c_hat(x, lambda), in carbon-stock units."""
        stock = np.asarray(stock)
        low, threshold, high = self.equilibria(driver)
        cubic = (stock-low)*(stock-threshold)*(stock-high)
        return stock-(self.rate/self.relaxation)*cubic

    def capacity_derivative(self, stock: float, driver: float) -> float:
        low, threshold, high = self.equilibria(driver)
        derivative = ((stock-threshold)*(stock-high)
                      +(stock-low)*(stock-high)
                      +(stock-low)*(stock-threshold))
        return float(1-(self.rate/self.relaxation)*derivative)

    def rhs(self, stock, driver):
        return self.relaxation*(self.capacity(stock, driver)-np.asarray(stock))

    def jacobian(self, stock: float, driver: float) -> float:
        """Full derivative M*(1-D capacity), not frozen M alone."""
        return self.relaxation*(self.capacity_derivative(stock, driver)-1)

    def potential(self, stock, driver):
        return self.capacity(stock, driver)-np.asarray(stock)

    def potential_energy_derivative(self, stock: float, driver: float) -> float:
        """d[potential^2/2]/dt for frozen driver."""
        potential = float(self.potential(stock, driver))
        return self.relaxation*(self.capacity_derivative(stock, driver)-1)*potential**2


def linear_ramp(t: float, duration: float) -> float:
    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("duration must be finite and positive")
    return float(np.clip(t/duration, 0.0, 1.0))


def simulate_capacity_ramp(
    duration: float,
    model: CubicCapacityFeedback = CubicCapacityFeedback(),
    *,
    hold: float = 300.0,
    rtol: float = 1e-10,
    max_step: float = 0.25,
    output_step: float = 0.25,
) -> dict[str, np.ndarray | str | float]:
    """Ramp a single driver and hold its endpoint to identify the attractor."""
    if hold < 0 or output_step <= 0:
        raise ValueError("hold must be nonnegative and output_step positive")
    initial = model.equilibria(0.0)[2]
    end = duration+hold

    def fun(time, state):
        return [model.rhs(state[0], linear_ramp(time, duration))]

    solution = solve_ivp(
        fun, (0.0, end), [initial], method="DOP853", rtol=rtol,
        atol=rtol*0.01, max_step=max_step, dense_output=True,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    time = np.unique(np.r_[np.arange(0.0, end, output_step), end])
    stock = solution.sol(time)[0]
    driver = np.minimum(time/duration, 1.0)
    capacity = np.array([model.capacity(x, z) for x, z in zip(stock, driver)])
    final_low, final_threshold, final_high = model.equilibria(1.0)
    scale = min(final_threshold-final_low, final_high-final_threshold)
    if abs(stock[-1]-final_low) < 1e-6*scale:
        outcome = "low"
    elif abs(stock[-1]-final_high) < 1e-6*scale:
        outcome = "high"
    else:
        outcome = "unresolved"
    return {
        "time": time,
        "stock": stock,
        "driver": driver,
        "capacity": capacity,
        "potential": capacity-stock,
        "outcome": outcome,
        "final_threshold": final_threshold,
    }


def critical_ramp_duration(
    model: CubicCapacityFeedback = CubicCapacityFeedback(),
    *,
    lower: float = 20.0,
    upper: float = 120.0,
    iterations: int = 20,
) -> tuple[float, float]:
    """Numerical local bracket; no global uniqueness claim is implied."""
    if iterations <= 0:
        raise ValueError("iterations must be positive")
    low_outcome = simulate_capacity_ramp(lower, model)["outcome"]
    high_outcome = simulate_capacity_ramp(upper, model)["outcome"]
    if (low_outcome, high_outcome) != ("low", "high"):
        raise ValueError("duration endpoints must bracket low and high outcomes")
    for _ in range(iterations):
        middle = (lower+upper)/2
        outcome = simulate_capacity_ramp(middle, model)["outcome"]
        if outcome == "unresolved":
            raise RuntimeError("increase hold time near the basin boundary")
        if outcome == "low":
            lower = middle
        else:
            upper = middle
    return lower, upper


def periodic_capacity(time, mean=5.0, amplitude=2.0, omega=2*np.pi):
    return mean+amplitude*np.sin(omega*np.asarray(time))


def periodic_attractor(time, relaxation=1.0, mean=5.0, amplitude=2.0,
                       omega=2*np.pi):
    """Exact attracting periodic orbit of dx/dt=k(c(t)-x)."""
    time = np.asarray(time)
    gain = relaxation/np.sqrt(relaxation**2+omega**2)
    phase = np.arctan2(omega, relaxation)
    return mean+amplitude*gain*np.sin(omega*time-phase)
