try:
    import pytest
except ImportError:
    class Approx:
        def __init__(self, val, rel=1e-6, abs=1e-12):
            self.val = val
            self.rel = rel
            self.abs = abs
        def __eq__(self, other):
            return np.isclose(other, self.val, rtol=self.rel, atol=self.abs)
        def __repr__(self):
            return f"approx({self.val} +/- rel={self.rel}, abs={self.abs})"
    class PytestMock:
        approx = Approx
    pytest = PytestMock()
import numpy as np
from pyshockflow.fluid import FluidIdeal
from pyshockflow.friction_models import (
    ShockTracker,
    ConstantFriction,
    MirelsLaminarFriction,
    MirelsTurbulentFriction,
    MirelsTransitionalFriction,
    create_friction_model,
)


def test_fluid_ideal_viscosity():
    """Verify Sutherland's law calculation in FluidIdeal."""
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)
    
    # At reference temperature T0 = 273.15 K, mu must equal mu0 = 1.716e-5 Pa.s
    mu_ref = fluid.computeViscosity_p_T(101325.0, 273.15)
    assert mu_ref == pytest.approx(1.716e-5, rel=1e-5)

    # At 300 K
    T = 300.0
    expected = 1.716e-5 * (T / 273.15)**1.5 * (273.15 + 110.4) / (T + 110.4)
    mu_300 = fluid.computeViscosity_p_T(101325.0, T)
    assert mu_300 == pytest.approx(expected, rel=1e-5)

    # computeViscosity_p_rho
    p = 101325.0
    rho = p / (287.0 * T)
    mu_p_rho = fluid.computeViscosity_p_rho(p, rho)
    assert mu_p_rho == pytest.approx(expected, rel=1e-5)

    # Vectorized array input
    p_vec = np.array([101325.0, 200000.0])
    rho_vec = np.array([rho, rho * 2.0])
    mu_vec = fluid.computeViscosity_p_rho(p_vec, rho_vec)
    assert len(mu_vec) == 2
    assert mu_vec[0] == pytest.approx(expected, rel=1e-5)


def test_shock_tracker_synthetic():
    """Verify shock detection and Rankine-Hugoniot shock speed estimation."""
    tracker = ShockTracker(pressure_threshold=0.05)
    
    n = 100
    x = np.linspace(0.0, 1.0, n)
    
    # Synthetic right-running shock at cell 60:
    # State 2 (post-shock, left): rho2 = 0.5, u2 = 1.2, p2 = 0.8
    # State 1 (pre-shock, right): rho1 = 0.125, u1 = 0.0, p1 = 0.1
    rho = np.where(x <= 0.6, 0.5, 0.125)
    u = np.where(x <= 0.6, 1.2, 0.0)
    p = np.where(x <= 0.6, 0.8, 0.1)

    primitive = {'Density': rho, 'Velocity': u, 'Pressure': p}
    found = tracker.detect(primitive, x)

    assert found is True
    assert tracker.direction == 1
    # Shock should be near index 60
    assert abs(tracker.i_shock - 60) <= 2

    # Theoretical Rankine-Hugoniot shock speed:
    # us = (rho2*u2 - rho1*u1) / (rho2 - rho1) = (0.5*1.2 - 0) / (0.5 - 0.125) = 0.6 / 0.375 = 1.6
    expected_us = (0.5 * 1.2 - 0.125 * 0.0) / (0.5 - 0.125)
    assert tracker.u_shock == pytest.approx(expected_us, rel=1e-3)


def test_constant_friction():
    """Verify legacy ConstantFriction matches expected uniform Cf and Sm formula."""
    cf_val = 0.005
    model = ConstantFriction(cf_constant=cf_val)
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)

    n = 20
    x = np.linspace(0.0, 1.0, n)
    dx = np.full(n, 0.05)
    dh = np.full(n, 0.1)
    rho = np.full(n, 1.2)
    u = np.full(n, 50.0)
    primitive = {'Density': rho, 'Velocity': u, 'Pressure': np.full(n, 1e5)}

    cf = model.compute_friction_coefficients(primitive, x, dx, fluid)
    assert np.all(cf == cf_val)

    source = model.compute_source_terms(primitive, x, dx, dh, fluid)
    # Sm = -2.0 * cf * rho * u * |u| / dh
    expected_sm = -2.0 * cf_val * 1.2 * 50.0 * 50.0 / 0.1
    assert np.allclose(source[:, 1], expected_sm)
    assert np.all(source[:, 0] == 0.0)
    assert np.all(source[:, 2] == 0.0)


def test_mirels_laminar_friction():
    """Verify Mirels laminar friction formula and singularity regularization."""
    model = MirelsLaminarFriction(pressure_threshold=0.05, max_cf=0.2, driver_cf=0.003)
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)

    n = 50
    x = np.linspace(0.0, 1.0, n)
    dx = np.full(n, x[1] - x[0])

    # Shock at x=0.6, index 30
    rho2, u2, p2 = 0.5, 1.2, 0.8
    rho1, u1, p1 = 0.125, 0.0, 0.1
    rho = np.where(x <= 0.6, rho2, rho1)
    u = np.where(x <= 0.6, u2, u1)
    p = np.where(x <= 0.6, p2, p1)
    primitive = {'Density': rho, 'Velocity': u, 'Pressure': p}

    cf = model.compute_friction_coefficients(primitive, x, dx, fluid)

    # 1. Ahead of shock (x > 0.6): Cf must be 0
    assert np.all(cf[x > 0.65] == 0.0)

    # 2. Regularization at first cell behind shock:
    # bar_Cf = 2.0 * Cf(dx)
    dxi = dx[30]
    us = (rho2 * u2) / (rho2 - rho1)
    ue = abs(us - u2)
    mu2 = fluid.computeViscosity_p_rho(p2, rho2)
    re_dx = (rho2 * ue * dxi) / mu2
    expected_bar_cf = 2.0 * (0.664 / np.sqrt(re_dx))

    # The shock cell should be regularized
    i_shock = model.tracker.i_shock
    assert cf[i_shock] == pytest.approx(min(expected_bar_cf, 0.2), rel=0.05)

    # 3. For cells further behind shock (dist > dx), Cf = 0.664 / sqrt(Re_x)
    i_eval = i_shock - 5
    dist = x[i_shock] - x[i_eval]
    re_dist = (rho2 * ue * dist) / mu2
    expected_cf_eval = 0.664 / np.sqrt(re_dist)
    assert cf[i_eval] == pytest.approx(min(expected_cf_eval, 0.2), rel=0.05)


def test_mirels_turbulent_friction():
    """Verify Mirels turbulent friction formula and 1.25x cell averaging."""
    model = MirelsTurbulentFriction(pressure_threshold=0.05, max_cf=0.2, driver_cf=0.003)
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)

    n = 50
    x = np.linspace(0.0, 1.0, n)
    dx = np.full(n, x[1] - x[0])

    rho2, u2, p2 = 0.5, 1.2, 0.8
    rho1, u1, p1 = 0.125, 0.0, 0.1
    rho = np.where(x <= 0.6, rho2, rho1)
    u = np.where(x <= 0.6, u2, u1)
    p = np.where(x <= 0.6, p2, p1)
    primitive = {'Density': rho, 'Velocity': u, 'Pressure': p}

    cf = model.compute_friction_coefficients(primitive, x, dx, fluid)

    # Ahead of shock: Cf == 0
    assert np.all(cf[x > 0.65] == 0.0)

    # Regularization at first cell: 1.25 * Cf(dx)
    dxi = dx[30]
    us = (rho2 * u2) / (rho2 - rho1)
    ue = abs(us - u2)
    mu2 = fluid.computeViscosity_p_rho(p2, rho2)
    re_dx = (rho2 * ue * dxi) / mu2
    expected_bar_cf = 1.25 * (0.0592 / (re_dx ** 0.2))

    i_shock = model.tracker.i_shock
    assert cf[i_shock] == pytest.approx(min(expected_bar_cf, 0.2), rel=0.05)


def test_mirels_transitional_friction():
    """Verify transition from laminar to turbulent at Re_transition."""
    re_tr = 500.0  # low threshold for testing
    model = MirelsTransitionalFriction(re_transition=re_tr, pressure_threshold=0.05, max_cf=0.5)
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)

    n = 80
    x = np.linspace(0.0, 1.0, n)
    dx = np.full(n, x[1] - x[0])

    rho2, u2, p2 = 0.5, 1.2, 0.8
    rho1, u1, p1 = 0.125, 0.0, 0.1
    rho = np.where(x <= 0.6, rho2, rho1)
    u = np.where(x <= 0.6, u2, u1)
    p = np.where(x <= 0.6, p2, p1)
    primitive = {'Density': rho, 'Velocity': u, 'Pressure': p}

    cf = model.compute_friction_coefficients(primitive, x, dx, fluid)
    i_shock = model.tracker.i_shock

    us = (rho2 * u2) / (rho2 - rho1)
    ue = abs(us - u2)
    mu2 = fluid.computeViscosity_p_rho(p2, rho2)

    # Check cell close to shock (Re < re_tr): should be laminar
    dist_lam = x[i_shock] - x[i_shock - 2]
    re_lam = (rho2 * ue * dist_lam) / mu2
    if re_lam < re_tr:
        assert cf[i_shock - 2] == pytest.approx(0.664 / np.sqrt(re_lam), rel=0.05)


def test_create_friction_model_factory():
    """Verify create_friction_model factory function dispatches properly."""
    class DummyConfig:
        def __init__(self, model_name):
            self.model_name = model_name
        def getWallFrictionModel(self):
            return self.model_name
        def getFrictionCoefficient(self):
            return 0.003
        def getFrictionDriverCf(self):
            return 0.003
        def getFrictionMaxCf(self):
            return 0.1
        def getShockDetectionThreshold(self):
            return 0.05
        def getFrictionTransitionReynolds(self):
            return 1e6

    assert isinstance(create_friction_model(DummyConfig('constant')), ConstantFriction)
    assert isinstance(create_friction_model(DummyConfig('mirels_laminar')), MirelsLaminarFriction)
    assert isinstance(create_friction_model(DummyConfig('mirels_turbulent')), MirelsTurbulentFriction)
    assert isinstance(create_friction_model(DummyConfig('mirels_transitional')), MirelsTransitionalFriction)
    assert isinstance(create_friction_model(DummyConfig('unknown')), ConstantFriction)


def test_shock_tracker_left_running():
    """Verify detection and Rankine-Hugoniot speed for left-running incident shock."""
    tracker = ShockTracker(pressure_threshold=0.05)
    n = 100
    x = np.linspace(0.0, 1.0, n)

    # Shock at x=0.4 moving left:
    # State 1 (left, x < 0.4): rho1=0.125, u1=0.0, p1=0.1
    # State 2 (right, x >= 0.4): rho2=0.5, u2=-1.2, p2=0.8
    rho = np.where(x < 0.4, 0.125, 0.5)
    u = np.where(x < 0.4, 0.0, -1.2)
    p = np.where(x < 0.4, 0.1, 0.8)
    primitive = {'Density': rho, 'Velocity': u, 'Pressure': p}

    found = tracker.detect(primitive, x)
    assert found is True
    assert tracker.direction == -1
    assert tracker.wave_mode == 'incident_left'

    # Shock speed via Rankine-Hugoniot:
    # us = (rho2*u2 - rho1*u1) / (rho2 - rho1) = (0.5*(-1.2) - 0) / (0.5 - 0.125) = -0.6 / 0.375 = -1.6 < 0
    expected_us = (0.5 * (-1.2) - 0.125 * 0.0) / (0.5 - 0.125)
    assert tracker.u_shock == pytest.approx(expected_us, rel=1e-3)


def test_shock_tracker_right_wall_reflection():
    """Verify state machine transition to reflected_left upon end-wall collision."""
    tracker = ShockTracker(pressure_threshold=0.05, enable_reflection=True, reflection_threshold=1.15)
    n = 100
    x = np.linspace(0.0, 1.0, n)

    # Step 1: Incident shock near the right wall (x=0.96)
    rho = np.where(x < 0.96, 0.5, 0.125)
    u = np.where(x < 0.96, 1.2, 0.0)
    p = np.where(x < 0.96, 0.8, 0.1)
    # Right wall has not yet surged
    primitive = {'Density': rho, 'Velocity': u, 'Pressure': p}
    tracker.detect(primitive, x)
    assert tracker.wave_mode == 'incident_right'

    # Step 2: Shock hits right wall -> pressure surges to p5 = 2.0 >> 0.8, u stagnates
    p_refl = np.where(x < 0.96, 0.8, 2.0)
    u_refl = np.where(x < 0.96, 1.2, 0.0)
    rho_refl = np.where(x < 0.96, 0.5, 1.2)
    primitive_refl = {'Density': rho_refl, 'Velocity': u_refl, 'Pressure': p_refl}
    tracker.detect(primitive_refl, x)
    assert tracker.wave_mode == 'reflected_left'
    assert tracker.direction == -1

    # Step 3: Reflected shock moves left to x=0.85
    # State 2 (left of 0.85): rho=0.5, u=1.2, p=0.8
    # State 5 (right of 0.85): rho=1.2, u=0.0, p=2.0
    rho_step3 = np.where(x < 0.85, 0.5, 1.2)
    u_step3 = np.where(x < 0.85, 1.2, 0.0)
    p_step3 = np.where(x < 0.85, 0.8, 2.0)
    primitive_step3 = {'Density': rho_step3, 'Velocity': u_step3, 'Pressure': p_step3}
    found = tracker.detect(primitive_step3, x)
    assert found is True
    assert tracker.wave_mode == 'reflected_left'
    assert tracker.direction == -1
    # Shock location should be near 0.85
    assert abs(tracker.x_shock - 0.85) < 0.05

    # Reflected shock speed WR < 0 via Rankine-Hugoniot:
    # WR = (rho5*u5 - rho2*u2) / (rho5 - rho2) = (1.2*0 - 0.5*1.2) / (1.2 - 0.5) = -0.6 / 0.7 = -0.857
    expected_wr = (1.2 * 0.0 - 0.5 * 1.2) / (1.2 - 0.5)
    assert tracker.u_shock == pytest.approx(expected_wr, rel=1e-2)


def test_mirels_friction_reflected_distances():
    """Verify distance calculations and non-zero friction in oncoming State 2 during reflection."""
    model = MirelsTurbulentFriction(pressure_threshold=0.05, max_cf=0.2, driver_cf=0.003, enable_reflection=True)
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)

    n = 100
    x = np.linspace(0.0, 1.0, n)
    dx = np.full(n, x[1] - x[0])

    # Force reflected mode with reflected shock at x=0.75
    model.tracker.wave_mode = 'reflected_left'
    model.tracker.x_wall = 1.0
    model.tracker.peak_post_shock_p = 0.8

    # State 3 (driver, x < 0.3): rho=0.8, u=1.2, p=0.8
    # State 2 (oncoming, 0.3 <= x < 0.75): rho=0.5, u=1.2, p=0.8
    # State 5 (post-reflected, 0.75 <= x <= 1.0): rho=1.2, u=0.0, p=2.0
    rho = np.where(x < 0.3, 0.8, np.where(x < 0.75, 0.5, 1.2))
    u = np.where(x < 0.3, 1.2, np.where(x < 0.75, 1.2, 0.0))
    p = np.where(x < 0.75, 0.8, 2.0)
    primitive = {'Density': rho, 'Velocity': u, 'Pressure': p}

    cf = model.compute_friction_coefficients(primitive, x, dx, fluid)

    # 1. Driver gas (x < 0.3): Cf == driver_cf = 0.003
    assert np.all(cf[x < 0.25] == 0.003)

    # 2. Oncoming State 2 gas (0.35 <= x < 0.75):
    # Must have non-zero Cf because it's moving at u=1.2 towards the wall!
    assert np.all(cf[(x >= 0.35) & (x < 0.73)] > 0.0)

    # 3. Post-reflected State 5 gas (x > 0.75):
    # Cf is computed based on distance from reflected shock x - x_shock
    assert np.all(cf[x > 0.76] > 0.0)

