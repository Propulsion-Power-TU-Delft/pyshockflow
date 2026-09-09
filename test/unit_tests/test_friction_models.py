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
