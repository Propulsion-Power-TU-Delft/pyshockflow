import pytest
import numpy as np
from pathlib import Path
from pyshockflow.fluid import FluidIdeal
from pyshockflow.tank_boundary import Tank0D


def test_tank_initial_state():
    """Verify thermodynamic initialization of the 0D tank."""
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)
    v_tank = 0.05
    p0 = 100000.0
    T0 = 300.0

    tank = Tank0D(volume=v_tank, p_initial=p0, T_initial=T0, fluid=fluid, location='right')

    expected_rho0 = p0 / (287.0 * T0)
    expected_m0 = expected_rho0 * v_tank
    expected_e0 = p0 / ((1.4 - 1.0) * expected_rho0)
    expected_U0 = expected_m0 * expected_e0

    assert tank.rho == pytest.approx(expected_rho0, rel=1e-5)
    assert tank.m == pytest.approx(expected_m0, rel=1e-5)
    assert tank.e == pytest.approx(expected_e0, rel=1e-5)
    assert tank.U == pytest.approx(expected_U0, rel=1e-5)
    assert tank.p == pytest.approx(p0, rel=1e-5)
    assert tank.T == pytest.approx(T0, rel=1e-5)


def test_tank_charging_analytical():
    """
    Verify adiabatic tank charging under constant mdot and htot matches
    the exact analytical solution: dp/dt = (gamma - 1) * mdot * htot / V.
    """
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)
    v_tank = 0.1
    p0 = 100000.0
    T0 = 300.0

    tank = Tank0D(volume=v_tank, p_initial=p0, T_initial=T0, fluid=fluid, location='right')

    mdot = 0.02  # kg/s
    cp = 1.4 * 287.0 / (1.4 - 1.0)
    Tt = 350.0   # Total temperature of charging gas
    htot = cp * Tt

    tank.current_mdot = mdot
    tank.current_htot = htot

    dt = 0.001
    n_steps = 100
    t_total = dt * n_steps

    for step in range(1, n_steps + 1):
        tank.step(dt, step * dt)

    # Analytical pressure for adiabatic tank:
    # U(t) = U0 + mdot * htot * t
    # p(t) = (gamma - 1) * U(t) / V = p0 + (gamma - 1) * mdot * htot * t / V
    dp_analytical = (1.4 - 1.0) * mdot * htot * t_total / v_tank
    expected_p = p0 + dp_analytical

    assert tank.p == pytest.approx(expected_p, rel=1e-4)
    # Mass must equal m0 + mdot * t
    expected_m = tank.m0 + mdot * t_total
    assert tank.m == pytest.approx(expected_m, rel=1e-5)


def test_tank_subsonic_boundary_coupling():
    """Verify subsonic outflow imposes tank static pressure onto halo cell."""
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)
    v_tank = 0.1
    p_tank = 120000.0
    T_tank = 300.0

    tank = Tank0D(volume=v_tank, p_initial=p_tank, T_initial=T_tank, fluid=fluid, location='right')

    # Subsonic outflow at right boundary: u > 0, M < 1
    # Internal cell is index -2, halo is index -1
    primitive = {
        'Density': np.array([1.2, 1.2, 1.2, 1.2]),
        'Velocity': np.array([50.0, 50.0, 50.0, 50.0]),
        'Pressure': np.array([150000.0, 150000.0, 150000.0, 150000.0]),
        'Energy': np.array([3e5, 3e5, 3e5, 3e5]),
    }
    area_tube = 0.01

    tank.apply_boundary_condition(primitive, iInternal=-2, iHalo=-1, area_tube=area_tube, fluid=fluid)

    # Halo static pressure must match tank pressure
    assert primitive['Pressure'][-1] == pytest.approx(p_tank, rel=1e-5)
    # Velocity and density extrapolated from internal cell
    assert primitive['Velocity'][-1] == pytest.approx(50.0, rel=1e-5)
    assert primitive['Density'][-1] == pytest.approx(1.2, rel=1e-5)
    # Check that mdot was recorded positively (flow into tank)
    assert tank.current_mdot == pytest.approx(1.2 * 50.0 * 0.01, rel=1e-5)


def test_tank_supersonic_boundary_coupling():
    """Verify supersonic outflow imposes transparent boundary on halo cell."""
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)
    v_tank = 0.1
    p_tank = 100000.0
    T_tank = 300.0

    tank = Tank0D(volume=v_tank, p_initial=p_tank, T_initial=T_tank, fluid=fluid, location='right')

    # Supersonic outflow: a = sqrt(1.4*287*300) = 347 m/s, u = 600 m/s => M > 1
    primitive = {
        'Density': np.array([1.2, 1.2, 1.2, 1.2]),
        'Velocity': np.array([600.0, 600.0, 600.0, 600.0]),
        'Pressure': np.array([200000.0, 200000.0, 200000.0, 200000.0]),
        'Energy': np.array([4e5, 4e5, 4e5, 4e5]),
    }
    area_tube = 0.01

    tank.apply_boundary_condition(primitive, iInternal=-2, iHalo=-1, area_tube=area_tube, fluid=fluid)

    # Halo pressure must match internal cell pressure (transparent)
    assert primitive['Pressure'][-1] == pytest.approx(200000.0, rel=1e-5)
    assert primitive['Velocity'][-1] == pytest.approx(600.0, rel=1e-5)
    # Mass flow into tank is still positive
    assert tank.current_mdot == pytest.approx(1.2 * 600.0 * 0.01, rel=1e-5)


def test_tank_supersonic_adverse_backpressure_coupling():
    """Verify supersonic outflow communicates tank pressure when tank backpressure exceeds jet pressure."""
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)
    v_tank = 0.1
    p_tank = 500000.0  # Tank is at 5 bar (adverse backpressure)
    T_tank = 300.0

    tank = Tank0D(volume=v_tank, p_initial=p_tank, T_initial=T_tank, fluid=fluid, location='right')

    # Jet is supersonic: u = 600 m/s, p = 200000 Pa (2 bar < 5 bar tank pressure)
    primitive = {
        'Density': np.array([1.2, 1.2, 1.2, 1.2]),
        'Velocity': np.array([600.0, 600.0, 600.0, 600.0]),
        'Pressure': np.array([200000.0, 200000.0, 200000.0, 200000.0]),
        'Energy': np.array([4e5, 4e5, 4e5, 4e5]),
    }
    area_tube = 0.01

    tank.apply_boundary_condition(primitive, iInternal=-2, iHalo=-1, area_tube=area_tube, fluid=fluid)

    # Halo pressure must match the adverse tank pressure to allow the Riemann solver to generate unchoking shock
    assert primitive['Pressure'][-1] == pytest.approx(500000.0, rel=1e-5)
    assert primitive['Velocity'][-1] == pytest.approx(600.0, rel=1e-5)
    assert primitive['Density'][-1] == pytest.approx(1.2, rel=1e-5)
    assert tank.current_mdot == pytest.approx(1.2 * 600.0 * 0.01, rel=1e-5)


def test_tank_isothermal_mode():
    """Verify isothermal mode holds temperature strictly constant."""
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)
    v_tank = 0.1
    p0 = 100000.0
    T0 = 300.0

    tank = Tank0D(volume=v_tank, p_initial=p0, T_initial=T0, fluid=fluid, location='right', thermal_mode='isothermal')

    tank.current_mdot = 0.05
    tank.current_htot = 1e5
    dt = 0.01

    for step in range(50):
        tank.step(dt, (step + 1) * dt)

    assert tank.T == pytest.approx(T0, rel=1e-6)
    # Pressure must strictly equal rho * R * T0
    assert tank.p == pytest.approx(tank.rho * 287.0 * T0, rel=1e-5)


def test_tank_save_history(tmp_path):
    """Verify CSV export of tank history."""
    fluid = FluidIdeal(gmma=1.4, Rgas=287.0)
    tank = Tank0D(volume=0.1, p_initial=1e5, T_initial=300.0, fluid=fluid, location='right')
    tank.current_mdot = 0.01
    tank.current_htot = 3e5
    tank.step(0.01, 0.01)

    tank.save_history(tmp_path, 'test_history.csv')
    csv_file = tmp_path / 'test_history.csv'
    assert csv_file.exists()

    lines = csv_file.read_text().strip().split('\n')
    assert len(lines) >= 3  # Header + t=0 + t=0.01
    assert 'Time_s,Pressure_Pa,Temperature_K,Density_kgm3,Mass_kg,MassFlowRate_kgs' in lines[0]
