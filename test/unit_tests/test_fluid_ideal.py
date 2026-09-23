import pytest
from pyshockflow.fluid import FluidIdeal

@pytest.fixture
def ideal_fluid():
    return FluidIdeal(gmma=1.0, Rgas=1.0)

def test_computeTemperature_p_rho(ideal_fluid):
    p = 100
    rho = 1
    T = ideal_fluid.computeTemperature_p_rho(p, rho)
    assert T == pytest.approx(100, abs=1e-6)

def test_computeInletFromRiemannInvariant_ideal():
    fluid = FluidIdeal(gmma=1.4, Rgas=287.05)
    totPressure = 101325.0
    totTemperature = 288.15

    # 1. Steady-state consistency check for left inlet (direction = +1)
    M = 0.35
    T_exact = totTemperature / (1.0 + 0.2 * M**2)
    p_exact = totPressure / (1.0 + 0.2 * M**2)**3.5
    rho_exact = p_exact / (287.05 * T_exact)
    a_exact = (1.4 * 287.05 * T_exact)**0.5
    u_exact = M * a_exact

    rho_b, u_b, p_b, e_b = fluid.computeInletFromRiemannInvariant(
        u_int=u_exact, p_int=p_exact, rho_int=rho_exact,
        totPressure=totPressure, totTemperature=totTemperature, direction=1
    )
    assert p_b == pytest.approx(p_exact, rel=1e-5)
    assert u_b == pytest.approx(u_exact, rel=1e-5)
    assert rho_b == pytest.approx(rho_exact, rel=1e-5)

    # 2. Steady-state consistency check for right inlet (direction = -1)
    rho_b_r, u_b_r, p_b_r, e_b_r = fluid.computeInletFromRiemannInvariant(
        u_int=-u_exact, p_int=p_exact, rho_int=rho_exact,
        totPressure=totPressure, totTemperature=totTemperature, direction=-1
    )
    assert p_b_r == pytest.approx(p_exact, rel=1e-5)
    assert u_b_r == pytest.approx(-u_exact, rel=1e-5)
    assert rho_b_r == pytest.approx(rho_exact, rel=1e-5)

    # 3. Acoustic wave invariant check
    # Along C- (u - a): J^- = u - 2*a/(gamma - 1)
    # Impose an acoustic perturbation from domain interior
    delta_u = -4.0
    delta_p = rho_exact * a_exact * delta_u  # du = dp / (rho*a)
    p_pert = p_exact + delta_p
    u_pert = u_exact + delta_u
    rho_pert = p_pert / (287.05 * (T_exact * (p_pert / p_exact)**(0.4 / 1.4)))

    rho_b_p, u_b_p, p_b_p, e_b_p = fluid.computeInletFromRiemannInvariant(
        u_int=u_pert, p_int=p_pert, rho_int=rho_pert,
        totPressure=totPressure, totTemperature=totTemperature, direction=1
    )
    a_b_p = (1.4 * p_b_p / rho_b_p)**0.5
    a_pert = (1.4 * p_pert / rho_pert)**0.5
    J_int = u_pert - 2.0 * a_pert / 0.4
    J_boundary = u_b_p - 2.0 * a_b_p / 0.4
    assert J_boundary == pytest.approx(J_int, rel=1e-4)
    