import pytest
import numpy as np
from pyshockflow.fluid import FluidReal

@pytest.fixture
def real_fluid():
    return FluidReal(fluid_name='air', fluid_library='CoolProp')
    
def test_computeTemperature_p_rho(real_fluid):
    p = 101325
    rho = 1.225
    T = real_fluid.computeTemperature_p_rho(p, rho)
    assert T == pytest.approx(288.15, abs=0.5)

def test_computeInletQuantities(real_fluid):
    pressure = np.linspace(1E5, 1E7, 10)
    totPressure = pressure*1.5
    totTemperature = 288.15
    direction = 1
    
    for i in range(len(pressure)):
        rho, u, e = real_fluid.computeInletQuantities(pressure[i], totPressure[i], totTemperature, direction)
        s_total = real_fluid.computeEntropy_p_T(totPressure[i], totTemperature)
        s_static = real_fluid.computeEntropy_p_rho(pressure[i], rho)
        assert s_static == pytest.approx(s_total, abs=1e-3)

def test_computeInletFromRiemannInvariant_real(real_fluid):
    import CoolProp.CoolProp as CP
    totPressure = 101325.0
    totTemperature = 288.15
    s_t = CP.PropsSI('S', 'P', totPressure, 'T', totTemperature, 'Air')
    h_t = CP.PropsSI('H', 'P', totPressure, 'T', totTemperature, 'Air')

    # 1. Steady-state state at p = 90000 Pa
    p_exact = 90000.0
    rho_exact = CP.PropsSI('D', 'P', p_exact, 'S', s_t, 'Air')
    h_exact = CP.PropsSI('H', 'P', p_exact, 'S', s_t, 'Air')
    u_exact = float(np.sqrt(max(0.0, 2.0 * (h_t - h_exact))))

    # Compute inlet from Riemann invariant with exact internal cell
    rho_b, u_b, p_b, e_b = real_fluid.computeInletFromRiemannInvariant(
        u_int=u_exact, p_int=p_exact, rho_int=rho_exact,
        totPressure=totPressure, totTemperature=totTemperature, direction=1
    )
    assert p_b == pytest.approx(p_exact, rel=1e-3)
    assert u_b == pytest.approx(u_exact, rel=1e-3)
    assert rho_b == pytest.approx(rho_exact, rel=1e-3)

    # Verify isentropic total enthalpy and entropy conservation
    h_b = CP.PropsSI('H', 'P', p_b, 'D', rho_b, 'Air')
    s_b = CP.PropsSI('S', 'P', p_b, 'D', rho_b, 'Air')
    h_tot_b = h_b + 0.5 * u_b**2
    assert h_tot_b == pytest.approx(h_t, rel=1e-3)
    assert s_b == pytest.approx(s_t, rel=1e-3)

    # 2. Right inlet check (direction = -1)
    rho_b_r, u_b_r, p_b_r, e_b_r = real_fluid.computeInletFromRiemannInvariant(
        u_int=-u_exact, p_int=p_exact, rho_int=rho_exact,
        totPressure=totPressure, totTemperature=totTemperature, direction=-1
    )
    assert p_b_r == pytest.approx(p_exact, rel=1e-3)
    assert u_b_r == pytest.approx(-u_exact, rel=1e-3)