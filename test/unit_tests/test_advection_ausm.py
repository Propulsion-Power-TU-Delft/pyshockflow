import pytest
import numpy as np
from pyshockflow.fluid import FluidIdeal, FluidReal
from pyshockflow.advection_ausm import AdvectionAUSMplusUP
from pyshockflow.ausm_vectorized import compute_ausm_flux_ideal, compute_ausm_flux_real
from pyshockflow.kernels_numba import (
    NUMBA_AVAILABLE,
    compute_ausm_flux_ideal_numba,
    compute_ausm_flux_real_numba
)
from pyshockflow import Driver, Config


@pytest.fixture
def ideal_fluid():
    return FluidIdeal(gmma=1.4, Rgas=287.05)


@pytest.fixture
def real_fluid():
    return FluidReal(fluid_name='CO2', fluid_library='CoolProp')


def test_ausm_ideal_consistency(ideal_fluid):
    """Test consistency between scalar, vectorized, and numba AUSM+-up for ideal gas."""
    rhoL = np.array([1.0, 1.0, 0.125, 0.5])
    rhoR = np.array([0.125, 1.0, 1.0, 0.25])
    uL = np.array([0.0, 100.0, -50.0, 0.0])
    uR = np.array([0.0, -100.0, 50.0, 10.0])
    pL = np.array([100000.0, 100000.0, 10000.0, 50000.0])
    pR = np.array([10000.0, 100000.0, 100000.0, 25000.0])

    n = len(rhoL)
    flux_scalar = np.zeros((n, 3))
    for i in range(n):
        ausm = AdvectionAUSMplusUP(rhoL[i], rhoR[i], uL[i], uR[i], pL[i], pR[i], ideal_fluid)
        flux_scalar[i, :] = ausm.computeFlux()

    flux_vec = compute_ausm_flux_ideal(rhoL, rhoR, uL, uR, pL, pR, ideal_fluid.gmma)
    np.testing.assert_allclose(flux_scalar, flux_vec, rtol=1e-8, atol=1e-6)

    if NUMBA_AVAILABLE:
        flux_numba = compute_ausm_flux_ideal_numba(rhoL, rhoR, uL, uR, pL, pR, ideal_fluid.gmma)
        np.testing.assert_allclose(flux_scalar, flux_numba, rtol=1e-8, atol=1e-6)


def test_ausm_real_consistency(real_fluid):
    """Test consistency between scalar, vectorized, and numba AUSM+-up for real fluid (CO2)."""
    rhoL = np.array([100.0, 150.0, 80.0])
    rhoR = np.array([50.0, 80.0, 120.0])
    uL = np.array([0.0, 20.0, -15.0])
    uR = np.array([0.0, -10.0, 25.0])
    pL = np.array([6.0e6, 7.5e6, 5.5e6])
    pR = np.array([3.0e6, 4.0e6, 6.5e6])

    n = len(rhoL)
    flux_scalar = np.zeros((n, 3))
    eL_arr = np.zeros(n)
    eR_arr = np.zeros(n)
    aL_arr = np.zeros(n)
    aR_arr = np.zeros(n)

    for i in range(n):
        ausm = AdvectionAUSMplusUP(rhoL[i], rhoR[i], uL[i], uR[i], pL[i], pR[i], real_fluid)
        flux_scalar[i, :] = ausm.computeFlux()
        eL_arr[i] = ausm.eL
        eR_arr[i] = ausm.eR
        aL_arr[i] = ausm.aL
        aR_arr[i] = ausm.aR

    flux_vec = compute_ausm_flux_real(rhoL, rhoR, uL, uR, pL, pR, eL_arr, eR_arr, aL_arr, aR_arr)
    np.testing.assert_allclose(flux_scalar, flux_vec, rtol=1e-8, atol=1e-6)

    if NUMBA_AVAILABLE:
        flux_numba = compute_ausm_flux_real_numba(rhoL, rhoR, uL, uR, pL, pR, eL_arr, eR_arr, aL_arr, aR_arr)
        np.testing.assert_allclose(flux_scalar, flux_numba, rtol=1e-8, atol=1e-6)


def test_driver_sod_ideal_ausm(tmp_path):
    """Run full Sod test simulation using AUSM+-up scheme."""
    ini_content = """[GEOMETRY]
LENGTH = 1.0
INTERFACE_LOCATION = 0.5

[FLUID]
FLUID_NAME = Air
FLUID_MODEL = ideal
FLUID_GAMMA = 1.4
GAS_R_CONSTANT = 287.05

[SIMULATION]
NUMBER_POINTS = 100
TIME_MAX = 0.0002
TIME_STEP_METHOD = cfl
CFL_MAX = 0.5
PRESSURE_LEFT = 100000.0
PRESSURE_RIGHT = 10000.0
DENSITY_LEFT = 1.0
DENSITY_RIGHT = 0.125
VELOCITY_LEFT = 0.0
VELOCITY_RIGHT = 0.0
TEMPERATURE_LEFT = 348.4
TEMPERATURE_RIGHT = 278.7
BOUNDARY_CONDITION_LEFT = transparent
BOUNDARY_CONDITION_RIGHT = transparent
NUMERICAL_SCHEME = ausm_plus_up
HIGH_ORDER = no

[OUTPUT]
FOLDER_NAME = output_ausm
FILE_NAME = sod_ausm.csv
SHOW_ANIMATION = no
FORMAT = csv
"""
    config_file = tmp_path / "input.ini"
    config_file.write_text(ini_content)

    config = Config(str(config_file))
    driver = Driver(config)
    driver.solve()

    rho = driver.solutionPrimitive['Density']
    p = driver.solutionPrimitive['Pressure']
    u = driver.solutionPrimitive['Velocity']

    assert np.all(rho > 0.0)
    assert np.all(p > 0.0)
    assert np.max(rho) <= 1.05
    assert np.min(rho) >= 0.12
    assert np.max(u) > 50.0


def test_driver_co2_real_ausm(tmp_path):
    """Run real gas CO2 shock tube simulation using AUSM+-up scheme."""
    ini_content = """[GEOMETRY]
LENGTH = 1.0
INTERFACE_LOCATION = 0.5

[FLUID]
FLUID_NAME = CO2
FLUID_MODEL = real
FLUID_LIBRARY = CoolProp

[SIMULATION]
NUMBER_POINTS = 50
TIME_MAX = 0.0001
TIME_STEP_METHOD = cfl
CFL_MAX = 0.5
PRESSURE_LEFT = 6000000.0
PRESSURE_RIGHT = 3000000.0
DENSITY_LEFT = 120.0
DENSITY_RIGHT = 55.0
VELOCITY_LEFT = 0.0
VELOCITY_RIGHT = 0.0
TEMPERATURE_LEFT = 300.0
TEMPERATURE_RIGHT = 300.0
BOUNDARY_CONDITION_LEFT = transparent
BOUNDARY_CONDITION_RIGHT = transparent
NUMERICAL_SCHEME = ausm_plus_up
HIGH_ORDER = no

[OUTPUT]
FOLDER_NAME = output_co2_ausm
FILE_NAME = co2_ausm.csv
SHOW_ANIMATION = no
FORMAT = csv
"""
    config_file = tmp_path / "input.ini"
    config_file.write_text(ini_content)

    config = Config(str(config_file))
    driver = Driver(config)
    driver.solve()

    rho = driver.solutionPrimitive['Density']
    p = driver.solutionPrimitive['Pressure']
    assert np.all(rho > 0.0)
    assert np.all(p > 0.0)
