# pyshockflow

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python: 3.10+](https://img.shields.io/badge/python-3.10+-brightgreen.svg)](https://www.python.org/)
[![Acceleration: Numba JIT](https://img.shields.io/badge/acceleration-Numba%20JIT%20%2B%20SIMD-orange.svg)](https://numba.pydata.org/)
[![Thermodynamics: CoolProp](https://img.shields.io/badge/thermodynamics-CoolProp%20HEOS%20%2B%20LuT-blueviolet.svg)](http://www.coolprop.org/)
[![Tests: Pytest Passing](https://img.shields.io/badge/tests-pytest%20passing-success.svg)](test/regression_tests/)
[![Platform: Linux | macOS | Windows](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](https://github.com/Propulsion-Power-TU-Delft/pyshockflow)

**`pyshockflow`** is an object-oriented, high-performance finite-volume solver for quasi-1D compressible flows of ideal and real fluids. Developed in the **Propulsion & Power group** at **Delft University of Technology (TU Delft)**, the code is tailored for the study of unsteady wave dynamics, shock tubes, Ludwieg tubes, non-ideal compressible fluid dynamics (NICFD), non-classical gas dynamics (such as dense gases, supercritical fluids, and BZT vapors), variable-area nozzle flows, and coupled tank discharge/charging systems.

---

## Table of Contents

- [Key Features](#key-features)
- [Performance & Acceleration](#performance--acceleration)
- [Installation & Setup](#installation--setup)
- [Quick Start](#quick-start)
  - [1. Running via CLI](#1-running-via-cli)
  - [2. Running via Python API](#2-running-via-python-api)
  - [3. Post-Processing & Results](#3-post-processing--results)
- [Input Configuration Reference (`input.ini`)](#input-configuration-reference-inputini)
- [Numerical Methods & Physical Models](#numerical-methods--physical-models)
  - [Governing Equations](#governing-equations)
  - [Riemann Solvers & Numerical Fluxes](#riemann-solvers--numerical-fluxes)
  - [High-Order MUSCL & Slope Limiters](#high-order-muscl--slope-limiters)
  - [Thermodynamic Models & Look-Up Table (LuT)](#thermodynamic-models--look-up-table-lut)
  - [Wall Friction & Mirels Boundary Layer Correlations](#wall-friction--mirels-boundary-layer-correlations)
  - [Bidirectional & Reflected Shock Tracking](#bidirectional--reflected-shock-tracking)
  - [0D Lumped-Parameter Tank Boundary Condition](#0d-lumped-parameter-tank-boundary-condition)
  - [Quasi-1D Nozzles & Heat Transfer](#quasi-1d-nozzles--heat-transfer)
- [Verification & Automated Test Suite](#verification--automated-test-suite)
- [Repository Structure](#repository-structure)
- [Authors & Contact](#authors--contact)
- [Citation & References](#citation--references)
- [License](#license)

---

## Key Features

- **Extensive Riemann Solver Library**:
  - **Godunov**: Exact Riemann solver for ideal gas flows.
  - **Roe (Ideal)**: Classic Roe approximate Riemann solver with Harten–Hyman entropy fix.
  - **Roe-Arabi (Real Gas)**: Real-gas extension by Arabi et al. (2017) with speed-of-sound averaging and dissipation correction.
  - **Roe-Vinokur (Real Gas)**: Generalized algebraic projection formulation by Vinokur & Montagné (1990) employing fundamental thermodynamic derivatives $\chi$ and $\kappa$.
  - **HLLC**: Harten-Lax-van Leer Contact solver for both ideal and non-ideal fluids.
  - **AUSM+-up**: Liou (2006) Advection Upstream Splitting Method for all-speed regimes and shock-dominated real-gas flows.
- **High-Order Spatial Reconstruction**:
  - Second-order MUSCL reconstruction on cell-averaged primitives.
  - TVD slope limiters: **Van Albada**, **Van Leer**, **Minmod**, and **Superbee**.
- **Thermodynamic Flexibility & High-Accuracy EOS**:
  - Caloric and thermal **Ideal Gas** ($P = \rho R T$).
  - **Real Gas Models** backed by low-level C++ **CoolProp** Helmholtz energy equations of state (`AbstractState('HEOS')`).
  - **2D Bicubic Spline Look-Up Tables (LuT)** on $(\log_{10} P, T)$ ensuring thermodynamic consistency ($a^2 \equiv \chi + \kappa h > 0$) with fast vectorized Newton–Raphson state recovery.
- **Advanced Wall Friction & Shock Boundary Layer Modeling**:
  - Constant Darcy/Fanning wall friction factor.
  - **Mirels Laminar Boundary Layer** ($C_f = 0.664 / \sqrt{Re_x}$, NACA TN 3401).
  - **Mirels Turbulent Boundary Layer** ($C_f = 0.0592 / Re_x^{0.2}$, NACA TN 3712 / AIAA J. 1964).
  - **Mirels Transitional Model** switching dynamically at a specified Reynolds number $Re_{tr}$.
  - Leading-edge finite-volume singularity cell regularization ($\bar{C}_f = \frac{1}{1-n} C_f(\Delta x)$).
- **Intelligent Shock Wave & Reflection Tracking (`ShockTracker`)**:
  - Four-mode finite-state machine (`incident_right`, `reflected_left`, `incident_left`, `reflected_right`).
  - Automatic detection of incident shock, contact discontinuity, wall collision, and reflected wave front.
- **Dynamic 0D Tank Boundary Condition (`Tank0D`)**:
  - Coupled reservoir ODEs for mass and energy conservation: $\frac{dm}{dt} = \dot{m}$, $\frac{dU}{dt} = \dot{m} h_{tot}$.
  - Simulates dynamic charging and blowdown/discharge with transient backpressure (adiabatic or isothermal).
- **Quasi-1D Geometry & Thermal Source Terms**:
  - Variable cross-section area profiles (de Laval convergent-divergent nozzles, thrusters) via CSV input (`NOZZLE_FILEPATH`).
  - Geometrical source terms following Vimercati & Guardone (2018).
  - Prescribed wall heat flux source term ($q_w$ $[\text{W/m}^2]$).
- **Mesh & Simulation Controls**:
  - Uniform Cartesian mesh or localized adaptive refinement (`ADAPT_MESH_REFINEMENT`).
  - Adaptive CFL time-stepping or fixed $\Delta t$.
  - Unsteady wave propagation or steady-state convergence tracking with residual monitoring ($\|R\|$).
  - Full simulation restart capability (`RESTART_FILE`).
- **Publication-Ready Post-Processing**:
  - Matplotlib styling modules (`nicfd_styles.py`, `thesis_plots.py`) matching conference guidelines (e.g., NICFD 2026 LaTeX template with single/two-column dimensions).
  - Animation and video generators for transient wave propagation.

---

## Performance & Acceleration

`pyshockflow` incorporates a two-tier computational acceleration architecture:

1. **SIMD NumPy Vectorization**: Interface fluxes, wave decompositions, and slope limiters are evaluated across all cell interfaces simultaneously.
2. **Numba LLVM JIT Kernel Compilation**: Critical loops in `kernels_numba.py` are decorated with `@njit(fastmath=True)`, fusing flux evaluation, MUSCL reconstruction, and residual assembly into single-pass loops held entirely in hardware CPU registers/L1 cache with **zero heap allocations**.

| Configuration / Kernel | Baseline Python Loops | Vectorized NumPy | Numba JIT Fused Loops | **Speedup** |
| :--- | :---: | :---: | :---: | :---: |
| **Sod Shock Tube (1,000 cells, Ideal)** | $11.43\text{ s}$ | $0.0272\text{ s}$ | **$0.0078\text{ s}$** | **$1,465\times$** |
| **Toro Ideal Roe Flux (per step)** | $64.7\ \mu\text{s}$ | $64.7\ \mu\text{s}$ | **$10.5\ \mu\text{s}$** | **$6.2\times$** |
| **Arabi Real Gas Roe Flux (per step)** | $69.4\ \mu\text{s}$ | $69.4\ \mu\text{s}$ | **$9.5\ \mu\text{s}$** | **$7.3\times$** |
| **Vinokur Real Gas Roe Flux (per step)** | $122.2\ \mu\text{s}$ | $122.2\ \mu\text{s}$ | **$17.1\ \mu\text{s}$** | **$7.2\times$** |
| **Real Gas Helmholtz vs 2D LuT** | $23.86\text{ s}$ (Exact EOS) | — | **$6.40\text{ s}$ (LuT Spline)** | **$3.73\times$** |

*Note: If Numba is not installed, the solver automatically falls back to NumPy vectorized execution with zero code changes.*

---

## Installation & Setup

### 1. Clone the Repository
```bash
git clone https://github.com/Propulsion-Power-TU-Delft/pyshockflow.git
cd pyshockflow
```

### 2. Create and Activate Conda Environment
```bash
conda env create -f environment.yml
conda activate pyshockflow
```

Or manually with Miniforge / Conda:
```bash
conda create -n pyshockflow python=3.12 numpy scipy matplotlib coolprop pytest -c conda-forge
conda activate pyshockflow
```

*(Optional, recommended for max performance)* Install Numba:
```bash
pip install numba
```

### 3. Install `pyshockflow` in Editable Mode
```bash
pip install -e .
```

### 4. Verify Installation
Run the regression test suite to ensure all solvers and thermodynamics pass:
```bash
pytest test/regression_tests/test_cases.py
pytest test/unit_tests/
```

---

## Quick Start

### 1. Running via CLI

Every simulation case is defined by an `input.ini` file. Navigate to any testcase directory and run:

```bash
cd testcases/validation_ideal/roe/
python main.py
```

### 2. Running via Python API

You can programmatically configure and execute simulations in scripts or Jupyter notebooks:

```python
from pyshockflow import Config, Driver

# Load configuration file
config = Config('input.ini')

# Initialize driver and execute finite-volume time integration
tube = Driver(config)
tube.solve()
```

### 3. Post-Processing & Results

Results are saved into the directory specified by `FOLDER_NAME` in `input.ini`. During execution, raw solution snapshots are dumped, and upon completion they are assembled into `Results.pik` containing raw NumPy arrays:

```python
import pickle
import matplotlib.pyplot as plt
from pyshockflow.nicfd_styles import set_nicfd_style, create_figure

# Enable publication-quality formatting
set_nicfd_style()

# Load assembled solution
with open('Results/Results.pik', 'rb') as f:
    data = pickle.load(f)

x = data['X Coords']            # Spatial nodes [m]
t = data['Time']                # Time history [s]
p = data['Primitive']['Pressure']  # Shape: (nNodes, nTimes) [Pa]
rho = data['Primitive']['Density'] # Shape: (nNodes, nTimes) [kg/m^3]
u = data['Primitive']['Velocity']  # Shape: (nNodes, nTimes) [m/s]

# Plot final time step
fig, ax = create_figure(fraction=0.8, aspect_ratio=1.6)
ax.plot(x, p[:, -1] / 1e5, 'k-', lw=1.5, label=f'$t = {t[-1]*1e3:.2f}$ ms')
ax.set_xlabel('$x$ [m]')
ax.set_ylabel('Pressure [bar]')
ax.legend()
plt.show()
```

---

## Input Configuration Reference (`input.ini`)

The simulation settings are configured via an INI file containing four main sections:

```ini
[GEOMETRY]
; Total tube / duct length [m]
LENGTH = 1.0

; Initial diaphragm / interface position [m]
INTERFACE_LOCATION = 0.5

; Domain topology: 'default' (constant cross-section) or 'nozzle' (variable area)
TOPOLOGY = default

; Path to CSV file specifying nozzle area profile x, A(x) (only if TOPOLOGY = nozzle)
NOZZLE_FILEPATH = nozzle.csv

; Reference cross-sectional area [m^2]
REFERENCE_AREA = 1.0e-4


[SIMULATION]
; Number of grid cells in the physical domain
NUMBER_POINTS = 500

; Simulation time limit [s] and maximum CFL condition (< 1.0)
TIME_MAX = 0.002
CFL_MAX = 0.85
TIME_STEP_METHOD = adaptive    ; 'adaptive' or 'constant'

; Initial states (Left: high-pressure driver, Right: low-pressure driven)
PRESSURE_LEFT = 1.0e6          ; Static pressure [Pa]
PRESSURE_RIGHT = 1.0e5
TEMPERATURE_LEFT = 300.0       ; Static temperature [K]
TEMPERATURE_RIGHT = 300.0
VELOCITY_LEFT = 0.0            ; Initial flow velocity [m/s]
VELOCITY_RIGHT = 0.0

; Boundary conditions: 'reflective', 'transparent', 'periodic', 'inlet', 'outlet', 'tank'
BOUNDARY_CONDITION_LEFT = reflective
BOUNDARY_CONDITION_RIGHT = transparent

; Numerical Scheme: 'godunov', 'roe', 'roe_arabi', 'roe_vinokur', 'hllc', 'ausm+up'
NUMERICAL_SCHEME = roe

; High-order spatial accuracy (2nd order MUSCL)
MUSCL_RECONSTRUCTION = yes     ; 'yes' or 'no'
FLUX_LIMITER = van albada      ; 'van albada', 'van leer', 'minmod', 'superbee'

; Entropy fix for Roe-type schemes (Harten-Hyman)
ENTROPY_FIX_ACTIVE = yes
ENTROPY_FIX_COEFFICIENT = 0.2

; Simulation type: 'unsteady' (wave dynamics) or 'steady' (residual convergence)
SIMULATION_TYPE = unsteady

; Optional localized mesh refinement
MESH_REFINEMENT = no
X_START_REFINEMENT = 0.45
X_END_REFINEMENT = 0.55
NUMBER_POINTS_REFINEMENT = 200
ADAPT_MESH_REFINEMENT = yes

; Wall friction options
WALL_FRICTION_ACTIVE = no      ; 'yes' or 'no'
WALL_FRICTION_MODEL = constant ; 'constant', 'mirels_laminar', 'mirels_turbulent', 'mirels_transitional'
FRICTION_COEFFICIENT = 0.003   ; Darcy-Weisbach / Fanning factor (for constant model)
FRICTION_MAX_CF = 0.1          ; Maximum skin friction cutoff
FRICTION_DRIVER_CF = 0.003     ; Skin friction in driver gas behind contact surface
FRICTION_TRANSITION_REYNOLDS = 1.0e6 ; Transition Reynolds number

; Shock tracking & reflection settings (for Mirels friction)
SHOCK_REFLECTION_TRACKING = yes
SHOCK_DETECTION_THRESHOLD = 0.05
WALL_REFLECTION_THRESHOLD = 1.15

; Wall heat transfer
WALL_HEAT_TRANSFER_ACTIVE = no
WALL_HEAT_FLUX = 0.0           ; Prescribed heat flux [W/m^2]

; 0D Tank Boundary Parameters (if BOUNDARY_CONDITION_LEFT or RIGHT = tank)
TANK_VOLUME = 0.05             ; Reservoir volume [m^3]
TANK_INITIAL_PRESSURE = 1.0e5  ; Initial tank pressure [Pa]
TANK_INITIAL_TEMPERATURE = 300.0 ; Initial tank temperature [K]
TANK_THERMAL_MODE = adiabatic  ; 'adiabatic' or 'isothermal'

; 2D Look-Up Table (LuT) Acceleration (for real gas)
USE_LUT = no                   ; 'yes' or 'no'
LUT_PRESSURE_MIN = 5.0e5       ; Min pressure bound [Pa]
LUT_PRESSURE_MAX = 1.0e8       ; Max pressure bound [Pa]
LUT_TEMPERATURE_MIN = 250.0    ; Min temperature bound [K]
LUT_TEMPERATURE_MAX = 1500.0   ; Max temperature bound [K]
LUT_GRID_SIZE = 250, 250       ; Grid resolution (nP, nT)


[FLUID]
; Fluid name in CoolProp database (e.g., 'air', 'CO2', 'N2', 'MDM', 'MM')
FLUID_NAME = air

; Thermodynamic model: 'ideal' or 'real'
FLUID_MODEL = ideal

; Ideal gas parameters (used when FLUID_MODEL = ideal)
FLUID_GAMMA = 1.4
GAS_R_CONSTANT = 287.05

; Real gas library (used when FLUID_MODEL = real)
FLUID_LIBRARY = CoolProp


[OUTPUT]
FOLDER_NAME = Results
FILE_NAME = SimulationRun
OUTPUT_FREQUENCY = 100         ; Dump snapshot every N iterations
SHOW_ANIMATION = no            ; Show live animation during run
```

---

## Numerical Methods & Physical Models

### Governing Equations

`pyshockflow` solves the unsteady quasi-one-dimensional Euler equations with variable cross-sectional area, wall shear stress, and wall heat transfer:

$$\frac{\partial (\rho A)}{\partial t} + \frac{\partial (\rho u A)}{\partial x} = 0$$

$$\frac{\partial (\rho u A)}{\partial t} + \frac{\partial \left[(\rho u^2 + p) A\right]}{\partial x} = p \frac{\partial A}{\partial x} - \frac{1}{2} \rho u |u| C_f P_w$$

$$\frac{\partial (\rho E A)}{\partial t} + \frac{\partial \left[(\rho E + p) u A\right]}{\partial x} = q_w P_w$$

where $A(x)$ is the tube cross-sectional area, $P_w(x)$ is the wetted perimeter, $E = e + \frac{1}{2} u^2$ is the specific total energy, $C_f$ is the local skin friction coefficient, and $q_w$ is the prescribed wall heat flux.

### Riemann Solvers & Numerical Fluxes

| Solver Name | Gas Model | Formulation & Characteristics | References |
| :--- | :---: | :--- | :--- |
| **`godunov`** | Ideal | Exact iterative Riemann solver based on Toro's star-region decomposition. | Toro (2013) |
| **`roe`** | Ideal | Standard Roe approximate solver using parameter vector averaging. | Roe (1981), Toro (2013) |
| **`roe_arabi`** | Real | Direct speed-of-sound averaging $\tilde{a} = \frac{\sqrt{\rho_L} a_L + \sqrt{\rho_R} a_R}{\sqrt{\rho_L} + \sqrt{\rho_R}}$ with modified energy dissipation jump. | Arabi et al. (2017) |
| **`roe_vinokur`**| Real | Generalized projection with pressure derivatives $\chi = \left.\frac{\partial p}{\partial \rho}\right\|_e - \frac{e}{\rho}\left.\frac{\partial p}{\partial e}\right\|_\rho$ and $\kappa = \frac{1}{\rho}\left.\frac{\partial p}{\partial e}\right\|_\rho$. | Vinokur & Montagné (1990) |
| **`hllc`** | Ideal & Real | Two-wave contact-restoring Harten–Lax–van Leer solver with acoustic wave speed estimates. | Toro (2013) |
| **`ausm+up`** | Ideal & Real | All-speed pressure-based flux splitting with numerical dissipation scaling to low and high Mach numbers. | Liou (2006) |

### High-Order MUSCL & Slope Limiters

Second-order spatial accuracy is achieved via Monotonic Upstream-Centered Scheme for Conservation Laws (MUSCL) reconstruction of primitive variables $(\rho, u, p)$ at cell interfaces:

$$U_{L, i+1/2} = U_i + \frac{1}{2} \psi(r_i) (U_{i+1} - U_i), \qquad U_{R, i+1/2} = U_{i+1} - \frac{1}{2} \psi(r_{i+1}) (U_{i+2} - U_{i+1})$$

where $\psi(r)$ is evaluated with TVD limiters:
- **Van Albada**: $\psi(r) = \frac{r^2 + r}{1 + r^2}$
- **Van Leer**: $\psi(r) = \frac{r + |r|}{1 + |r|}$
- **Minmod**: $\psi(r) = \max(0, \min(1, r))$
- **Superbee**: $\psi(r) = \max(0, \max(\min(2r, 1), \min(r, 2)))$

### Thermodynamic Models & Look-Up Table (LuT)

Evaluating multi-parameter Helmholtz energy formulations (HEOS) in CoolProp at every cell interface can be computationally demanding for non-ideal fluids. `pyshockflow` features an accelerated **2D Look-Up Table (LuT)** architecture:
- **Logarithmic Tensor Grid**: Discretized on $(\log_{10} P, T)$ with uniform spacing.
- **Bicubic Spline Interpolation**: Evaluates $\rho, e, a, \chi, \kappa$ via `scipy.interpolate.RectBivariateSpline`.
- **Thermodynamic Consistency Guarantee**: Enforces $a_{i,j}^2 \equiv \chi_{i,j} + \kappa_{i,j} h_{i,j} > 0$ on all table nodes to prevent imaginary acoustic speeds and divergence across strong shock waves.
- **Vectorized Newton–Raphson Inversion**: Fast 3-iteration state inversion for $(P, \rho) \to T$ and $(\rho, e) \to P$.

### Wall Friction & Mirels Boundary Layer Correlations

Behind a moving shock wave in a shock tube, the boundary layer grows from the shock front toward the contact surface. `pyshockflow` models this spatially developing boundary layer via **Mirels correlations**:
- **Laminar (Mirels / NACA TN 3401)**:
  $$C_f(x) = \frac{0.664}{\sqrt{Re_x}}, \qquad Re_x = \frac{\rho_2 |u_s - u_2| x}{\mu_2}$$
- **Turbulent (Mirels / NACA TN 3712 / AIAA J. 1964)**:
  $$C_f(x) = \frac{0.0592}{Re_x^{0.2}}$$
- **Singularity Regularization**: Over the cell containing the shock front $[0, \Delta x]$, the integrable boundary-layer singularity is regularized analytically:
  $$\bar{C}_f = \frac{1}{\Delta x} \int_0^{\Delta x} C_0 x^{-n} \, dx = \frac{1}{1 - n} C_f(\Delta x)$$
  yielding $\bar{C}_f = 2.0 C_f(\Delta x)$ for laminar flow and $\bar{C}_f = 1.25 C_f(\Delta x)$ for turbulent flow.

### Bidirectional & Reflected Shock Tracking

The integrated `ShockTracker` monitors the flow field and operates a finite-state machine with four tracking modes:
1. `incident_right`: Incident shock traveling right toward $x = L$.
2. `reflected_left`: Shock reflected from the right solid end-wall, propagating left into post-shock State 2 gas.
3. `incident_left`: Incident shock traveling left toward $x = 0$.
4. `reflected_right`: Shock reflected from the left solid end-wall, propagating right.

The tracker simultaneously identifies the contact surface and measures shock speed $u_s$, upstream state (State 1), post-shock state (State 2), and post-reflected state (State 5).

### 0D Lumped-Parameter Tank Boundary Condition

The `Tank0D` module couples a finite-volume reservoir of volume $V$, initial pressure $p_0$, and temperature $T_0$ to either boundary of the shock tube. It dynamically models charging (tube emptying into tank) and discharge (tank blowdown into tube) by solving:

$$\frac{dm}{dt} = \dot{m}_{\text{interface}}, \qquad \frac{dU}{dt} = \dot{m}_{\text{interface}} h_{tot}$$

Exported time histories (`TankHistory_Left.csv`, `TankHistory_Right.csv`) record reservoir pressure, temperature, density, mass, and mass flow rate over time.

### Quasi-1D Nozzles & Heat Transfer

For variable cross-section geometries (such as convergent-divergent rocket nozzles or supersonic diffusers), `pyshockflow` computes the geometrical source terms directly from the area gradient $\frac{dA}{dx}$:

$$S_{\text{geom}} = \left[ -\rho u \frac{1}{A}\frac{dA}{dx}, \ -\rho u^2 \frac{1}{A}\frac{dA}{dx}, \ -u (\rho E + p) \frac{1}{A}\frac{dA}{dx} \right]^T$$

Wall heat flux $q_w$ $[\text{W/m}^2]$ can be coupled to simulate heated or cooled channel flows.

---

## Verification & Automated Test Suite

`pyshockflow` is continuously validated using an automated test framework:

```bash
# Run all regression benchmarks against reference solutions
pytest test/regression_tests/test_cases.py -v

# Run targeted unit tests
pytest test/unit_tests/ -v
```

The regression suite tests end-to-end physics against verified reference solutions with strict relative tolerances:
- **`sod_test_roe`**: Standard Sod shock tube with Roe solver.
- **`sod_test_godunov`**: Sod shock tube with exact Godunov Riemann solver.
- **`mesh_refinement`**: Localized grid refinement across discontinuities.
- **`nozzle_steady`**: Transonic de Laval nozzle with normal shock capture.
- **`thruster`**: Supersonic rocket thruster expansion.
- **`nonideal_CO2`**: Supercritical CO₂ shock tube comparing Arabi and Vinokur formulations.
- **`duct_with_friction`**: Duct flow with wall shear stress source term.
- **`duct_with_heatflux`**: Duct flow with wall thermal energy source term.

Unit tests verify individual components including the AUSM+-up and HLLC solvers, ideal and real fluid property routines, Mirels friction models, and 0D tank ODE integration.

---

## Repository Structure

```
pyshockflow/
├── docs/                               # Detailed technical notes and optimization reports
│   └── gemini_assisted_additions.md
├── images/                             # Verification and benchmark figures
│   ├── godunov_idealgas.png
│   └── co2_validation.png
├── src/
│   └── pyshockflow/                    # Core solver package
│       ├── __init__.py                 # Exported classes and functions
│       ├── config.py                   # INI configuration parser
│       ├── driver.py                   # Main finite-volume solver & time integration
│       ├── fluid.py                    # Fluid thermodynamics, CoolProp & 2D LuT
│       ├── riemann_problem.py          # Exact Riemann problem solver
│       ├── math_utils.py               # Primitive-conservative conversions
│       ├── roe_vectorized.py           # Vectorized Roe solvers (Toro, Arabi, Vinokur)
│       ├── hllc_vectorized.py          # Vectorized HLLC solver
│       ├── ausm_vectorized.py          # Vectorized AUSM+-up solver
│       ├── kernels_numba.py            # Numba LLVM JIT single-pass fused kernels
│       ├── friction_models.py          # Constant & Mirels friction, ShockTracker
│       ├── tank_boundary.py            # Coupled 0D lumped-parameter tank model
│       ├── output.py                   # Snapshot assembler & animation tools
│       ├── nicfd_styles.py             # Publication plotting style (NICFD template)
│       └── thesis_plots.py             # Academic plotting templates
├── test/
│   ├── regression_tests/               # Automated regression suite
│   │   └── test_cases.py
│   └── unit_tests/                     # Unit tests for schemes, fluids, friction, tank
├── testcases/                          # Ready-to-run simulation cases
│   ├── nonideal_CO2_lut_comparison/    # Exact Helmholtz vs 2D LuT benchmark
│   ├── nonideal_C02/                   # Real-gas scheme comparison (Arabi, Vinokur, HLLC, AUSM)
│   ├── nozzle_steady/                  # Steady nozzle with normal shock
│   ├── run_tank/                       # Shock tube coupled to 0D tank
│   ├── shock_tube_with_friction/       # Constant vs Mirels shock-induced friction
│   ├── thruster/                       # Supersonic thruster expansion
│   └── validation_ideal/               # Toro 5 tests validation
├── environment.yml                     # Conda environment specification
├── pyproject.toml                      # Package installation metadata
├── LICENSE                             # GNU General Public License v3.0
└── README.md                           # Main documentation
```

---

## Authors & Contact

- **Francesco Neri** — *PhD Candidate*, Propulsion and Power, TU Delft  
  Email: `f.neri@tudelft.nl` | [GitHub](https://github.com/francesconeri95)
- **Matteo Pini** — *Associate Professor*, Propulsion and Power, TU Delft  
  Email: `m.pini@tudelft.nl`

**Affiliation**:  
*Propulsion and Power Group, Faculty of Aerospace Engineering*  
*Delft University of Technology (TU Delft), Kluyverweg 1, 2629 HS Delft, The Netherlands*

---

## Citation & References

If you use `pyshockflow` in your research, academic thesis, or publications, please cite the repository:

```bibtex
@software{pyshockflow2026,
  author = {Neri, Francesco and Pini, Matteo},
  title = {pyshockflow: A High-Performance Finite-Volume Solver for Quasi-1D Ideal and Non-Ideal Compressible Fluid Dynamics},
  url = {https://github.com/Propulsion-Power-TU-Delft/pyshockflow},
  version = {0.1.0},
  year = {2026}
}
```

### Key Literature

1. **Toro, E. F.** (2013). *Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction*. Springer Science & Business Media.
2. **Arabi, S., Trépanier, J.-Y., & Camarero, R.** (2017). "A simple extension of Roe's scheme for real gases." *Journal of Computational Physics*, 329, 16–28.
3. **Vinokur, M., & Montagné, J.-L.** (1990). "Generalized Flux-Vector Splitting and Roe-Average for an Arbitrary Gas." *Journal of Computational Physics*, 89(2), 276–300.
4. **Liou, M.-S.** (2006). "A sequel to AUSM, Part II: $\text{AUSM}^+$-up for all speeds." *Journal of Computational Physics*, 214(1), 137–170.
5. **Mirels, H.** (1955). *Laminar boundary layer behind shock advancing into stationary fluid*. NACA TN 3401.
6. **Mirels, H.** (1956). *Boundary layer behind shock or thin expansion wave moving into stationary fluid*. NACA TN 3712.
7. **Mirels, H.** (1964). "Shock tube test time limitation due to turbulent-wall boundary layer." *AIAA Journal*, 2(1), 84–93.
8. **Vimercati, D., & Guardone, A.** (2018). "On the numerical simulation of non-classical quasi-1D steady nozzle flows: Capturing sonic shocks." *Computers & Fluids*, 172, 608–621.
9. **Blazek, J.** (2015). *Computational Fluid Dynamics: Principles and Applications*. Butterworth-Heinemann.
10. **Bell, I. H., et al.** (2014). "Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp." *Industrial & Engineering Chemistry Research*, 53(6), 2498–2508.

---

## License

This project is licensed under the **GNU General Public License v3.0** — see the [LICENSE](LICENSE) file for details.