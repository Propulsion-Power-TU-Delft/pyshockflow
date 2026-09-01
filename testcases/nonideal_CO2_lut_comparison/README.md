# Real-Gas CO2 Shock Tube: Look-Up Table (LuT) vs Exact Helmholtz EoS

This test case compares the performance and solution accuracy of a **1000-point real-gas shock tube simulation** for supercritical carbon dioxide ($CO_2$) solved with:
1. **Exact Analytical Helmholtz EoS** (`input_exact.ini`, `USE_LUT = no`)
2. **2D Bicubic Spline Look-Up Table** (`input_lut.ini`, `USE_LUT = yes`)

---

## Shock Tube Initial Conditions
- **Left State** (High-Pressure Supercritical Driver):
  - $P_L = 73.76\text{ MPa}$ ($737.6\text{ bar}$)
  - $\rho_L = 348.8\text{ kg/m}^3$
  - $u_L = 0\text{ m/s}$
- **Right State** (Low-Pressure Driven Gas):
  - $P_R = 0.7376\text{ MPa}$ ($7.376\text{ bar}$)
  - $\rho_R = 3.488\text{ kg/m}^3$
  - $u_R = 0\text{ m/s}$
- **Grid Resolution**: 1000 uniform cells ($L = 10\text{ m}$)
- **Time Horizon**: $t_{\text{final}} = 2.75\text{ ms}$
- **Numerical Scheme**: Arabi et al. (2017) real-gas Roe solver

---

## Input File Configuration

### Case 1: Exact Helmholtz (`input_exact.ini`)
```ini
[SIMULATION]
NUMBER_POINTS = 1000
USE_LUT = no
...
```

### Case 2: 2D Bicubic Spline Look-Up Table (`input_lut.ini`)
```ini
[SIMULATION]
NUMBER_POINTS = 1000
USE_LUT = yes

[FLUID]
FLUID_NAME = CO2
FLUID_MODEL = real
FLUID_LIBRARY = CoolProp
LUT_PRESSURE_MIN = 5e5         ; Minimum pressure [Pa]
LUT_PRESSURE_MAX = 1e8         ; Maximum pressure [Pa]
LUT_TEMPERATURE_MIN = 250      ; Minimum temperature [K]
LUT_TEMPERATURE_MAX = 1500     ; Maximum temperature [K]
LUT_GRID_SIZE = 250, 250       ; Grid resolution: 250 x 250 nodes
```

---

## Running the Benchmark

To run both simulations, print the timing and speedup metrics, and check solution accuracy:

```bash
python main.py
```

To generate comparison plots:

```bash
python plot_comparison.py
```
