"""
0D Lumped Tank Boundary Condition for 1D Shock Tube Flows.
Integrates mass and internal energy conservation ODEs for a finite volume reservoir
connected to the tube boundary, dynamically providing time-evolving backpressure.
"""

import numpy as np
import csv
from pathlib import Path


class Tank0D:
    """
    0D lumped parameter rigid tank model coupled to tube boundaries.

    Governing ODEs:
      dm/dt = m_dot
      dU/dt = m_dot * h_tot

    Parameters
    ----------
    volume : float
        Tank internal volume [m^3].
    p_initial : float
        Initial static pressure in the tank [Pa].
    T_initial : float
        Initial temperature in the tank [K].
    fluid : FluidIdeal or FluidReal
        Fluid thermodynamic model.
    location : str, optional
        Boundary location: 'right' (default) or 'left'.
    thermal_mode : str, optional
        'adiabatic' (default) or 'isothermal'.
    """
    def __init__(self, volume, p_initial, T_initial, fluid, location='right', thermal_mode='adiabatic'):
        self.volume = float(volume)
        self.p0 = float(p_initial)
        self.T0 = float(T_initial)
        self.fluid = fluid
        self.location = str(location).lower()
        self.thermal_mode = str(thermal_mode).lower()

        # Compute initial thermodynamic state
        self.rho0 = float(self.fluid.computeDensity_p_T(self.p0, self.T0))
        self.e0 = float(self.fluid.computeStaticEnergy_p_rho(self.p0, self.rho0))
        self.m0 = self.rho0 * self.volume
        self.U0 = self.m0 * self.e0

        # Current state
        self.p = self.p0
        self.T = self.T0
        self.rho = self.rho0
        self.e = self.e0
        self.m = self.m0
        self.U = self.U0

        # Interface flux storage for the current time step
        self.current_mdot = 0.0
        self.current_htot = 0.0

        # History records
        self.time_history = [0.0]
        self.pressure_history = [self.p0]
        self.temperature_history = [self.T0]
        self.density_history = [self.rho0]
        self.mass_history = [self.m0]
        self.mdot_history = [0.0]

    def apply_boundary_condition(self, solutionPrimitive, iInternal, iHalo, area_tube, fluid):
        """
        Evaluate and set the halo cell primitive state based on the current tank state,
        local Mach number, and flow direction.

        Parameters
        ----------
        solutionPrimitive : dict
            Primitive variables dictionary.
        iInternal : int
            Index of the adjacent physical cell (1 for left, -2 for right).
        iHalo : int
            Index of the halo cell (0 for left, -1 for right).
        area_tube : float or np.ndarray
            Cross-sectional area at the interface [m^2].
        fluid : FluidIdeal or FluidReal
            Fluid model.
        """
        u_int = solutionPrimitive['Velocity'][iInternal]
        p_int = solutionPrimitive['Pressure'][iInternal]
        rho_int = solutionPrimitive['Density'][iInternal]
        e_int = solutionPrimitive['Energy'][iInternal]

        area = float(area_tube[iInternal] if hasattr(area_tube, '__len__') else area_tube)
        mach_int = fluid.computeMach_u_p_rho(u_int, p_int, rho_int)

        # Determine whether flow is entering tank (tube charging tank)
        # or exiting tank into tube (tank discharging into tube)
        if self.location == 'right':
            is_charging = (p_int >= self.p) or (u_int > 0.0)
        else:
            is_charging = (p_int >= self.p) or (u_int < 0.0)

        if is_charging:
            # Flow is leaving the tube and charging the tank
            if mach_int >= 1.0 and self.p <= p_int:
                # Under-expanded supersonic outflow: transparent condition
                # (all characteristics leave domain, downstream pressure cannot propagate upstream)
                p_halo = p_int
                rho_halo = rho_int
                u_halo = u_int
                e_halo = e_int
            else:
                # Subsonic outflow OR adverse backpressure on supersonic flow (self.p > p_int):
                # The tank imposes its static backpressure on the halo cell.
                # When self.p overcomes the supersonic dynamic head, the Riemann solver
                # naturally generates an upstream-propagating shock wave, unchoking the nozzle.
                p_halo = self.p
                rho_halo = rho_int
                u_halo = u_int
                e_halo = fluid.computeStaticEnergy_p_rho(p_halo, rho_halo)

            mdot = rho_int * abs(u_int) * area if ((self.location == 'right' and u_int > 0) or (self.location == 'left' and u_int < 0)) else 0.0
            h_tot = e_int + p_int / max(rho_int, 1e-12) + 0.5 * (u_int ** 2)
            self.current_mdot = mdot
            self.current_htot = h_tot

        else:
            # Backflow: flow is entering the tube from the tank (tank discharges)
            direction = -1.0 if self.location == 'right' else 1.0
            p_tank_tot = self.p
            T_tank_tot = self.T

            try:
                rho_halo, u_halo, e_halo = fluid.computeInletQuantities(
                    p_int, p_tank_tot, T_tank_tot, direction
                )
                p_halo = p_int
            except Exception:
                p_halo = self.p
                rho_halo = self.rho
                u_halo = 0.0
                e_halo = self.e

            # Mass leaving tank is negative for dm/dt = m_dot
            mdot = -rho_halo * abs(u_halo) * area
            h_tot = self.e + self.p / max(self.rho, 1e-12)
            self.current_mdot = mdot
            self.current_htot = h_tot

        # Set halo states
        solutionPrimitive['Density'][iHalo] = rho_halo
        solutionPrimitive['Velocity'][iHalo] = u_halo
        solutionPrimitive['Pressure'][iHalo] = p_halo
        solutionPrimitive['Energy'][iHalo] = e_halo

    def step(self, dt, current_time, mass_flux=None, energy_flux=None):
        """
        Advance the 0D tank ODEs by time step dt using Forward Euler.

        Parameters
        ----------
        dt : float
            Time step duration [s].
        current_time : float
            Current simulation time instant [s].
        mass_flux : float, optional
            Exact interface mass flow rate entering the tank [kg/s].
            If None, uses current_mdot computed from cell states.
        energy_flux : float, optional
            Exact interface total energy flow rate entering the tank [J/s].
            If None, uses current_mdot * current_htot.
        """
        if mass_flux is not None:
            self.current_mdot = float(mass_flux)

        # Mass update
        m_new = max(self.m + self.current_mdot * dt, 1e-9)
        rho_new = m_new / self.volume

        # Energy & thermodynamic update
        if self.thermal_mode == 'isothermal':
            T_new = self.T0
            if hasattr(self.fluid, 'gmma'):
                p_new = rho_new * self.fluid.Rgas * T_new
                e_new = self.fluid.computeStaticEnergy_p_rho(p_new, rho_new)
            else:
                p_new = self.fluid.computePressure_rho_e(rho_new, self.e0)
                e_new = self.e0
            U_new = m_new * e_new
        else:
            # Adiabatic charging/discharging
            if energy_flux is not None:
                U_new = max(self.U + float(energy_flux) * dt, 1e-9)
            else:
                U_new = max(self.U + self.current_mdot * self.current_htot * dt, 1e-9)
            e_new = U_new / m_new

            if hasattr(self.fluid, 'gmma'):
                # Ideal gas EOS: p = (gamma - 1) * rho * e
                p_new = (self.fluid.gmma - 1.0) * rho_new * e_new
                T_new = p_new / max(rho_new * self.fluid.Rgas, 1e-12)
            else:
                # Real gas EOS
                p_new = self.fluid.computePressure_rho_e(rho_new, e_new)
                try:
                    T_new = self.fluid.computeTemperature_p_rho(p_new, rho_new)
                except Exception:
                    T_new = self.T

        # Safeguards
        p_new = max(float(p_new), 1e-3)
        T_new = max(float(T_new), 1.0)

        # Update state
        self.m = m_new
        self.rho = rho_new
        self.U = U_new
        self.e = e_new
        self.p = p_new
        self.T = T_new

        # Append to time-series history
        self.time_history.append(current_time)
        self.pressure_history.append(self.p)
        self.temperature_history.append(self.T)
        self.density_history.append(self.rho)
        self.mass_history.append(self.m)
        self.mdot_history.append(self.current_mdot)

    def save_history(self, output_folder, file_name='TankHistory.csv'):
        """
        Export tank history time series to CSV.

        Parameters
        ----------
        output_folder : str or Path
            Destination directory.
        file_name : str, optional
            Output CSV filename.
        """
        out_path = Path(output_folder) / file_name
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(out_path, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Time_s', 'Pressure_Pa', 'Temperature_K', 'Density_kgm3', 'Mass_kg', 'MassFlowRate_kgs'])
            for t, p, T, rho, m, mdot in zip(
                self.time_history,
                self.pressure_history,
                self.temperature_history,
                self.density_history,
                self.mass_history,
                self.mdot_history
            ):
                writer.writerow([f"{t:.8e}", f"{p:.8e}", f"{T:.8e}", f"{rho:.8e}", f"{m:.8e}", f"{mdot:.8e}"])
