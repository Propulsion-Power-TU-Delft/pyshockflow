"""
Friction models for 1D compressible shock tube flows.
Implements constant pipe friction and Mirels shock-boundary-layer correlations
(NACA TN 3401, NACA TN 3712, AIAA J. 1964).
"""

import numpy as np


class ShockTracker:
    """
    Detects and tracks moving shock waves and contact surfaces in 1D flow fields.
    Computes shock location, direction, instantaneous shock speed via Rankine-Hugoniot,
    and post-shock (State 2) / pre-shock (State 1) properties.
    """
    def __init__(self, pressure_threshold=0.05):
        self.pressure_threshold = pressure_threshold
        self.shock_detected = False
        self.x_shock = 0.0
        self.i_shock = 0
        self.direction = 1  # +1: right-running, -1: left-running
        self.u_shock = 0.0
        self.x_contact = 0.0
        self.i_contact = 0
        self.rho1 = 1.0
        self.u1 = 0.0
        self.p1 = 1.0
        self.rho2 = 1.0
        self.u2 = 0.0
        self.p2 = 1.0

    def detect(self, primitive, x_nodes):
        """
        Analyze primitive solution arrays to detect shock front and contact surface.

        Parameters
        ----------
        primitive : dict
            Primitive variables with keys 'Density', 'Velocity', 'Pressure'.
        x_nodes : np.ndarray
            Spatial node coordinates (including halo cells).

        Returns
        -------
        bool
            True if a shock was successfully identified, False otherwise.
        """
        rho = primitive['Density']
        u = primitive['Velocity']
        p = primitive['Pressure']
        n = len(p)

        if n < 6:
            self.shock_detected = False
            return False

        # Compute gradients in the internal physical domain (avoiding halo boundary cells)
        dp = np.diff(p)
        du = np.diff(u)
        dx = np.diff(x_nodes)
        dx = np.where(dx == 0.0, 1e-12, dx)

        # Compression indicator: negative velocity gradient accompanied by pressure jump
        compression = -du / dx
        p_avg = 0.5 * (p[:-1] + p[1:])
        p_jump = np.abs(dp) / np.maximum(p_avg, 1e-12)

        # Filter by pressure threshold to distinguish shocks from acoustic noise
        valid_shock_cells = (p_jump > self.pressure_threshold) & (compression > 0.0)

        if not np.any(valid_shock_cells):
            self.shock_detected = False
            return False

        # Shock index is at the peak compression among valid shock cells
        metric = compression * p_jump
        metric[~valid_shock_cells] = -1.0
        i_jump = int(np.argmax(metric))

        # Determine propagation direction from pressure jump sign
        # dp[i] = p[i+1] - p[i]
        # For right-running shock: p is higher on the left, so p[i] > p[i+1] -> dp[i] < 0
        if dp[i_jump] < 0:
            self.direction = 1  # right-running
            i_shock = i_jump + 1
            i_behind = max(1, i_shock - 2)
            i_ahead = min(n - 2, i_shock + 2)
            self.rho2, self.u2, self.p2 = rho[i_behind], u[i_behind], p[i_behind]
            self.rho1, self.u1, self.p1 = rho[i_ahead], u[i_ahead], p[i_ahead]
        else:
            self.direction = -1  # left-running
            i_shock = i_jump
            i_behind = min(n - 2, i_shock + 2)
            i_ahead = max(1, i_shock - 2)
            self.rho2, self.u2, self.p2 = rho[i_behind], u[i_behind], p[i_behind]
            self.rho1, self.u1, self.p1 = rho[i_ahead], u[i_ahead], p[i_ahead]

        self.i_shock = i_shock
        self.x_shock = x_nodes[i_shock]

        # Instantaneous shock velocity from Rankine-Hugoniot mass conservation jump:
        # rho1 * (u1 - us) = rho2 * (u2 - us) => us = (rho2*u2 - rho1*u1) / (rho2 - rho1)
        d_rho = self.rho2 - self.rho1
        if abs(d_rho) > 1e-6:
            self.u_shock = (self.rho2 * self.u2 - self.rho1 * self.u1) / d_rho
        else:
            self.u_shock = self.u2

        # Detect contact surface (density gradient where pressure and velocity gradients are small)
        self._detect_contact_surface(rho, u, p, x_nodes)

        self.shock_detected = True
        return True

    def _detect_contact_surface(self, rho, u, p, x_nodes):
        """Identify contact surface separating State 2 (test gas) from State 3 (driver gas)."""
        if self.direction == 1:
            # Search upstream of the shock (cells 1 to i_shock - 2)
            search_end = max(2, self.i_shock - 2)
            if search_end > 2:
                drho = np.abs(np.diff(rho[1:search_end]))
                i_contact = 1 + int(np.argmax(drho))
                self.i_contact = i_contact
                self.x_contact = x_nodes[i_contact]
            else:
                self.i_contact = 1
                self.x_contact = x_nodes[1]
        else:
            # Search downstream of the shock for left-running waves
            search_start = min(len(rho) - 3, self.i_shock + 2)
            if search_start < len(rho) - 2:
                drho = np.abs(np.diff(rho[search_start:-1]))
                i_contact = search_start + int(np.argmax(drho))
                self.i_contact = i_contact
                self.x_contact = x_nodes[i_contact]
            else:
                self.i_contact = len(rho) - 2
                self.x_contact = x_nodes[-2]


class FrictionModel:
    """Base abstract class for wall friction models."""
    def compute_friction_coefficients(self, primitive, x_nodes, dx, fluid):
        """Compute local skin friction coefficient Cf across all nodes."""
        raise NotImplementedError

    def compute_source_terms(self, primitive, x_nodes, dx, hydraulic_diameter, fluid):
        """Compute momentum source term Sm = -4*tau_w / Dh [N/m^3] for all nodes."""
        cf = self.compute_friction_coefficients(primitive, x_nodes, dx, fluid)
        rho = primitive['Density']
        u = primitive['Velocity']
        # Sm = -4 * (0.5 * cf * rho * u * |u|) / Dh = -2.0 * cf * rho * u * |u| / Dh
        sm = -2.0 * cf * rho * u * np.abs(u) / hydraulic_diameter
        source = np.zeros((len(u), 3))
        source[:, 1] = sm
        return source


class ConstantFriction(FrictionModel):
    """
    Standard constant friction factor model (legacy behavior).
    Cf is uniform across the entire domain.
    """
    def __init__(self, cf_constant=0.003):
        self.cf_constant = cf_constant

    def compute_friction_coefficients(self, primitive, x_nodes, dx, fluid):
        n = len(primitive['Velocity'])
        return np.full(n, self.cf_constant)


class MirelsFrictionBase(FrictionModel):
    """
    Base class for Mirels shock boundary layer correlations.
    Handles shock tracking, distance calculation, region-of-validity masking,
    and leading-edge finite-volume singularity cell averaging.
    """
    def __init__(self, pressure_threshold=0.05, max_cf=0.1, driver_cf=0.003):
        self.tracker = ShockTracker(pressure_threshold=pressure_threshold)
        self.max_cf = max_cf
        self.driver_cf = driver_cf

    def _compute_reynolds(self, x_dist, rho_post, u_post, u_shock, mu_post):
        """
        Shock-attached Reynolds number:
        Re_x = rho_2 * |u_shock - u_post| * x / mu_2
        """
        u_rel = abs(u_shock - u_post)
        mu_safe = max(mu_post, 1e-12)
        return (rho_post * u_rel * x_dist) / mu_safe

    def _regularize_singularity(self, cf_edge, n_exponent):
        """
        Finite-volume cell average of Cf(x) = C0 * x^(-n) over [0, dx]:
        bar_Cf = 1 / (1 - n) * Cf(dx)
        For laminar (n=0.5): bar_Cf = 2.0 * Cf(dx)
        For turbulent (n=0.2): bar_Cf = 1.25 * Cf(dx)
        """
        factor = 1.0 / (1.0 - n_exponent)
        return factor * cf_edge

    def compute_friction_coefficients(self, primitive, x_nodes, dx, fluid):
        """Compute spatially developing skin friction coefficient Cf."""
        n = len(x_nodes)
        cf = np.zeros(n)

        shock_found = self.tracker.detect(primitive, x_nodes)
        if not shock_found:
            # Fallback to driver/baseline friction if no shock detected yet
            return np.full(n, self.driver_cf)

        # Post-shock properties
        rho2 = self.tracker.rho2
        u2 = self.tracker.u2
        p2 = self.tracker.p2
        u_shock = self.tracker.u_shock
        direction = self.tracker.direction
        x_shock = self.tracker.x_shock
        x_contact = self.tracker.x_contact

        # Dynamic viscosity of State 2
        try:
            mu2 = fluid.computeViscosity_p_rho(p2, rho2)
        except Exception:
            mu2 = 1.8e-5

        # Classify cells into ahead, test gas (Region 2), and driver gas (Region 3)
        for i in range(n):
            xi = x_nodes[i]
            dxi = dx[i]

            if direction == 1:
                # Right-running shock: shock at x_shock, moving right
                if xi > x_shock:
                    # Undisturbed quiescent gas ahead of the shock
                    cf[i] = 0.0
                elif xi >= x_contact:
                    # Test gas behind shock: distance x measured from shock front
                    dist = x_shock - xi
                    cf[i] = self._evaluate_cf_local(dist, dxi, rho2, u2, u_shock, mu2)
                else:
                    # Driver gas behind contact surface
                    cf[i] = self.driver_cf
            else:
                # Left-running shock: shock at x_shock, moving left
                if xi < x_shock:
                    cf[i] = 0.0
                elif xi <= x_contact:
                    dist = xi - x_shock
                    cf[i] = self._evaluate_cf_local(dist, dxi, rho2, u2, u_shock, mu2)
                else:
                    cf[i] = self.driver_cf

        # Apply maximum Cf numerical limiter clamp
        cf = np.clip(cf, 0.0, self.max_cf)
        return cf

    def _evaluate_cf_local(self, dist, dxi, rho2, u2, u_shock, mu2):
        """Override in subclasses to evaluate specific Cf correlations."""
        raise NotImplementedError


class MirelsLaminarFriction(MirelsFrictionBase):
    """
    Mirels Laminar Boundary Layer Model (NACA TN 3401 / Blasius approximation):
    Cf(x) = 0.664 / sqrt(Re_x)
    Singularity regularization over shock cell [0, dx]:
    bar_Cf = 2.0 * Cf(dx)
    """
    def _evaluate_cf_local(self, dist, dxi, rho2, u2, u_shock, mu2):
        n_exp = 0.5
        if dist < dxi:
            # First cell behind the shock: analytical finite-volume cell average
            re_dx = self._compute_reynolds(dxi, rho2, u2, u_shock, mu2)
            re_safe = max(re_dx, 1e-4)
            cf_edge = 0.664 / np.sqrt(re_safe)
            return self._regularize_singularity(cf_edge, n_exp)
        else:
            # Standard cell center evaluation
            re_x = self._compute_reynolds(dist, rho2, u2, u_shock, mu2)
            re_safe = max(re_x, 1e-4)
            return 0.664 / np.sqrt(re_safe)


class MirelsTurbulentFriction(MirelsFrictionBase):
    """
    Mirels Turbulent Boundary Layer Model (NACA TN 3712 / AIAA J. 1964):
    Cf(x) = 0.0592 / Re_x^0.2
    Singularity regularization over shock cell [0, dx]:
    bar_Cf = 1.25 * Cf(dx)
    """
    def _evaluate_cf_local(self, dist, dxi, rho2, u2, u_shock, mu2):
        n_exp = 0.2
        if dist < dxi:
            # First cell behind the shock: analytical finite-volume cell average
            re_dx = self._compute_reynolds(dxi, rho2, u2, u_shock, mu2)
            re_safe = max(re_dx, 1e-4)
            cf_edge = 0.0592 / (re_safe ** 0.2)
            return self._regularize_singularity(cf_edge, n_exp)
        else:
            # Standard cell center evaluation
            re_x = self._compute_reynolds(dist, rho2, u2, u_shock, mu2)
            re_safe = max(re_x, 1e-4)
            return 0.0592 / (re_safe ** 0.2)


class MirelsTransitionalFriction(MirelsFrictionBase):
    """
    Mirels Transitional Boundary Layer Model:
    Switches from laminar (0.664 / sqrt(Re_x)) to turbulent (0.0592 / Re_x^0.2)
    when Re_x >= re_transition.
    """
    def __init__(self, re_transition=1.0e6, pressure_threshold=0.05, max_cf=0.1, driver_cf=0.003):
        super().__init__(pressure_threshold=pressure_threshold, max_cf=max_cf, driver_cf=driver_cf)
        self.re_transition = re_transition

    def _evaluate_cf_local(self, dist, dxi, rho2, u2, u_shock, mu2):
        if dist < dxi:
            # Shock cell is always laminar (initiates at x=0)
            re_dx = self._compute_reynolds(dxi, rho2, u2, u_shock, mu2)
            re_safe = max(re_dx, 1e-4)
            cf_edge = 0.664 / np.sqrt(re_safe)
            return self._regularize_singularity(cf_edge, 0.5)
        else:
            re_x = self._compute_reynolds(dist, rho2, u2, u_shock, mu2)
            re_safe = max(re_x, 1e-4)
            if re_safe < self.re_transition:
                return 0.664 / np.sqrt(re_safe)
            else:
                return 0.0592 / (re_safe ** 0.2)


def create_friction_model(config):
    """
    Factory function to instantiate the configured wall friction model.

    Parameters
    ----------
    config : Config
        Simulation configuration object.

    Returns
    -------
    FrictionModel
        Instantiated friction model instance.
    """
    model_name = config.getWallFrictionModel().lower().strip()
    driver_cf = config.getFrictionDriverCf()
    max_cf = config.getFrictionMaxCf()
    p_threshold = config.getShockDetectionThreshold()

    if model_name in ['mirels_laminar', 'laminar']:
        return MirelsLaminarFriction(pressure_threshold=p_threshold, max_cf=max_cf, driver_cf=driver_cf)
    elif model_name in ['mirels_turbulent', 'turbulent']:
        return MirelsTurbulentFriction(pressure_threshold=p_threshold, max_cf=max_cf, driver_cf=driver_cf)
    elif model_name in ['mirels_transitional', 'transitional']:
        re_tr = config.getFrictionTransitionReynolds()
        return MirelsTransitionalFriction(re_transition=re_tr, pressure_threshold=p_threshold, max_cf=max_cf, driver_cf=driver_cf)
    else:
        # Default fallback to constant friction for backward compatibility
        cf_const = config.getFrictionCoefficient()
        return ConstantFriction(cf_constant=cf_const)
