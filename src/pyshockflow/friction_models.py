"""
Friction models for 1D compressible shock tube flows.
Implements constant pipe friction and Mirels shock-boundary-layer correlations
(NACA TN 3401, NACA TN 3712, AIAA J. 1964) with bidirectional and reflected wave tracking.
"""

import numpy as np


class ShockTracker:
    """
    Tracks incident and reflected shock waves and contact surfaces in 1D flow fields.
    Maintains an internal wave state machine to robustly handle end-wall reflections:
      - 'incident_right': Shock moving right into quiescent gas toward x = L.
      - 'reflected_left': Shock reflected off right end-wall, moving left into State 2.
      - 'incident_left': Shock moving left into quiescent gas toward x = 0.
      - 'reflected_right': Shock reflected off left end-wall, moving right into State 2.
    """
    def __init__(self, pressure_threshold=0.05, enable_reflection=True, reflection_threshold=1.15):
        self.pressure_threshold = pressure_threshold
        self.enable_reflection = enable_reflection
        self.reflection_threshold = reflection_threshold

        self.shock_detected = False
        self.wave_mode = 'incident_right'
        self.x_shock = 0.0
        self.i_shock = 0
        self.direction = 1  # +1: right-running, -1: left-running
        self.u_shock = 0.0
        self.x_contact = 0.0
        self.i_contact = 0

        # State 1: Pre-incident shock (quiescent test gas)
        self.rho1 = 1.0
        self.u1 = 0.0
        self.p1 = 1.0

        # State 2: Post-incident shock / oncoming test gas slug
        self.rho2 = 1.0
        self.u2 = 0.0
        self.p2 = 1.0

        # State 5: Post-reflected shock (compressed, near-stagnant gas at wall)
        self.rho5 = 1.0
        self.u5 = 0.0
        self.p5 = 1.0

        # Wall location and peak tracking
        self.x_wall = 0.0
        self.peak_post_shock_p = 0.0

    def detect(self, primitive, x_nodes):
        """
        Analyze primitive solution arrays to track shock wave and contact surface,
        automatically transitioning states upon end-wall reflection.

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

        dp = np.diff(p)
        du = np.diff(u)
        dx = np.diff(x_nodes)
        dx = np.where(dx == 0.0, 1e-12, dx)

        compression = -du / dx
        p_avg = 0.5 * (p[:-1] + p[1:])
        p_jump = np.abs(dp) / np.maximum(p_avg, 1e-12)

        # Dispatch based on current state machine mode
        if self.wave_mode == 'incident_right':
            return self._detect_incident_right(rho, u, p, dp, compression, p_jump, x_nodes, n)
        elif self.wave_mode == 'reflected_left':
            return self._detect_reflected_left(rho, u, p, dp, compression, p_jump, x_nodes, n)
        elif self.wave_mode == 'incident_left':
            return self._detect_incident_left(rho, u, p, dp, compression, p_jump, x_nodes, n)
        elif self.wave_mode == 'reflected_right':
            return self._detect_reflected_right(rho, u, p, dp, compression, p_jump, x_nodes, n)
        else:
            return self._detect_incident_right(rho, u, p, dp, compression, p_jump, x_nodes, n)

    def _detect_incident_right(self, rho, u, p, dp, compression, p_jump, x_nodes, n):
        """Track incident right-running shock and monitor for right end-wall collision."""
        # Check if right wall reflection already occurred
        if self.enable_reflection and self.peak_post_shock_p > 0.0:
            p_right_wall = np.max(p[-4:])
            if p_right_wall >= self.reflection_threshold * self.peak_post_shock_p:
                self.wave_mode = 'reflected_left'
                self.x_wall = x_nodes[-2]
                self.direction = -1
                return self._detect_reflected_left(rho, u, p, dp, compression, p_jump, x_nodes, n)

        valid = (p_jump > self.pressure_threshold) & (compression > 0.0) & (dp < 0.0)

        # Check for initial orientation: if shock is actually left-running at the start
        if not np.any(valid):
            left_running_valid = (p_jump > self.pressure_threshold) & (compression > 0.0) & (dp > 0.0)
            if np.any(left_running_valid) and self.peak_post_shock_p == 0.0:
                self.wave_mode = 'incident_left'
                return self._detect_incident_left(rho, u, p, dp, compression, p_jump, x_nodes, n)
            self.shock_detected = False
            return False

        metric = compression * p_jump
        metric[~valid] = -1.0
        i_jump = int(np.argmax(metric))

        self.direction = 1
        i_shock = i_jump + 1
        i_behind = max(1, i_shock - 2)
        i_ahead = min(n - 2, i_shock + 2)

        self.rho2, self.u2, self.p2 = rho[i_behind], u[i_behind], p[i_behind]
        self.rho1, self.u1, self.p1 = rho[i_ahead], u[i_ahead], p[i_ahead]
        self.i_shock = i_shock
        self.x_shock = x_nodes[i_shock]
        self.peak_post_shock_p = max(self.peak_post_shock_p, self.p2)

        d_rho = self.rho2 - self.rho1
        if abs(d_rho) > 1e-6:
            self.u_shock = (self.rho2 * self.u2 - self.rho1 * self.u1) / d_rho
        else:
            self.u_shock = self.u2

        self._detect_contact_surface(rho, x_nodes, n)
        self.shock_detected = True

        # Check for right end-wall reflection
        if self.enable_reflection:
            near_right_wall = (i_shock >= n - 4) or (x_nodes[i_shock] >= x_nodes[-2] - 3.0 * (x_nodes[-2] - x_nodes[-3]))
            p_right_wall = p[-2]
            surge = p_right_wall >= self.reflection_threshold * self.peak_post_shock_p
            if near_right_wall and surge:
                self.wave_mode = 'reflected_left'
                self.x_wall = x_nodes[-2]
                self.direction = -1

        return True

    def _detect_reflected_left(self, rho, u, p, dp, compression, p_jump, x_nodes, n):
        """Track reflected shock traveling to the left into State 2."""
        # For left-running reflected shock:
        # Pressure is higher on the right (p5 > p2 => dp = p[i+1] - p[i] > 0)
        # Velocity drops from u2 > 0 (left) to u5 approx 0 (right) => compression > 0
        valid = (p_jump > self.pressure_threshold) & (compression > 0.0) & (dp > 0.0)

        # Restrict search between contact surface and the right end-wall
        search_min = max(1, self.i_contact)
        search_max = n - 2

        mask = np.zeros(len(dp), dtype=bool)
        mask[search_min:search_max] = True
        valid_refl = valid & mask

        if not np.any(valid_refl):
            # If temporarily inside boundary layer smearing, keep current shock location
            self.shock_detected = True
            return True

        metric = compression * p_jump
        metric[~valid_refl] = -1.0
        i_jump = int(np.argmax(metric))

        self.direction = -1
        i_shock = i_jump
        i_ahead = max(1, i_shock - 2)      # Oncoming State 2 (left)
        i_behind = min(n - 2, i_shock + 2)  # Stagnant State 5 (right)

        self.rho2, self.u2, self.p2 = rho[i_ahead], u[i_ahead], p[i_ahead]
        self.rho5, self.u5, self.p5 = rho[i_behind], u[i_behind], p[i_behind]
        self.i_shock = i_shock
        self.x_shock = x_nodes[i_shock]

        # Instantaneous reflected shock speed W_R < 0 via Rankine-Hugoniot
        d_rho = self.rho5 - self.rho2
        if abs(d_rho) > 1e-6:
            self.u_shock = (self.rho5 * self.u5 - self.rho2 * self.u2) / d_rho
        else:
            self.u_shock = -abs(self.u2)

        self._detect_contact_surface(rho, x_nodes, n)
        self.shock_detected = True
        return True

    def _detect_incident_left(self, rho, u, p, dp, compression, p_jump, x_nodes, n):
        """Track incident left-running shock and monitor for left end-wall collision."""
        # Check if left wall reflection already occurred
        if self.enable_reflection and self.peak_post_shock_p > 0.0:
            p_left_wall = np.max(p[:4])
            if p_left_wall >= self.reflection_threshold * self.peak_post_shock_p:
                self.wave_mode = 'reflected_right'
                self.x_wall = x_nodes[1]
                self.direction = 1
                return self._detect_reflected_right(rho, u, p, dp, compression, p_jump, x_nodes, n)

        valid = (p_jump > self.pressure_threshold) & (compression > 0.0) & (dp > 0.0)

        if not np.any(valid):
            self.shock_detected = False
            return False

        metric = compression * p_jump
        metric[~valid] = -1.0
        i_jump = int(np.argmax(metric))

        self.direction = -1
        i_shock = i_jump
        i_behind = min(n - 2, i_shock + 2)
        i_ahead = max(1, i_shock - 2)

        self.rho2, self.u2, self.p2 = rho[i_behind], u[i_behind], p[i_behind]
        self.rho1, self.u1, self.p1 = rho[i_ahead], u[i_ahead], p[i_ahead]
        self.i_shock = i_shock
        self.x_shock = x_nodes[i_shock]
        self.peak_post_shock_p = max(self.peak_post_shock_p, self.p2)

        d_rho = self.rho2 - self.rho1
        if abs(d_rho) > 1e-6:
            self.u_shock = (self.rho2 * self.u2 - self.rho1 * self.u1) / d_rho
        else:
            self.u_shock = self.u2

        self._detect_contact_surface(rho, x_nodes, n)
        self.shock_detected = True

        # Check for left end-wall reflection
        if self.enable_reflection:
            near_left_wall = (i_shock <= 3) or (x_nodes[i_shock] <= x_nodes[1] + 3.0 * (x_nodes[2] - x_nodes[1]))
            p_left_wall = p[1]
            surge = p_left_wall >= self.reflection_threshold * self.peak_post_shock_p
            if near_left_wall and surge:
                self.wave_mode = 'reflected_right'
                self.x_wall = x_nodes[1]
                self.direction = 1

        return True

    def _detect_reflected_right(self, rho, u, p, dp, compression, p_jump, x_nodes, n):
        """Track reflected shock traveling to the right after reflecting off left wall."""
        valid = (p_jump > self.pressure_threshold) & (compression > 0.0) & (dp < 0.0)

        search_min = 1
        search_max = min(n - 2, max(self.i_contact, 2))

        mask = np.zeros(len(dp), dtype=bool)
        mask[search_min:search_max] = True
        valid_refl = valid & mask

        if not np.any(valid_refl):
            self.shock_detected = True
            return True

        metric = compression * p_jump
        metric[~valid_refl] = -1.0
        i_jump = int(np.argmax(metric))

        self.direction = 1
        i_shock = i_jump + 1
        i_ahead = min(n - 2, i_shock + 2)
        i_behind = max(1, i_shock - 2)

        self.rho2, self.u2, self.p2 = rho[i_ahead], u[i_ahead], p[i_ahead]
        self.rho5, self.u5, self.p5 = rho[i_behind], u[i_behind], p[i_behind]
        self.i_shock = i_shock
        self.x_shock = x_nodes[i_shock]

        d_rho = self.rho5 - self.rho2
        if abs(d_rho) > 1e-6:
            self.u_shock = (self.rho5 * self.u5 - self.rho2 * self.u2) / d_rho
        else:
            self.u_shock = abs(self.u2)

        self._detect_contact_surface(rho, x_nodes, n)
        self.shock_detected = True
        return True

    def _detect_contact_surface(self, rho, x_nodes, n):
        """Identify contact surface separating test gas from driver gas."""
        if self.wave_mode in ['incident_right', 'reflected_left']:
            search_end = max(2, self.i_shock - 2)
            if search_end > 2:
                drho = np.abs(np.diff(rho[1:search_end]))
                if len(drho) > 0:
                    i_contact = 1 + int(np.argmax(drho))
                    self.i_contact = i_contact
                    self.x_contact = x_nodes[i_contact]
            else:
                self.i_contact = 1
                self.x_contact = x_nodes[1]
        else:
            search_start = min(n - 3, self.i_shock + 2)
            if search_start < n - 2:
                drho = np.abs(np.diff(rho[search_start:-1]))
                if len(drho) > 0:
                    i_contact = search_start + int(np.argmax(drho))
                    self.i_contact = i_contact
                    self.x_contact = x_nodes[i_contact]
            else:
                self.i_contact = n - 2
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
    Handles shock tracking, distance calculation, wave state machine (incident/reflected),
    and leading-edge finite-volume singularity cell averaging.
    """
    def __init__(self, pressure_threshold=0.05, max_cf=0.1, driver_cf=0.003,
                 enable_reflection=True, reflection_threshold=1.15):
        self.tracker = ShockTracker(
            pressure_threshold=pressure_threshold,
            enable_reflection=enable_reflection,
            reflection_threshold=reflection_threshold
        )
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
            return np.full(n, self.driver_cf)

        wave_mode = self.tracker.wave_mode
        rho2 = self.tracker.rho2
        u2 = self.tracker.u2
        p2 = self.tracker.p2
        u_shock = self.tracker.u_shock
        x_shock = self.tracker.x_shock
        x_contact = self.tracker.x_contact
        x_wall = self.tracker.x_wall

        # Dynamic viscosity of State 2
        try:
            mu2 = fluid.computeViscosity_p_rho(p2, rho2)
        except Exception:
            mu2 = 1.8e-5

        # Dynamic viscosity of State 5 (if in reflected mode)
        if wave_mode in ['reflected_left', 'reflected_right']:
            try:
                mu5 = fluid.computeViscosity_p_rho(self.tracker.p5, self.tracker.rho5)
            except Exception:
                mu5 = mu2
            rho5 = self.tracker.rho5
            u5 = self.tracker.u5

        for i in range(n):
            xi = x_nodes[i]
            dxi = dx[i]

            if wave_mode == 'incident_right':
                if xi > x_shock:
                    cf[i] = 0.0
                elif xi >= x_contact:
                    dist = x_shock - xi
                    cf[i] = self._evaluate_cf_local(dist, dxi, rho2, u2, u_shock, mu2)
                else:
                    cf[i] = self.driver_cf

            elif wave_mode == 'reflected_left':
                if xi < x_contact:
                    cf[i] = self.driver_cf
                elif xi < x_shock:
                    # Oncoming State 2 test gas, moving toward wall at u2 > 0
                    dist = max(dxi, x_wall - xi)
                    cf[i] = self._evaluate_cf_local(dist, dxi, rho2, u2, abs(u_shock - u2), mu2)
                elif xi <= x_wall:
                    # Post-reflected State 5 gas behind the left-running reflected shock
                    dist = xi - x_shock
                    cf[i] = self._evaluate_cf_local(dist, dxi, rho5, u5, u_shock, mu5)
                else:
                    cf[i] = 0.0

            elif wave_mode == 'incident_left':
                if xi < x_shock:
                    cf[i] = 0.0
                elif xi <= x_contact:
                    dist = xi - x_shock
                    cf[i] = self._evaluate_cf_local(dist, dxi, rho2, u2, u_shock, mu2)
                else:
                    cf[i] = self.driver_cf

            elif wave_mode == 'reflected_right':
                if xi > x_contact:
                    cf[i] = self.driver_cf
                elif xi > x_shock:
                    # Oncoming State 2 test gas, moving left toward wall at u2 < 0
                    dist = max(dxi, xi - x_wall)
                    cf[i] = self._evaluate_cf_local(dist, dxi, rho2, u2, abs(u_shock - u2), mu2)
                elif xi >= x_wall:
                    # Post-reflected State 5 gas behind right-running reflected shock
                    dist = x_shock - xi
                    cf[i] = self._evaluate_cf_local(dist, dxi, rho5, u5, u_shock, mu5)
                else:
                    cf[i] = 0.0

        cf = np.clip(cf, 0.0, self.max_cf)
        return cf

    def _evaluate_cf_local(self, dist, dxi, rho_eval, u_eval, u_shock_eval, mu_eval):
        """Override in subclasses to evaluate specific Cf correlations."""
        raise NotImplementedError


class MirelsLaminarFriction(MirelsFrictionBase):
    """
    Mirels Laminar Boundary Layer Model (NACA TN 3401 / Blasius approximation):
    Cf(x) = 0.664 / sqrt(Re_x)
    Singularity regularization over shock cell [0, dx]:
    bar_Cf = 2.0 * Cf(dx)
    """
    def _evaluate_cf_local(self, dist, dxi, rho_eval, u_eval, u_shock_eval, mu_eval):
        n_exp = 0.5
        if dist < dxi:
            re_dx = self._compute_reynolds(dxi, rho_eval, u_eval, u_shock_eval, mu_eval)
            re_safe = max(re_dx, 1e-4)
            cf_edge = 0.664 / np.sqrt(re_safe)
            return self._regularize_singularity(cf_edge, n_exp)
        else:
            re_x = self._compute_reynolds(dist, rho_eval, u_eval, u_shock_eval, mu_eval)
            re_safe = max(re_x, 1e-4)
            return 0.664 / np.sqrt(re_safe)


class MirelsTurbulentFriction(MirelsFrictionBase):
    """
    Mirels Turbulent Boundary Layer Model (NACA TN 3712 / AIAA J. 1964):
    Cf(x) = 0.0592 / Re_x^0.2
    Singularity regularization over shock cell [0, dx]:
    bar_Cf = 1.25 * Cf(dx)
    """
    def _evaluate_cf_local(self, dist, dxi, rho_eval, u_eval, u_shock_eval, mu_eval):
        n_exp = 0.2
        if dist < dxi:
            re_dx = self._compute_reynolds(dxi, rho_eval, u_eval, u_shock_eval, mu_eval)
            re_safe = max(re_dx, 1e-4)
            cf_edge = 0.0592 / (re_safe ** 0.2)
            return self._regularize_singularity(cf_edge, n_exp)
        else:
            re_x = self._compute_reynolds(dist, rho_eval, u_eval, u_shock_eval, mu_eval)
            re_safe = max(re_x, 1e-4)
            return 0.0592 / (re_safe ** 0.2)


class MirelsTransitionalFriction(MirelsFrictionBase):
    """
    Mirels Transitional Boundary Layer Model:
    Switches from laminar (0.664 / sqrt(Re_x)) to turbulent (0.0592 / Re_x^0.2)
    when Re_x >= re_transition.
    """
    def __init__(self, re_transition=1.0e6, pressure_threshold=0.05, max_cf=0.1, driver_cf=0.003,
                 enable_reflection=True, reflection_threshold=1.15):
        super().__init__(
            pressure_threshold=pressure_threshold,
            max_cf=max_cf,
            driver_cf=driver_cf,
            enable_reflection=enable_reflection,
            reflection_threshold=reflection_threshold
        )
        self.re_transition = re_transition

    def _evaluate_cf_local(self, dist, dxi, rho_eval, u_eval, u_shock_eval, mu_eval):
        if dist < dxi:
            re_dx = self._compute_reynolds(dxi, rho_eval, u_eval, u_shock_eval, mu_eval)
            re_safe = max(re_dx, 1e-4)
            cf_edge = 0.664 / np.sqrt(re_safe)
            return self._regularize_singularity(cf_edge, 0.5)
        else:
            re_x = self._compute_reynolds(dist, rho_eval, u_eval, u_shock_eval, mu_eval)
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
    enable_refl = getattr(config, 'isShockReflectionTrackingActive', lambda: True)()
    refl_thresh = getattr(config, 'getWallReflectionThreshold', lambda: 1.15)()

    if model_name in ['mirels_laminar', 'laminar']:
        return MirelsLaminarFriction(
            pressure_threshold=p_threshold,
            max_cf=max_cf,
            driver_cf=driver_cf,
            enable_reflection=enable_refl,
            reflection_threshold=refl_thresh
        )
    elif model_name in ['mirels_turbulent', 'turbulent']:
        return MirelsTurbulentFriction(
            pressure_threshold=p_threshold,
            max_cf=max_cf,
            driver_cf=driver_cf,
            enable_reflection=enable_refl,
            reflection_threshold=refl_thresh
        )
    elif model_name in ['mirels_transitional', 'transitional']:
        re_tr = config.getFrictionTransitionReynolds()
        return MirelsTransitionalFriction(
            re_transition=re_tr,
            pressure_threshold=p_threshold,
            max_cf=max_cf,
            driver_cf=driver_cf,
            enable_reflection=enable_refl,
            reflection_threshold=refl_thresh
        )
    else:
        cf_const = config.getFrictionCoefficient()
        return ConstantFriction(cf_constant=cf_const)
