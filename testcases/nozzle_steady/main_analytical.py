import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


# ============================================================
# QUASI-1D COMPRESSIBLE NOZZLE SOLVER
# ============================================================
#
# Perfect gas, steady, quasi-1D, adiabatic flow.
#
# Nozzle geometry:
#
#     A(x) = alpha * (x^2 - x) + delta
#
#     x in [0, 1]
#
# Inlet total conditions:
#
#     Pt_in = 101325 Pa
#     Tt_in = 288.15 K
#
# Outlet pressures:
#
#     45, 70, 90, 97 kPa
#
# The solver handles:
#
#     1. Fully subsonic flow
#     2. Choked flow, supersonic to exit
#     3. Choked flow with an internal normal shock
#
# ============================================================


# ============================================================
# GAS PROPERTIES
# ============================================================

gamma = 1.4
R = 287.05                         # [J/(kg K)]


# ============================================================
# INLET TOTAL CONDITIONS
# ============================================================

Pt_in = 101325.0                   # [Pa]
Tt_in = 288.15                     # [K]
p_out_cases = [45e3, 75e3, 90e3, 95e3, 97e3]

# ============================================================
# NOZZLE GEOMETRY
# ============================================================

alpha = 1.364e-3
delta = 7.211e-4


def area(x):
    """
    Nozzle area distribution:

        A(x) = alpha * (x^2 - x) + delta
    """

    return alpha * (x**2 - x) + delta


# Computational grid
x = np.linspace(0.0, 1.0, 2001)

A = area(x)

# Throat location
x_throat = 0.5

# Throat and exit areas
A_throat = area(x_throat)
A_exit = area(1.0)

# Exit-to-throat area ratio
A_ratio_exit = A_exit / A_throat


# ============================================================
# ISENTROPIC RELATIONS
# ============================================================

def area_mach(M):
    """
    Isentropic area-Mach relation:

        A/A* =
        1/M *
        [
            2/(gamma+1)
            *
            (1 + (gamma-1)/2 M^2)
        ]
        ^[(gamma+1)/(2(gamma-1))]
    """

    return (
        1.0 / M
        *
        (
            2.0 / (gamma + 1.0)
            *
            (
                1.0
                + 0.5 * (gamma - 1.0) * M**2
            )
        )
        ** (
            (gamma + 1.0)
            / (2.0 * (gamma - 1.0))
        )
    )


def pressure_from_mach(M, Pt):
    """
    Static pressure from Mach number and total pressure.
    """

    return Pt * (
        1.0
        + 0.5 * (gamma - 1.0) * M**2
    ) ** (
        -gamma / (gamma - 1.0)
    )


def temperature_from_mach(M, Tt):
    """
    Static temperature from Mach number and total temperature.
    """

    return Tt / (
        1.0
        + 0.5 * (gamma - 1.0) * M**2
    )


def density_from_p_T(p, T):
    """
    Perfect-gas density.
    """

    return p / (R * T)


def velocity_from_M_T(M, T):
    """
    Flow velocity.
    """

    return M * np.sqrt(gamma * R * T)


def mach_from_pressure(p, Pt):
    """
    Mach number from static and total pressure.

        p/Pt =
        [1 + (gamma-1)/2 M^2]^[-gamma/(gamma-1)]
    """

    value = (
        (Pt / p)
        ** ((gamma - 1.0) / gamma)
        - 1.0
    )

    return np.sqrt(
        2.0 / (gamma - 1.0)
        * value
    )


# ============================================================
# SOLVE AREA-MACH RELATION
# ============================================================

def mach_from_area_ratio(A_ratio, branch):
    """
    Obtain Mach number from A/A*.

    branch = 'subsonic' or 'supersonic'
    """

    if A_ratio < 1.0 - 1e-10:

        raise ValueError(
            f"A/A* = {A_ratio:.8e} < 1"
        )

    if np.isclose(
        A_ratio,
        1.0,
        atol=1e-10
    ):

        return 1.0

    def residual(M):

        return area_mach(M) - A_ratio

    if branch == "subsonic":

        return brentq(
            residual,
            1e-8,
            1.0 - 1e-10
        )

    elif branch == "supersonic":

        return brentq(
            residual,
            1.0 + 1e-10,
            100.0
        )

    else:

        raise ValueError(
            "branch must be 'subsonic' or 'supersonic'"
        )


# ============================================================
# CRITICAL CONDITIONS
# ============================================================

def critical_pressure():

    return Pt_in * (
        2.0 / (gamma + 1.0)
    ) ** (
        gamma / (gamma - 1.0)
    )


def choked_mass_flow_rate():

    return (
        A_throat
        * Pt_in
        / np.sqrt(Tt_in)
        * np.sqrt(gamma / R)
        * (
            (gamma + 1.0) / 2.0
        )
        ** (
            -(gamma + 1.0)
            / (2.0 * (gamma - 1.0))
        )
    )


# ============================================================
# FULLY SUBSONIC SOLUTION
# ============================================================

def solve_fully_subsonic(p_out):

    # --------------------------------------------------------
    # Exit Mach number
    # --------------------------------------------------------

    M_exit = mach_from_pressure(
        p_out,
        Pt_in
    )

    if M_exit >= 1.0:

        raise ValueError(
            "Exit Mach number is >= 1."
        )

    # --------------------------------------------------------
    # Effective critical area
    # --------------------------------------------------------

    A_star = (
        A_exit
        / area_mach(M_exit)
    )

    # Check throat condition
    throat_ratio = A_throat / A_star

    if throat_ratio < 1.0:

        raise ValueError(
            "The nozzle is choked."
        )

    # --------------------------------------------------------
    # Mach distribution
    # --------------------------------------------------------

    M = np.zeros_like(x)

    for i in range(len(x)):

        M[i] = mach_from_area_ratio(
            A[i] / A_star,
            "subsonic"
        )

    # --------------------------------------------------------
    # Thermodynamic properties
    # --------------------------------------------------------

    p = pressure_from_mach(
        M,
        Pt_in
    )

    T = temperature_from_mach(
        M,
        Tt_in
    )

    rho = density_from_p_T(
        p,
        T
    )

    u = velocity_from_M_T(
        M,
        T
    )

    mdot = rho * u * A

    return {
        "regime": "Fully subsonic",
        "x": x,
        "A": A,
        "M": M,
        "p": p,
        "T": T,
        "rho": rho,
        "u": u,
        "mdot": mdot,
        "M_exit": M[-1],
        "p_exit": p[-1],
        "A_star": A_star
    }


# ============================================================
# CHOKED ISENTROPIC SOLUTION
# ============================================================

def solve_choked_isentropic():

    M = np.zeros_like(x)

    for i in range(len(x)):

        A_ratio = A[i] / A_throat

        if x[i] < x_throat:

            M[i] = mach_from_area_ratio(
                A_ratio,
                "subsonic"
            )

        elif np.isclose(
            x[i],
            x_throat
        ):

            M[i] = 1.0

        else:

            M[i] = mach_from_area_ratio(
                A_ratio,
                "supersonic"
            )

    p = pressure_from_mach(
        M,
        Pt_in
    )

    T = temperature_from_mach(
        M,
        Tt_in
    )

    rho = density_from_p_T(
        p,
        T
    )

    u = velocity_from_M_T(
        M,
        T
    )

    mdot = rho * u * A

    return {
        "regime": "Choked, supersonic to exit",
        "x": x,
        "A": A,
        "M": M,
        "p": p,
        "T": T,
        "rho": rho,
        "u": u,
        "mdot": mdot,
        "M_exit": M[-1],
        "p_exit": p[-1],
        "A_star": A_throat
    }


# ============================================================
# NORMAL SHOCK RELATIONS
# ============================================================

def normal_shock_M2(M1):

    return np.sqrt(
        (
            1.0
            + 0.5 * (gamma - 1.0) * M1**2
        )
        /
        (
            gamma * M1**2
            - 0.5 * (gamma - 1.0)
        )
    )


def normal_shock_p2_p1(M1):

    return (
        1.0
        + 2.0 * gamma / (gamma + 1.0)
        * (M1**2 - 1.0)
    )


def normal_shock_rho2_rho1(M1):

    return (
        (gamma + 1.0) * M1**2
        /
        (
            2.0
            + (gamma - 1.0) * M1**2
        )
    )


def normal_shock_T2_T1(M1):

    return (
        normal_shock_p2_p1(M1)
        /
        normal_shock_rho2_rho1(M1)
    )


def normal_shock_Pt2_Pt1(M1):

    """
    Total pressure ratio across a normal shock.
    """

    return (
        (
            (gamma + 1.0) * M1**2
            /
            (
                (gamma - 1.0) * M1**2
                + 2.0
            )
        )
        ** (gamma / (gamma - 1.0))
        *
        (
            (gamma + 1.0)
            /
            (
                2.0 * gamma * M1**2
                - (gamma - 1.0)
            )
        )
        ** (1.0 / (gamma - 1.0))
    )


# ============================================================
# EXIT PRESSURE FOR A GIVEN SHOCK LOCATION
# ============================================================

def exit_pressure_for_shock(x_shock):
    """
    Calculate the exit pressure produced by a normal shock
    located at x_shock.
    """

    # --------------------------------------------------------
    # Area at shock
    # --------------------------------------------------------

    A_shock = area(x_shock)

    # --------------------------------------------------------
    # Upstream Mach number
    #
    # Upstream flow is choked and therefore:
    #
    #     A* = A_throat
    #
    # --------------------------------------------------------

    M1 = mach_from_area_ratio(
        A_shock / A_throat,
        "supersonic"
    )

    # --------------------------------------------------------
    # Normal shock
    # --------------------------------------------------------

    M2 = normal_shock_M2(M1)

    Pt2_Pt1 = normal_shock_Pt2_Pt1(M1)

    Pt2 = Pt_in * Pt2_Pt1

    # --------------------------------------------------------
    # New downstream critical area
    #
    # The downstream flow has Mach M2 at A_shock.
    #
    # Therefore:
    #
    #     A_shock / A*_2 = f(M2)
    #
    # --------------------------------------------------------

    A_star_2 = (
        A_shock
        / area_mach(M2)
    )

    # --------------------------------------------------------
    # Exit Mach number
    # --------------------------------------------------------

    A_exit_ratio = (
        A_exit
        / A_star_2
    )

    M_exit = mach_from_area_ratio(
        A_exit_ratio,
        "subsonic"
    )

    # --------------------------------------------------------
    # Exit pressure
    # --------------------------------------------------------

    p_exit = pressure_from_mach(
        M_exit,
        Pt2
    )

    return {
        "x_shock": x_shock,
        "A_shock": A_shock,
        "M1": M1,
        "M2": M2,
        "Pt2": Pt2,
        "A_star_2": A_star_2,
        "M_exit": M_exit,
        "p_exit": p_exit
    }


# ============================================================
# COMPLETE SHOCK SOLUTION
# ============================================================

def solve_with_shock(p_out):
    """
    Find shock location satisfying:

        p_exit = p_out.
    """

    # Shock must be downstream of the throat.
    x_min = x_throat + 1e-6
    x_max = 1.0 - 1e-6

    # --------------------------------------------------------
    # Residual
    # --------------------------------------------------------

    def residual(x_shock):

        return (
            exit_pressure_for_shock(
                x_shock
            )["p_exit"]
            - p_out
        )

    # --------------------------------------------------------
    # Check the pressure range
    # --------------------------------------------------------

    p_min = exit_pressure_for_shock(
        x_min
    )["p_exit"]

    p_max = exit_pressure_for_shock(
        x_max
    )["p_exit"]

    print()
    print(
        f"Shock pressure range: "
        f"{min(p_min, p_max)/1e3:.6f} "
        f"to "
        f"{max(p_min, p_max)/1e3:.6f} kPa"
    )

    # --------------------------------------------------------
    # Root solve
    # --------------------------------------------------------

    if residual(x_min) * residual(x_max) > 0.0:

        raise ValueError(
            f"No internal shock solution for "
            f"p_out = {p_out/1e3:.3f} kPa.\n"
            f"Allowed pressure range is approximately "
            f"{min(p_min,p_max)/1e3:.3f}--"
            f"{max(p_min,p_max)/1e3:.3f} kPa."
        )

    x_shock = brentq(
        residual,
        x_min,
        x_max,
        xtol=1e-12
    )

    shock = exit_pressure_for_shock(
        x_shock
    )

    # --------------------------------------------------------
    # Extract shock quantities
    # --------------------------------------------------------

    A_shock = shock["A_shock"]
    M1 = shock["M1"]
    M2 = shock["M2"]

    Pt2 = shock["Pt2"]
    A_star_2 = shock["A_star_2"]

    # --------------------------------------------------------
    # Static conditions immediately before shock
    # --------------------------------------------------------

    p1 = pressure_from_mach(
        M1,
        Pt_in
    )

    T1 = temperature_from_mach(
        M1,
        Tt_in
    )

    # --------------------------------------------------------
    # Static conditions immediately after shock
    # --------------------------------------------------------

    p2 = (
        p1
        * normal_shock_p2_p1(M1)
    )

    T2 = (
        T1
        * normal_shock_T2_T1(M1)
    )

    # --------------------------------------------------------
    # Build solution
    # --------------------------------------------------------

    M = np.zeros_like(x)
    p = np.zeros_like(x)
    T = np.zeros_like(x)

    for i in range(len(x)):

        # ----------------------------------------------------
        # Converging section
        # ----------------------------------------------------

        if x[i] < x_throat:

            M[i] = mach_from_area_ratio(
                A[i] / A_throat,
                "subsonic"
            )

            p[i] = pressure_from_mach(
                M[i],
                Pt_in
            )

            T[i] = temperature_from_mach(
                M[i],
                Tt_in
            )

        # ----------------------------------------------------
        # Throat
        # ----------------------------------------------------

        elif x[i] < x_shock:

            M[i] = mach_from_area_ratio(
                A[i] / A_throat,
                "supersonic"
            )

            p[i] = pressure_from_mach(
                M[i],
                Pt_in
            )

            T[i] = temperature_from_mach(
                M[i],
                Tt_in
            )

        # ----------------------------------------------------
        # Shock
        # ----------------------------------------------------

        elif np.isclose(
            x[i],
            x_shock,
            atol=0.5 * (x[1] - x[0])
        ):

            M[i] = M2
            p[i] = p2
            T[i] = T2

        # ----------------------------------------------------
        # Downstream of shock
        # ----------------------------------------------------

        else:

            M[i] = mach_from_area_ratio(
                A[i] / A_star_2,
                "subsonic"
            )

            p[i] = pressure_from_mach(
                M[i],
                Pt2
            )

            T[i] = temperature_from_mach(
                M[i],
                Tt_in
            )

    # --------------------------------------------------------
    # Other quantities
    # --------------------------------------------------------

    rho = density_from_p_T(
        p,
        T
    )

    u = velocity_from_M_T(
        M,
        T
    )

    mdot = rho * u * A

    return {
        "regime": "Choked with internal normal shock",
        "x": x,
        "A": A,
        "M": M,
        "p": p,
        "T": T,
        "rho": rho,
        "u": u,
        "mdot": mdot,
        "x_shock": x_shock,
        "A_shock": A_shock,
        "M1_shock": M1,
        "M2_shock": M2,
        "p1_shock": p1,
        "p2_shock": p2,
        "Pt2": Pt2,
        "A_star_2": A_star_2,
        "M_exit": M[-1],
        "p_exit": p[-1]
    }


# ============================================================
# AUTOMATIC SOLVER
# ============================================================

def solve_nozzle(p_out):

    # ========================================================
    # Characteristic pressures
    # ========================================================

    p_critical = critical_pressure()

    # Fully choked isentropic solution
    choked = solve_choked_isentropic()

    # Pressure obtained at the exit if the flow remains
    # isentropic and supersonic throughout the nozzle
    p_exit_iso = choked["p_exit"]

    # --------------------------------------------------------
    # Pressure corresponding to a normal shock exactly at
    # the nozzle exit.
    #
    # This is the lowest back pressure for which the shock
    # can still be located inside the nozzle.
    # --------------------------------------------------------

    shock_at_exit = exit_pressure_for_shock(
        1.0 - 1e-8
    )

    p_shock_exit = shock_at_exit["p_exit"]

    # --------------------------------------------------------
    # Pressure corresponding to a normal shock immediately
    # downstream of the throat.
    #
    # This is the highest back pressure for which an internal
    # shock can exist.
    # --------------------------------------------------------

    shock_at_throat = exit_pressure_for_shock(
        x_throat + 1e-8
    )

    p_shock_throat = shock_at_throat["p_exit"]

    # ========================================================
    # REGIME 1: Fully subsonic
    # ========================================================

    if p_out >= p_shock_throat:

        return solve_fully_subsonic(
            p_out
        )

    # ========================================================
    # REGIME 2: Choked + internal normal shock
    # ========================================================

    elif p_out >= p_shock_exit:

        return solve_with_shock(
            p_out
        )

    # ========================================================
    # REGIME 3: Choked + supersonic to exit
    #
    # The shock has moved outside the nozzle.
    # ========================================================

    else:

        choked["regime"] = (
            "Choked, supersonic to exit"
        )

        return choked


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    # --------------------------------------------------------
    # Basic geometry
    # --------------------------------------------------------

    print()
    print("=" * 75)
    print("QUASI-1D COMPRESSIBLE NOZZLE SOLVER")
    print("=" * 75)

    print()
    print("NOZZLE GEOMETRY")
    print("-" * 75)

    print(
        f"alpha             = {alpha:.6e}"
    )

    print(
        f"delta             = {delta:.6e}"
    )

    print(
        f"x_throat          = {x_throat:.6f}"
    )

    print(
        f"A_throat          = {A_throat:.8e}"
    )

    print(
        f"A_exit            = {A_exit:.8e}"
    )

    print(
        f"A_exit/A_throat   = {A_ratio_exit:.8f}"
    )

    # --------------------------------------------------------
    # Critical conditions
    # --------------------------------------------------------

    print()
    print("CRITICAL CONDITIONS")
    print("-" * 75)

    print(
        f"p_critical        = "
        f"{critical_pressure()/1e3:.6f} kPa"
    )

    print(
        f"mdot_choked       = "
        f"{choked_mass_flow_rate():.8e} kg/s"
    )

    # --------------------------------------------------------
    # Choked isentropic exit
    # --------------------------------------------------------

    choked = solve_choked_isentropic()

    print(
        f"isentropic M_exit = "
        f"{choked['M_exit']:.6f}"
    )

    print(
        f"isentropic p_exit = "
        f"{choked['p_exit']/1e3:.6f} kPa"
    )

    # --------------------------------------------------------
    # Solve cases
    # --------------------------------------------------------

    solutions = {}
    for p_out in p_out_cases:

        solutions[p_out] = solve_nozzle(
            p_out
        )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 75)
    print("OPERATING CONDITIONS")
    print("=" * 75)

    print(
        f"{'p_out [kPa]':>14}"
        f"{'Regime':>40}"
        f"{'M_exit':>14}"
        f"{'p_exit [kPa]':>16}"
    )

    print("-" * 84)

    for p_out in p_out_cases:

        sol = solutions[p_out]

        print(
            f"{p_out/1e3:14.3f}"
            f"{sol['regime']:>40}"
            f"{sol['M_exit']:14.6f}"
            f"{sol['p_exit']/1e3:16.6f}"
        )

    # --------------------------------------------------------
    # Shock information
    # --------------------------------------------------------

    print()
    print("=" * 75)
    print("SHOCK INFORMATION")
    print("=" * 75)

    for p_out in p_out_cases:

        sol = solutions[p_out]

        if "x_shock" in sol:

            print()
            print(
                f"p_out = "
                f"{p_out/1e3:.3f} kPa"
            )

            print(
                f"x_shock       = "
                f"{sol['x_shock']:.8f}"
            )

            print(
                f"A_shock       = "
                f"{sol['A_shock']:.8e}"
            )

            print(
                f"M1             = "
                f"{sol['M1_shock']:.6f}"
            )

            print(
                f"M2             = "
                f"{sol['M2_shock']:.6f}"
            )

            print(
                f"p1             = "
                f"{sol['p1_shock']/1e3:.6f} kPa"
            )

            print(
                f"p2             = "
                f"{sol['p2_shock']/1e3:.6f} kPa"
            )

            print(
                f"Pt2/Pt1         = "
                f"{sol['Pt2']/Pt_in:.8f}"
            )

    # --------------------------------------------------------
    # Mass flow check
    # --------------------------------------------------------

    print()
    print("=" * 75)
    print("MASS FLOW RATE CONSISTENCY")
    print("=" * 75)

    for p_out in p_out_cases:

        sol = solutions[p_out]

        mdot = sol["mdot"]

        mean_mdot = np.mean(mdot)

        variation = (
            np.max(mdot)
            - np.min(mdot)
        ) / mean_mdot

        print(
            f"p_out = {p_out/1e3:6.1f} kPa : "
            f"mean mdot = {mean_mdot:.8e} kg/s, "
            f"relative variation = {variation:.3e}"
        )

    # ========================================================
    # PLOTS
    # ========================================================

    # --------------------------------------------------------
    # Mach number
    # --------------------------------------------------------

    plt.figure(
        figsize=(8, 5)
    )

    for p_out in p_out_cases:

        sol = solutions[p_out]

        plt.plot(
            sol["x"],
            sol["M"],
            label=f"$p_{{out}}={p_out/1e3:.0f}$ kPa"
        )

        if "x_shock" in sol:

            plt.axvline(
                sol["x_shock"],
                linestyle=":",
                linewidth=1.0
            )

    plt.axhline(
        1.0,
        linestyle="--",
        linewidth=1.0
    )

    plt.axvline(
        x_throat,
        linestyle=":",
        linewidth=1.0
    )

    plt.xlabel("$x$")
    plt.ylabel("$M$")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # --------------------------------------------------------
    # Static pressure
    # --------------------------------------------------------

    plt.figure(
        figsize=(8, 5)
    )

    for p_out in p_out_cases:

        sol = solutions[p_out]

        plt.plot(
            sol["x"],
            sol["p"] / 1e3,
            label=f"$p_{{out}}={p_out/1e3:.0f}$ kPa"
        )

        if "x_shock" in sol:

            plt.axvline(
                sol["x_shock"],
                linestyle=":",
                linewidth=1.0
            )

    plt.axvline(
        x_throat,
        linestyle="--",
        linewidth=1.0
    )

    plt.xlabel("$x$")
    plt.ylabel("$p$ [kPa]")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    
    # write the solution in a single pickle file, containing a dict
    import pickle
    with open('nozzle_solutions_analytical.pkl', 'wb') as f:
        pickle.dump(solutions, f)
    print('solutions written to nozzle_solutions_analytical.pkl')

    plt.show()