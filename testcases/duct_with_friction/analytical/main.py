#!/usr/bin/env python3

"""
1-D Fanno-flow solver for a constant-area circular duct.

Assumptions
-----------
    - steady 1-D compressible flow
    - constant cross-sectional area
    - circular duct
    - adiabatic wall (q_w = 0)
    - no shaft work
    - constant Fanning friction coefficient Cf
    - calorically perfect gas
    - subsonic Fanno-flow branch

Inputs
------
    L          : duct length [m]
    A          : cross-sectional area [m^2]
    pt_in      : inlet total pressure [Pa]
    Tt_in      : inlet total temperature [K]
    p_out      : outlet static pressure [Pa]
    Cf         : Fanning friction coefficient [-]

The solution is obtained from the analytical Fanno-flow relations.

For a circular duct:

    D = sqrt(4 A / pi)

and the Fanno parameter is

    4 Cf L / D

The wall shear stress is

    tau_w = 0.5 * rho * u^2 * Cf

For an adiabatic Fanno flow:

    Tt = constant
    ht = constant
    pt decreases
    entropy increases

The present implementation solves the SUBSONIC branch.
Friction drives the flow toward M = 1.
"""


import numpy as np
import matplotlib.pyplot as plt
import pickle

from scipy.optimize import brentq


# =====================================================================
# INPUT DATA
# =====================================================================

L = 10.0                    # duct length [m]
A = 0.01                    # cross-sectional area [m^2]

pt_in = 101.3e3             # inlet total pressure [Pa]
Tt_in = 288.15              # inlet total temperature [K]

p_out = 90.0e3              # outlet static pressure [Pa]

Cfs = [
    0.000,
    0.001,
    0.003,
    0.010,
    0.030,
]


# =====================================================================
# GAS PROPERTIES
# =====================================================================

gamma = 1.4
R = 287.05                  # [J/(kg K)]

cp = gamma * R / (gamma - 1.0)


# =====================================================================
# NUMERICAL PARAMETERS
# =====================================================================

N = 500

XTOL = 1e-12
RTOL = 1e-12


# =====================================================================
# GEOMETRY
# =====================================================================

# Circular duct:
D = np.sqrt(4.0 * A / np.pi)

# Wetted perimeter:
P = np.pi * D

# P/A = 4/D:
P_over_A = P / A


# =====================================================================
# ISENTROPIC RELATIONS
# =====================================================================

def static_temperature(M):
    """
    Static temperature from the constant total temperature.

        T = Tt / [1 + (gamma - 1)/2 M^2]

    In Fanno flow:

        Tt = constant.
    """

    return Tt_in / (
        1.0
        +
        0.5 * (gamma - 1.0) * M**2
    )


def speed_of_sound(T):
    """
    Speed of sound.
    """

    return np.sqrt(gamma * R * T)


def total_pressure_from_static(p, M):
    """
    Total pressure corresponding to a local static state.

        pt = p [1 + (gamma - 1)/2 M^2]^(gamma/(gamma - 1))

    This is the stagnation pressure associated with the local
    thermodynamic state through an isentropic deceleration.

    IMPORTANT:
        In Fanno flow pt is NOT constant.
    """

    return p * (
        1.0
        +
        0.5 * (gamma - 1.0) * M**2
    ) ** (gamma / (gamma - 1.0))


def entropy_from_state(p, T):
    """
    Entropy relative to an arbitrary reference state.

    Only entropy differences are relevant here.

        ds = cp ln(T2/T1) - R ln(p2/p1)

    The reference constants therefore cancel.
    """

    return (
        cp * np.log(T / Tt_in)
        -
        R * np.log(p / pt_in)
    )


# =====================================================================
# FANNO FUNCTION
# =====================================================================

def fanno_function(M):
    """
    Standard Fanno-flow function:

        F(M) =
            (1 - M^2)/(gamma M^2)
            +
            (gamma + 1)/(2 gamma)
            ln[
                (gamma + 1) M^2
                /
                (2 + (gamma - 1) M^2)
            ]

    For the subsonic branch:

        M -> 0  : F -> infinity
        M -> 1  : F -> 0

    Therefore F(M) is monotonically decreasing for 0 < M < 1.
    """

    term1 = (
        1.0 - M**2
    ) / (
        gamma * M**2
    )

    argument = (
        (gamma + 1.0) * M**2
        /
        (
            2.0
            +
            (gamma - 1.0) * M**2
        )
    )

    term2 = (
        (gamma + 1.0)
        /
        (2.0 * gamma)
        *
        np.log(argument)
    )

    return term1 + term2


# =====================================================================
# FANNO STAR RELATIONS
# =====================================================================

def p_over_pstar(M):
    """
    Fanno relation:

        p / p*

            =
            1/M *
            sqrt[
                (gamma + 1)
                /
                (2 + (gamma - 1) M^2)
            ]

    where p* is the static pressure at the sonic Fanno state.
    """

    return (
        1.0 / M
        *
        np.sqrt(
            (gamma + 1.0)
            /
            (
                2.0
                +
                (gamma - 1.0) * M**2
            )
        )
    )


def T_over_Tstar(M):
    """
    Fanno relation:

        T / T*
            =
            [(gamma + 1)/2]
            /
            [1 + (gamma - 1)/2 M^2]

    Equivalent form:

        T/T*
            =
            (gamma + 1)
            /
            [2 + (gamma - 1) M^2]
    """

    return (
        (gamma + 1.0)
        /
        (
            2.0
            +
            (gamma - 1.0) * M**2
        )
    )


def rho_over_rhostar(M):
    """
    Fanno relation:

        rho / rho*
            =
            1/M *
            sqrt[
                (2 + (gamma - 1) M^2)
                /
                (gamma + 1)
            ]
    """

    return (
        1.0 / M
        *
        np.sqrt(
            (
                2.0
                +
                (gamma - 1.0) * M**2
            )
            /
            (gamma + 1.0)
        )
    )


def velocity_over_ustar(M):
    """
    Since mass flow is constant in a constant-area duct:

        rho u = rho* u*

    therefore

        u/u* = rho*/rho.
    """

    return 1.0 / rho_over_rhostar(M)


# =====================================================================
# MACH NUMBER EVOLUTION
# =====================================================================

def mach_out_from_mach_in(M_in, Cf):
    """
    Compute M_out from M_in using

        F(M_in) - F(M_out) = 4 Cf L / D

    on the subsonic branch.
    """

    if Cf == 0.0:
        return M_in

    friction_parameter = (
        4.0 * Cf * L / D
    )

    F_in = fanno_function(M_in)

    F_out = (
        F_in
        -
        friction_parameter
    )

    if F_out < 0.0:

        raise ValueError(
            "The specified inlet Mach number would choke "
            "before reaching the end of the duct."
        )

    if np.isclose(
        F_out,
        0.0,
        atol=1e-14
    ):
        return 1.0

    return brentq(
        lambda M:
            fanno_function(M) - F_out,
        M_in,
        1.0 - 1e-12,
        xtol=XTOL,
        rtol=RTOL,
    )


# =====================================================================
# OUTLET STATE FOR A GIVEN INLET MACH NUMBER
# =====================================================================

def outlet_pressure(M_in, Cf):
    """
    Calculate the outlet static pressure corresponding to a trial
    inlet Mach number.

    The procedure is:

        1. Calculate p_in from pt_in and M_in.
        2. Determine p* from the Fanno p/p* relation.
        3. Determine M_out from the Fanno length relation.
        4. Calculate p_out from p* and M_out.

    This is the key correction relative to the previous solver.

    The local static pressure is NOT calculated using the inlet
    total pressure throughout the duct.
    """

    # -------------------------------------------------------------
    # Inlet static pressure
    # -------------------------------------------------------------

    p_in = pt_in / (
        1.0
        +
        0.5 * (gamma - 1.0) * M_in**2
    ) ** (
        gamma / (gamma - 1.0)
    )

    # -------------------------------------------------------------
    # Outlet Mach number
    # -------------------------------------------------------------

    M_out = mach_out_from_mach_in(
        M_in,
        Cf
    )

    # -------------------------------------------------------------
    # Fanno star pressure
    # -------------------------------------------------------------

    p_star = (
        p_in
        /
        p_over_pstar(M_in)
    )

    # -------------------------------------------------------------
    # Outlet static pressure
    # -------------------------------------------------------------

    p_out_calc = (
        p_star
        *
        p_over_pstar(M_out)
    )

    return p_out_calc, M_out


# =====================================================================
# SOLVE INLET MACH NUMBER
# =====================================================================

def solve_inlet_mach(Cf):
    """
    Determine M_in from the prescribed outlet static pressure.

    The nonlinear equation is

        p_out_calculated(M_in) - p_out = 0

    The subsonic solution is selected.
    """

    # -------------------------------------------------------------
    # No-friction case
    # -------------------------------------------------------------

    if Cf == 0.0:

        M_in = brentq(
            lambda M:
                (
                    pt_in / (
                        1.0
                        +
                        0.5 * (gamma - 1.0) * M**2
                    ) ** (
                        gamma / (gamma - 1.0)
                    )
                    -
                    p_out
                ),
            1e-10,
            1.0 - 1e-10,
            xtol=XTOL,
            rtol=RTOL,
        )

        return M_in


    # -------------------------------------------------------------
    # Fanno length parameter
    # -------------------------------------------------------------

    friction_parameter = (
        4.0 * Cf * L / D
    )


    # -------------------------------------------------------------
    # Maximum inlet Mach number before choking
    # -------------------------------------------------------------

    M_in_max = brentq(
        lambda M:
            fanno_function(M)
            -
            friction_parameter,
        1e-8,
        1.0 - 1e-10,
        xtol=XTOL,
        rtol=RTOL,
    )


    # -------------------------------------------------------------
    # Root function
    # -------------------------------------------------------------

    def residual(M_in):

        p_calc, _ = outlet_pressure(
            M_in,
            Cf
        )

        return p_calc - p_out


    # -------------------------------------------------------------
    # Bracket the solution
    # -------------------------------------------------------------

    M_lo = 1e-8

    M_hi = (
        M_in_max
        *
        (1.0 - 1e-10)
    )

    f_lo = residual(M_lo)
    f_hi = residual(M_hi)


    if f_lo * f_hi > 0.0:

        raise RuntimeError(
            f"Could not bracket the inlet Mach number "
            f"for Cf = {Cf:.6g}.\n"
            f"Residual at lower bound = "
            f"{f_lo:.6e} Pa\n"
            f"Residual at upper bound = "
            f"{f_hi:.6e} Pa"
        )


    # -------------------------------------------------------------
    # Robust scalar solution
    # -------------------------------------------------------------

    M_in = brentq(
        residual,
        M_lo,
        M_hi,
        xtol=XTOL,
        rtol=RTOL,
    )

    return M_in


# =====================================================================
# COMPLETE STATE FROM MACH NUMBER
# =====================================================================

def state_from_M(
    M,
    p_star,
    Cf
):
    """
    Calculate the complete local state from Mach number and
    the Fanno star pressure.

    The static pressure is calculated from:

        p = p* (p/p*)

    NOT from the inlet total pressure.

    Total pressure is then calculated from the local static state.
    """

    # -------------------------------------------------------------
    # Static temperature
    # -------------------------------------------------------------

    T = static_temperature(M)


    # -------------------------------------------------------------
    # Static pressure
    # -------------------------------------------------------------

    p = (
        p_star
        *
        p_over_pstar(M)
    )


    # -------------------------------------------------------------
    # Thermodynamic quantities
    # -------------------------------------------------------------

    a = speed_of_sound(T)

    u = M * a

    rho = p / (
        R * T
    )

    mdot = rho * u * A


    # -------------------------------------------------------------
    # Total pressure
    # -------------------------------------------------------------

    pt = total_pressure_from_static(
        p,
        M
    )


    # -------------------------------------------------------------
    # Total temperature
    # -------------------------------------------------------------

    Tt = T * (
        1.0
        +
        0.5 * (gamma - 1.0) * M**2
    )


    # -------------------------------------------------------------
    # Total enthalpy
    # -------------------------------------------------------------

    h = cp * T

    ht = (
        h
        +
        0.5 * u**2
    )


    # -------------------------------------------------------------
    # Entropy
    # -------------------------------------------------------------

    s = entropy_from_state(
        p,
        T
    )


    # -------------------------------------------------------------
    # Wall shear stress
    # -------------------------------------------------------------

    tau_w = (
        0.5
        * rho
        * u**2
        * Cf
    )


    return {
        "M": M,
        "p": p,
        "T": T,
        "rho": rho,
        "a": a,
        "u": u,
        "mdot": mdot,
        "pt": pt,
        "Tt": Tt,
        "h": h,
        "ht": ht,
        "s": s,
        "tau_w": tau_w,
    }


# =====================================================================
# SOLVE COMPLETE CASE
# =====================================================================

def solve_case(Cf):

    # -------------------------------------------------------------
    # Determine inlet Mach number
    # -------------------------------------------------------------

    M_in = solve_inlet_mach(Cf)


    # -------------------------------------------------------------
    # Determine outlet Mach number
    # -------------------------------------------------------------

    M_out = mach_out_from_mach_in(
        M_in,
        Cf
    )


    # -------------------------------------------------------------
    # Determine inlet static pressure
    # -------------------------------------------------------------

    p_in = pt_in / (
        1.0
        +
        0.5 * (gamma - 1.0) * M_in**2
    ) ** (
        gamma / (gamma - 1.0)
    )


    # -------------------------------------------------------------
    # Determine Fanno star pressure
    # -------------------------------------------------------------

    p_star = (
        p_in
        /
        p_over_pstar(M_in)
    )


    # -------------------------------------------------------------
    # Analytical profile
    # -------------------------------------------------------------

    x = np.linspace(
        0.0,
        L,
        N
    )


    M = np.empty_like(x)


    if Cf == 0.0:

        M[:] = M_in

    else:

        F_in = fanno_function(
            M_in
        )

        F_target = (
            F_in
            -
            4.0 * Cf * x / D
        )


        for i, F_target_i in enumerate(
            F_target
        ):

            if F_target_i <= 1e-14:

                M[i] = 1.0

            else:

                M[i] = brentq(
                    lambda Mi:
                        fanno_function(Mi)
                        -
                        F_target_i,
                    M_in,
                    1.0 - 1e-12,
                    xtol=XTOL,
                    rtol=RTOL,
                )


    # -------------------------------------------------------------
    # Calculate complete profile
    # -------------------------------------------------------------

    T = np.empty_like(M)
    p = np.empty_like(M)
    a = np.empty_like(M)
    u = np.empty_like(M)
    rho = np.empty_like(M)
    mdot = np.empty_like(M)
    pt = np.empty_like(M)
    Tt = np.empty_like(M)
    h = np.empty_like(M)
    ht = np.empty_like(M)
    s = np.empty_like(M)
    tau_w = np.empty_like(M)


    for i, Mi in enumerate(M):

        state = state_from_M(
            Mi,
            p_star,
            Cf
        )

        T[i] = state["T"]
        p[i] = state["p"]
        a[i] = state["a"]
        u[i] = state["u"]
        rho[i] = state["rho"]
        mdot[i] = state["mdot"]
        pt[i] = state["pt"]
        Tt[i] = state["Tt"]
        h[i] = state["h"]
        ht[i] = state["ht"]
        s[i] = state["s"]
        tau_w[i] = state["tau_w"]


    # -------------------------------------------------------------
    # Check outlet pressure
    # -------------------------------------------------------------

    p_out_calculated = p[-1]

    pressure_error = (
        p_out_calculated
        -
        p_out
    )


    # -------------------------------------------------------------
    # Package results
    # -------------------------------------------------------------

    return {

        "Cf": Cf,

        "x": x,

        "M": M,
        "p": p,
        "T": T,
        "rho": rho,
        "a": a,
        "u": u,
        "mdot": mdot,

        "pt": pt,
        "Tt": Tt,

        "h": h,
        "ht": ht,
        "s": s,

        "tau_w": tau_w,

        "p_star": p_star,

        "M_in": M_in,
        "M_out": M_out,

        "p_in": p[0],
        "p_out": p[-1],

        "T_in": T[0],
        "T_out": T[-1],

        "rho_in": rho[0],
        "rho_out": rho[-1],

        "u_in": u[0],
        "u_out": u[-1],

        "mdot_in": mdot[0],
        "mdot_out": mdot[-1],

        "pt_in_calculated": pt[0],
        "pt_out_calculated": pt[-1],

        "Tt_in_calculated": Tt[0],
        "Tt_out_calculated": Tt[-1],

        "ht_in_calculated": ht[0],
        "ht_out_calculated": ht[-1],

        "s_in": s[0],
        "s_out": s[-1],

        "pressure_error": pressure_error,
    }


# =====================================================================
# MAIN
# =====================================================================

def main():

    print()
    print("=" * 110)
    print("1-D FANNO FLOW")
    print("=" * 110)


    # -----------------------------------------------------------------
    # Geometry
    # -----------------------------------------------------------------

    print()
    print("Geometry")
    print("-" * 110)

    print(
        f"L               = {L:.6f} m"
    )

    print(
        f"A               = {A:.6f} m^2"
    )

    print(
        f"D               = {D:.6f} m"
    )

    print(
        f"P               = {P:.6f} m"
    )

    print(
        f"P/A             = {P_over_A:.6f} 1/m"
    )


    # -----------------------------------------------------------------
    # Boundary conditions
    # -----------------------------------------------------------------

    print()
    print("Inlet / outlet conditions")
    print("-" * 110)

    print(
        f"pt_in           = {pt_in / 1e3:.6f} kPa"
    )

    print(
        f"Tt_in           = {Tt_in:.6f} K"
    )

    print(
        f"p_out           = {p_out / 1e3:.6f} kPa"
    )


    # -----------------------------------------------------------------
    # Gas
    # -----------------------------------------------------------------

    print()
    print("Gas")
    print("-" * 110)

    print(
        f"gamma           = {gamma:.6f}"
    )

    print(
        f"R               = {R:.6f} J/(kg K)"
    )

    print(
        f"cp              = {cp:.6f} J/(kg K)"
    )


    print()
    print("=" * 110)


    # =================================================================
    # SOLVE ALL CASES
    # =================================================================

    results = []


    for Cf in Cfs:

        try:

            result = solve_case(Cf)

            results.append(result)

        except Exception as error:

            print()
            print(
                f"ERROR for Cf = {Cf:.6g}:"
            )

            print(error)


    # =================================================================
    # SUMMARY
    # =================================================================

    print()
    print("=" * 110)
    print("SUMMARY")
    print("=" * 110)

    print()

    print(
        f"{'Cf':>8s}"
        f"{'M_in':>12s}"
        f"{'M_out':>12s}"
        f"{'p_in [kPa]':>15s}"
        f"{'p_out [kPa]':>15s}"
        f"{'pt_out [kPa]':>16s}"
        f"{'T_in [K]':>14s}"
        f"{'T_out [K]':>15s}"
        f"{'mdot [kg/s]':>16s}"
    )

    print("-" * 126)


    for result in results:

        print(
            f"{result['Cf']:8.4f}"
            f"{result['M_in']:12.7f}"
            f"{result['M_out']:12.7f}"
            f"{result['p_in']/1e3:15.6f}"
            f"{result['p_out']/1e3:15.6f}"
            f"{result['pt_out_calculated']/1e3:16.6f}"
            f"{result['T_in']:14.6f}"
            f"{result['T_out']:15.6f}"
            f"{result['mdot_in']:16.8f}"
        )


    # =================================================================
    # CONSISTENCY CHECKS
    # =================================================================

    print()
    print("=" * 110)
    print("CONSISTENCY CHECKS")
    print("=" * 110)


    for result in results:

        # -------------------------------------------------------------
        # Pressure error
        # -------------------------------------------------------------

        pressure_error = np.max(
            np.abs(
                result["p"][-1]
                -
                p_out
            )
        )


        # -------------------------------------------------------------
        # Mass-flow variation
        # -------------------------------------------------------------

        mass_flow_variation = (
            np.max(result["mdot"])
            -
            np.min(result["mdot"])
        )


        # -------------------------------------------------------------
        # Total-temperature variation
        # -------------------------------------------------------------

        total_temperature_variation = (
            np.max(result["Tt"])
            -
            np.min(result["Tt"])
        )


        # -------------------------------------------------------------
        # Total-enthalpy variation
        # -------------------------------------------------------------

        total_enthalpy_variation = (
            np.max(result["ht"])
            -
            np.min(result["ht"])
        )


        # -------------------------------------------------------------
        # Total-pressure change
        # -------------------------------------------------------------

        total_pressure_change = (
            result["pt"][-1]
            -
            result["pt"][0]
        )


        # -------------------------------------------------------------
        # Entropy change
        # -------------------------------------------------------------

        entropy_change = (
            result["s"][-1]
            -
            result["s"][0]
        )


        print()
        print(
            f"Cf = {result['Cf']:.4f}"
        )

        print(
            f"  outlet pressure error       = "
            f"{pressure_error:.6e} Pa"
        )

        print(
            f"  mass-flow variation         = "
            f"{mass_flow_variation:.6e} kg/s"
        )

        print(
            f"  total-temperature variation = "
            f"{total_temperature_variation:.6e} K"
        )

        print(
            f"  total-enthalpy variation    = "
            f"{total_enthalpy_variation:.6e} J/kg"
        )

        print(
            f"  total-pressure change       = "
            f"{total_pressure_change:.6e} Pa"
        )

        print(
            f"  entropy change              = "
            f"{entropy_change:.6e} J/(kg K)"
        )


        # -------------------------------------------------------------
        # Physical checks
        # -------------------------------------------------------------

        if Cf > 0.0:

            if result["pt"][-1] < result["pt"][0]:

                print(
                    "  OK: total pressure decreases."
                )

            else:

                print(
                    "  WARNING: total pressure did not decrease."
                )


            if result["s"][-1] > result["s"][0]:

                print(
                    "  OK: entropy increases."
                )

            else:

                print(
                    "  WARNING: entropy did not increase."
                )


        if total_temperature_variation < 1e-8:

            print(
                "  OK: total temperature is constant."
            )

        else:

            print(
                "  WARNING: total temperature is not constant."
            )


        if mass_flow_variation < 1e-10:

            print(
                "  OK: mass flow is constant."
            )

        else:

            print(
                "  WARNING: mass flow is not constant."
            )


    # =================================================================
    # PLOTS
    # =================================================================

    # -----------------------------------------------------------------
    # Mach number
    # -----------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    for result in results:

        ax.plot(
            result["x"],
            result["M"],
            label=fr"$C_f={result['Cf']}$",
        )

    ax.set_xlabel(
        r"$x$ [m]"
    )

    ax.set_ylabel(
        r"$M$"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()


    # -----------------------------------------------------------------
    # Static pressure
    # -----------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    for result in results:

        ax.plot(
            result["x"],
            result["p"] / 1e3,
            label=fr"$C_f={result['Cf']}$",
        )

    ax.set_xlabel(
        r"$x$ [m]"
    )

    ax.set_ylabel(
        r"$p$ [kPa]"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()


    # -----------------------------------------------------------------
    # Total pressure
    # -----------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    for result in results:

        ax.plot(
            result["x"],
            result["pt"] / 1e3,
            label=fr"$C_f={result['Cf']}$",
        )

    ax.set_xlabel(
        r"$x$ [m]"
    )

    ax.set_ylabel(
        r"$p_t$ [kPa]"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()


    # -----------------------------------------------------------------
    # Static temperature
    # -----------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    for result in results:

        ax.plot(
            result["x"],
            result["T"],
            label=fr"$C_f={result['Cf']}$",
        )

    ax.set_xlabel(
        r"$x$ [m]"
    )

    ax.set_ylabel(
        r"$T$ [K]"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()


    # -----------------------------------------------------------------
    # Total temperature
    # -----------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    for result in results:

        ax.plot(
            result["x"],
            result["Tt"],
            label=fr"$C_f={result['Cf']}$",
        )

    ax.set_xlabel(
        r"$x$ [m]"
    )

    ax.set_ylabel(
        r"$T_t$ [K]"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()


    # -----------------------------------------------------------------
    # Entropy
    # -----------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    for result in results:

        ax.plot(
            result["x"],
            result["s"],
            label=fr"$C_f={result['Cf']}$",
        )

    ax.set_xlabel(
        r"$x$ [m]"
    )

    ax.set_ylabel(
        r"$s$ [J/(kg K)]"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()


    # =================================================================
    # SAVE RESULTS
    # =================================================================

    with open(
        "results.pkl",
        "wb"
    ) as file:

        pickle.dump(
            results,
            file
        )


    print()
    print("=" * 110)
    print("Results saved to results.pkl")
    print("=" * 110)


# =====================================================================
# RUN
# =====================================================================

if __name__ == "__main__":

    main()

    plt.show()