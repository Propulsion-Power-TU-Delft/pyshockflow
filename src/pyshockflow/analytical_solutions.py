import numpy as np
from CoolProp.CoolProp import PropsSI, PhaseSI
from scipy.integrate import solve_ivp


def solve_nonideal_expansion(densityStart, temperatureStart, densityFinal, numberPoints, fluid):
    """Solve Equation 11 of Steady, 'isentropic flows of dense gases' by Cramer. The result
    is the Mach number distribution as a function of density

    Args:
        densityStart (float): initial density 
        temperatureStart (float): initial value
        densityFinal (float): final value at the end of expansion process
        numberPoints (int): discretization
        fluid (float): coolprop name

    Returns:
        solution arrays of rho and M
    """
    s = PropsSI('S', 'T', temperatureStart, 'D', densityStart, fluid)
    rho = np.linspace(densityStart, densityFinal, numberPoints)
    M = np.zeros(numberPoints)
    
    s = np.zeros_like(rho)
    bigGamma = np.zeros_like(rho)
    for i in range(numberPoints):
        s[i] = PropsSI('S', 'T', temperatureStart, 'D', rho[0], fluid)
        bigGamma[i] = PropsSI('FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS', 'S', s[i], 'D', rho[i], fluid)
    
    # integrate equations with trapz method. The coefficients are obtained discretizing the gov. equation
    drho = rho[1] - rho[0]
    for i in range(1, numberPoints):
        alpha = 1.0 - 0.5*(bigGamma[i] + bigGamma[i-1])
        beta = drho / (rho[i]+rho[i-1])
        A = 1.0 - alpha*beta
        B = -2.0*alpha*beta*M[i-1]
        C = -M[i-1]**2 * (1+alpha*beta) + 4.0*beta
        M[i] = (-B + np.sqrt(B**2 - 4.0*A*C)) / (2.0*A)
    return rho, M



def find_area_ratio_distribution(density, mach):
    """Integrate area distribution equation. Implicit scheme

    Args:
        density (array)
        mach (arrays)

    Returns:
        area: area distribution
    """
    area = np.zeros_like(density)+1
    nPoints = len(density)
    for i in range(5, nPoints):
        drho = density[i] - density[i-1]
        alpha = (1-mach[i]**2)/mach[i]**2 * drho/density[i]
        area[i] = area[i-1] / (1-alpha)
        
    return area
    