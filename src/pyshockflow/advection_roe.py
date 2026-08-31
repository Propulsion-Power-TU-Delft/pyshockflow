import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from numpy import sqrt
from pyshockflow import FluidIdeal
from pyshockflow.math_utils import *


class AdvectionRoeBase:
    def __init__(self, rhoL, rhoR, uL, uR, pL, pR, fluid):
        """
        Roe scheme numerics for ideal gas. Parameters are left and right values of density, velocity and pressure, and the fluid object.
        Formulation based on x-split Riemann Solver in the book by Toro.
        """
        self.rhoL = rhoL
        self.rhoR = rhoR
        self.uL = uL
        self.uR = uR
        self.pL = pL
        self.pR = pR
        self.fluid = fluid
        if isinstance(fluid, FluidIdeal):
            self.gmma = fluid.gmma
        self.eL = fluid.computeStaticEnergy_p_rho(pL, rhoL)
        self.eR = fluid.computeStaticEnergy_p_rho(pR, rhoR)
        self.htL = self.computeTotalEnthalpy(rhoL, uL, pL, self.eL)
        self.htR = self.computeTotalEnthalpy(rhoR, uR, pR, self.eR)
        self.u1L, self.u2L, self.u3L = getConservativesFromPrimitives(rhoL, uL, pL, self.fluid)
        self.u1R, self.u2R, self.u3R = getConservativesFromPrimitives(rhoR, uR, pR, self.fluid)
        self.aL = self.fluid.computeSoundSpeed_p_rho(self.pL, self.rhoL)
        self.aR = self.fluid.computeSoundSpeed_p_rho(self.pR, self.rhoR)


    def computeRoeAvg(self, fL, fR):
        """
        Roe Averaging Operator
        """
        favg = (sqrt(self.rhoL)*fL + sqrt(self.rhoR)*fR)/(sqrt(self.rhoL)+ sqrt(self.rhoR))
        return favg

    
    def computeAveragedVariables(self):
        """
        compute the Roe averaged variables for the 1D Euler equations
        """
        self.rhoAVG = sqrt(self.rhoL*self.rhoR)
        self.uAVG = self.computeRoeAvg(self.uL, self.uR)
        self.hAVG = self.computeRoeAvg(self.htL, self.htR)
        self.aAVG = sqrt((self.gmma-1)*(self.hAVG-0.5*self.uAVG**2))
    
    
    def computeTotalEnthalpy(self, rho, u, p, e):
        et = 0.5*u**2 + e
        ht = et+p/rho
        return ht
    

    def computeAveragedEigenvalues(self):
        """
        compute eigenvalues of the averaged Jacobian
        """
        self.lambda_vec = np.array([self.uAVG-self.aAVG, 
                                    self.uAVG, 
                                    self.uAVG+self.aAVG])
    

    def computeAveragedEigenvectors(self):
        """
        compute eigenvector matrix of the averaged flux Jacobian
        """
        self.eigenvector_mat = np.zeros((3, 3))
        
        self.eigenvector_mat[0, 0] = 1
        self.eigenvector_mat[1, 0] = self.uAVG-self.aAVG
        self.eigenvector_mat[2, 0] = self.hAVG-self.uAVG*self.aAVG

        self.eigenvector_mat[0, 1] = 1
        self.eigenvector_mat[1, 1] = self.uAVG
        self.eigenvector_mat[2, 1] = 0.5*self.uAVG**2

        self.eigenvector_mat[0, 2] = 1
        self.eigenvector_mat[1, 2] = self.uAVG+self.aAVG
        self.eigenvector_mat[2, 2] = self.hAVG+self.uAVG*self.aAVG
    

    def computeWaveStrengths(self):
        """
        Characteristic jumps due to initial conditions
        """
        self.alphas = np.zeros(3)
        self.alphas[0] = 1/2/self.aAVG**2 *(self.pR-self.pL-self.rhoAVG*self.aAVG*(self.uR-self.uL))
        self.alphas[1] = self.rhoR-self.rhoL - (self.pR-self.pL)/self.aAVG**2
        self.alphas[2] = 1/2/self.aAVG**2*(self.pR-self.pL + self.rhoAVG*self.aAVG*(self.uR-self.uL))
        

    def computeFlux(self, entropyFixActive, fixCoefficient):
        """
        compute the Roe flux. The flux is computed for 1D problems.
        """
        self.computeAveragedVariables()
        self.computeAveragedEigenvalues()
        self.computeAveragedEigenvectors()
        self.computeWaveStrengths()
        
        fluxL = self.EulerFlux(self.u1L, self.u2L, self.u3L)
        fluxR = self.EulerFlux(self.u1R, self.u2R, self.u3R)
        fluxRoe = 0.5*(fluxL+fluxR)

        # compute the entropy fixed abs eigenvalues
        if entropyFixActive==False:
            absEig = np.abs(self.lambda_vec)
        else:
            absEig = applyEntropyFix(self.lambda_vec, self.aAVG, fixCoefficient)

        for iDim in range(3):
            for jVec in range(3):
                fluxRoe[iDim] -= 0.5*self.alphas[jVec]*absEig[jVec]*self.eigenvector_mat[iDim, jVec]
        
        return fluxRoe
        
    def EulerFlux(self, u1, u2, u3):
        """
        Get the Euler flux starting from conservative variables. 
        """
        flux1D = computeAdvectionFluxFromConservatives(u1, u2, u3, self.fluid)
        return flux1D



class AdvectionRoeArabi(AdvectionRoeBase):
    """
    Generalised Roe Scheme for real gases, taken from the article 'A simple extension of Roe scheme for real gases', Arabi et al. 
    Journal of Computational Physics 2017. Formulation based on 1D problem.
    """
    def __init__(self, rhoL, rhoR, uL, uR, pL, pR, fluid):
        super().__init__(rhoL, rhoR, uL, uR, pL, pR, fluid)
        self.deltaP = (self.pR-self.pL)
        self.deltaU = (self.uR - self.uL)
        self.deltaRho = (self.rhoR - self.rhoL)
    
    
    def computeAveragedVariables(self):
        """
        compute the Roe averaged variables for the 1D Euler equations
        """
        self.rhoAVG = sqrt(self.rhoL*self.rhoR)
        self.uAVG = self.computeRoeAvg(self.uL, self.uR)
        self.hAVG = self.computeRoeAvg(self.htL, self.htR)
        self.aAVG = self.computeRoeAvg(self.aL, self.aR)


    def computeWaveStrengths(self):
        self.alphas = np.zeros(3)
        self.alphas[0] = 1/2/self.aAVG**2*(self.deltaP+self.rhoAVG*self.aAVG*self.deltaU)
        self.alphas[1] = 1/2/self.aAVG**2*(self.deltaP-self.rhoAVG*self.aAVG*self.deltaU)
        self.alphas[2] = self.deltaRho-self.deltaP/self.aAVG**2
    

    def computeAveragedEigenvalues(self):
        self.lambda_vec = np.array([self.uAVG+self.aAVG, 
                                    self.uAVG-self.aAVG,
                                    self.uAVG])
    

    def computeFlux(self, entropyFixActive, fixCoefficient):
        """
        Assemble the global flux, average + dissipation, following the approach of the article
        """
        self.computeAveragedVariables()
        self.computeAveragedEigenvalues()
        self.computeAveragedEigenvectors()
        self.computeWaveStrengths()
        
        fluxL = self.EulerFlux(self.u1L, self.u2L, self.u3L)
        fluxR = self.EulerFlux(self.u1R, self.u2R, self.u3R)

        # compute the entropy fixed abs eigenvalues
        if entropyFixActive==False:
            absEig = np.abs(self.lambda_vec)
        else:
            absEig = applyEntropyFix(self.lambda_vec, self.aAVG, fixCoefficient)

        deltaF = np.zeros(3)
        deltaF[0] = absEig[0]*self.alphas[0] + absEig[1]*self.alphas[1] + absEig[2]*self.alphas[2]
        deltaF[1] = (self.uAVG+self.aAVG)*absEig[0]*self.alphas[0] + (self.uAVG-self.aAVG)*absEig[1]*self.alphas[1] + self.uAVG*absEig[2]*self.alphas[2]

        X = (self.rhoR*self.uR*self.htR)-(self.rhoL*self.uL*self.htL)- \
            (self.hAVG+self.uAVG*self.aAVG)*(self.uAVG+self.aAVG)*(1/2/self.aAVG**2*(self.deltaP+self.rhoAVG*self.aAVG*self.deltaU)) - \
            (self.hAVG-self.uAVG*self.aAVG)*(self.uAVG-self.aAVG)*(1/2/self.aAVG**2*(self.deltaP-self.rhoAVG*self.aAVG*self.deltaU))
        
        if (self.uAVG>=0):
            pass
        else:
            X *= -1
        
        deltaF[2] = (self.hAVG+self.uAVG*self.aAVG)*absEig[0]*(self.alphas[0]) + \
                         (self.hAVG-self.uAVG*self.aAVG)*absEig[1]*(self.alphas[1]) + X
                         
        fluxRoe = 0.5*(fluxL+fluxR) - 0.5*deltaF
        return fluxRoe
    
    

class AdvectionRoeVinokur(AdvectionRoeBase):
    """
    Generalised Roe Scheme for real gases, 
    where the Roe avg state is taken from the article 'Generalized flux-vector splitting and Roe average for an equilibrium real gas', Vinokur and Montagnè 
    Journal of Computational Physics 1990. Formulation based on 1D problem.
    """
    def __init__(self, rhoL, rhoR, uL, uR, pL, pR, fluid):
        super().__init__(rhoL, rhoR, uL, uR, pL, pR, fluid)
        self.deltaP = (self.pR-self.pL)
        self.deltaU = (self.uR - self.uL)
        self.deltaRho = (self.rhoR - self.rhoL)
    
    
    def computeAveragedVariables(self):
        """
        compute the Roe averaged state following the approach described in the articleof Vinokur
        """
        alpha = np.sqrt(self.rhoL) / (np.sqrt(self.rhoL)+np.sqrt(self.rhoR))
        self.uAVG = alpha*self.uL + (1-alpha)*self.uR
        self.htAVG = alpha*self.htL + (1-alpha)*self.htR
        self.hL = self.htL - 0.5*self.uL**2
        self.hR = self.htR - 0.5*self.uR**2
        self.hAVG = alpha*self.hL + (1-alpha)*self.hR + 0.5*alpha*(1.0-alpha)*self.deltaU**2
        
        # compute mean initial guess state
        p_mean = 0.5*(self.pL+self.pR)
        eL = self.fluid.computeStaticEnergy_p_rho(self.pL, self.rhoL)
        rho_mean = 0.5*(self.rhoL+self.rhoR)
        rhoeL = self.rhoL*eL
        eR = self.fluid.computeStaticEnergy_p_rho(self.pR, self.rhoR)
        rhoeR = self.rhoR*eR
        rhoe_mean = 0.5*(rhoeL+rhoeR)
        e_mean = rhoe_mean/rho_mean
        
        chiL, kappaL = self.fluid.computeChiKappa_VinokurScheme_p_rho(self.pL, self.rhoL)
        chiR, kappaR = self.fluid.computeChiKappa_VinokurScheme_p_rho(self.pR, self.rhoR)
        chiM, kappaM = self.fluid.computeChiKappa_VinokurScheme_p_rho(p_mean, rho_mean)
        chiHat = (chiL + chiR + 4.0*chiM) / 6.0
        kappaHat = (kappaL + kappaR + 4.0*kappaM) / 6.0
        delta_rhoe = (rhoeR - rhoeL)
        
        # projection procedure to compute the average state starting fro the initial guess (hat values)
        error_term = self.deltaP - chiHat*self.deltaRho - kappaHat*delta_rhoe
        hM = 0.5*(self.hL+self.hR)
        kappah_hat = (kappaL*self.hL + kappaR*self.hR + 4.0*kappaM*hM) / 6.0
        csquare_L = chiL + kappaL*self.hL
        csquare_R = chiR + kappaR*self.hR
        csquare_M = chiM + kappaM*hM
        sHat = (csquare_L + csquare_R + 4.0*csquare_M) / 6.0
        D_term = (sHat*self.deltaRho)**2 + (self.deltaP)**2
        if self.deltaRho==0:
            self.chiAVG = chiHat
        else:
            self.chiAVG = (D_term * chiHat + sHat**2 * self.deltaRho * error_term) / (D_term - self.deltaP*error_term)
        
        if self.deltaP==0:
            self.kappaAVG = kappaHat
        else:
            self.kappaAVG = (D_term * kappaHat) / (D_term - self.deltaP*error_term)
        
        self.aAVG = np.sqrt(self.chiAVG + self.kappaAVG*self.hAVG)
    
    
    def computeFlux(self, entropyFixActive, fixCoefficient):
        """
        compute the global flux, average + dissipation
        """
        fluxL = self.EulerFlux(self.u1L, self.u2L, self.u3L)
        fluxR = self.EulerFlux(self.u1R, self.u2R, self.u3R)

        # compute the Eigenvectors matrices
        k1 = 0.5*self.kappaAVG*self.uAVG**2 + self.kappaAVG
        k2 = 0.5*self.uAVG**2 - self.chiAVG/self.kappaAVG
        
        # right eigenvectors matrix
        matrixR = np.array([[1, 1, 1],
                            [self.uAVG, self.uAVG+self.aAVG, self.uAVG-self.aAVG],
                            [k2, self.htAVG + self.aAVG*self.uAVG, self.htAVG - self.aAVG*self.uAVG]])
        
        # left eigenvectors matrix
        matrixRinv = np.array([[1-k1/self.aAVG**2, self.kappaAVG*self.uAVG/self.aAVG**2, -self.kappaAVG/self.aAVG**2],
                               [0.5*(k1/self.aAVG**2-self.uAVG/self.aAVG), -0.5*(self.kappaAVG*self.uAVG/self.aAVG**2-1/self.aAVG), 0.5*self.kappaAVG/self.aAVG**2],
                               [0.5*(k1/self.aAVG**2+self.uAVG/self.aAVG), -0.5*(self.kappaAVG*self.uAVG/self.aAVG**2+1/self.aAVG), 0.5*self.kappaAVG/self.aAVG**2]])
                
        # eigenvalues, to fix
        eigsAVG = np.array([self.uAVG, self.uAVG+self.aAVG, self.uAVG-self.aAVG])
        if entropyFixActive==False:
            absEig = np.abs(eigsAVG)
        else:
            absEig = applyEntropyFix(eigsAVG, self.aAVG, fixCoefficient)
        
        # eigenvalues matrix
        matrixLambda = np.diag(absEig)
        
        # compute the Flux
        deltaU = np.array([self.u1R-self.u1L, self.u2R-self.u2L, self.u3R-self.u3L]).reshape(3,1)
        deltaFlux = matrixR @ matrixLambda @ matrixRinv @ deltaU
        fluxRoe = 0.5*(fluxL+fluxR) - 0.5*deltaFlux.flatten()
        return fluxRoe
        

def applyEntropyFix(eigs, aAVG, kappa):
    """
    Apply Harten entropy fix to eigenvalues.
    
    eigs : ndarray of shape (3,)
        Raw Roe eigenvalues [u, u+a, u-a].
    aAVG : float
        Roe-averaged sound speed.
    kappa : float
        Fix coefficient (default 0.2).
    """
    delta = kappa * aAVG
    abs_eigs = np.abs(eigs)
    return np.where(abs_eigs < delta,
                    0.5 * (eigs**2 / delta + delta),
                    abs_eigs)
















