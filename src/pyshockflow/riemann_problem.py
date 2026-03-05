import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os
import pickle
from pyshockflow.plot_styles import *


class RiemannProblem:

    def __init__(self, x, t):
        """
        Initialize the Riemann Problem
        :param x: array with x-coordinates (centered around zero)
        :param t: arrays of time-instants
        """
        self.equations = 'Euler'
        self.x = x
        self.nx = len(self.x)
        self.t = t
        self.nt = len(self.t)
        self.gmma = 1.4

    def initializeState(self, inState):
        """
        Initialize the left and right states
        :param inState: array with values (density left, density right, velocity left, velocity right, pressure left,
        pressure right)
        """
        self.rhoL = inState[0]
        self.rhoR = inState[1]
        self.uL = inState[2]
        self.uR = inState[3]
        self.pL = inState[4]
        self.pR = inState[5]
        self.eL = self.computeEnergy_pRho(self.pL, self.rhoL)
        self.eR = self.computeEnergy_pRho(self.pR, self.rhoR)
        self.aL = np.sqrt(self.gmma * self.pL / self.rhoL)
        self.aR = np.sqrt(self.gmma * self.pR / self.rhoR)

    def computeEnergy_pRho(self, p, rho):
        """
        compute energy for an ideal gas given pressure and density
        :param p: pressure
        :param rho: density
        :return: energy
        """
        return p / (self.gmma - 1) / rho

    def initializeSolutionArrays(self):
        """
        Initialize the containers for the solution
        :return:
        """
        nx = len(self.x)
        nt = len(self.t)

        self.rho = np.zeros((nx, nt))
        self.u = np.zeros((nx, nt))
        self.p = np.zeros((nx, nt))
        self.e = np.zeros((nx, nt))

        self.rho[:, 0] = self.copyInitialState(self.rhoL, self.rhoR)
        self.u[:, 0] = self.copyInitialState(self.uL, self.uR)
        self.p[:, 0] = self.copyInitialState(self.pL, self.pR)
        self.e[:, 0] = self.copyInitialState(self.eL, self.eR)

    def copyInitialState(self, fL, fR):
        """
        Given left and right values, copy these values along the x-axis
        :param fL:
        :param fR:
        :return:
        """
        xmean = 0.5 * (self.x[-1] + self.x[0])
        f = np.zeros_like(self.x)
        for i in range(len(self.x)):
            if self.x[i] <= xmean:
                f[i] = fL
            else:
                f[i] = fR
        return f

    def plotSolution(self, iTime, folder_name = None, file_name = None):
        """
        Plot the solution at time instant element iTime
        :param iTime: element index in time array
        :param folder_name:
        :param file_name:
        :return:
        """
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        ax[0, 0].plot(self.x, self.rho[:, iTime], '-C0o', ms=2)
        ax[0, 0].set_ylabel(r'Density')

        ax[0, 1].plot(self.x, self.u[:, iTime], '-C1o', ms=2)
        ax[0, 1].set_ylabel(r'Velocity')

        ax[1, 0].plot(self.x, self.p[:, iTime], '-C2o', ms=2)
        ax[1, 0].set_ylabel(r'Pressure')

        ax[1, 1].plot(self.x, self.e[:, iTime], '-C3o', ms=2)
        ax[1, 1].set_ylabel(r'Energy')

        fig.suptitle('Time %.3f' % self.t[iTime])

        for row in ax:
            for col in row:
                col.set_xlabel('x')
                col.grid(alpha=.3)

        if file_name is not None and folder_name is not None:
            os.makedirs(folder_name, exist_ok=True)
            plt.savefig(folder_name + '/' + file_name + '.pdf', bbox_inches='tight')


    def computeStarRegion(self):
        """
        Given left and right states, compute the star pressure and velocity values.
        :return:
        """
        self.pStar = self.computeStarPressure()[0]
        fL, fR = self.compute_fL_fR(self.pStar)
        self.uStar = 0.5 * (self.uL + self.uR) + 0.5 * (fR - fL)

        # right wave
        if self.pStar > self.pR:
            # print('The right wave is a shock wave')
            self.right_wave = 'shock'
        else:
            # print('The right wave is a rarefaction wave')
            self.right_wave = 'rarefaction'

        # left wave
        if self.pStar > self.pL:
            # print('The left wave is a shock wave')
            self.left_wave = 'shock'
        else:
            # print('The left wave is a rarefaction wave')
            self.left_wave = 'rarefaction'

    def solve(self, space_domain='global', time_domain='global'):
        """
        Sample the solution in the x,t domains. Follows methodology given in the book Riemann solvers by Toro.
        :return:
        """
        if space_domain=='global':
            self.ix_values = range(self.nx)
        elif space_domain=='interface':
            self.ix_values = [self.nx//2]
        else:
            raise ValueError('Unknown value of space_domain')
        
        if time_domain=='global':
            self.it_values = range(1, self.nt)
        elif time_domain=='interface':
            self.it_values = [self.nt//2]
        else:
            raise ValueError('Unknown value of time_domain')


        for it in self.it_values:
            for ix in self.ix_values:
                S = self.x[ix] / self.t[it]
                if S < self.uStar:  # we are at the left of the contact discontinuity
                    if self.pStar > self.pL:  # the left wave is a shock wave
                        rhoStarL = self.rhoL * ((self.pStar / self.pL) + (self.gmma - 1) / (self.gmma + 1)) / (
                                (self.gmma - 1) / (self.gmma + 1) * self.pStar / self.pL + 1)
                        self.Sl = self.uL - self.aL * ((self.gmma + 1) / 2 / self.gmma * self.pStar / self.pL + (
                                    self.gmma - 1) / 2 / self.gmma) ** 0.5
                        if S < self.Sl:  # we are at the left of the shock, in the left zone
                            self.rho[ix, it] = self.rhoL
                            self.u[ix, it] = self.uL
                            self.p[ix, it] = self.pL
                            self.e[ix, it] = self.eL
                        else:
                            self.rho[ix, it] = rhoStarL
                            self.u[ix, it] = self.uStar
                            self.p[ix, it] = self.pStar
                            self.e[ix, it] = self.pStar / (self.gmma - 1) / rhoStarL
                    else:  # the left wave is a rarefaction wave
                        self.Shl = self.uL - self.aL
                        if S < self.Shl:  # we are at the left of the rarefaction
                            self.rho[ix, it] = self.rhoL
                            self.u[ix, it] = self.uL
                            self.p[ix, it] = self.pL
                            self.e[ix, it] = self.eL
                        else:
                            rhoStarL = self.rhoL * (self.pStar / self.pL) ** (1 / self.gmma)
                            aStarL = self.aL * (self.pStar / self.pL) ** ((self.gmma - 1) / 2 / self.gmma)
                            self.Stl = self.uStar - aStarL
                            if S > self.Stl:  # we are in the star region, between the rarefaction and the contact
                                self.rho[ix, it] = rhoStarL
                                self.u[ix, it] = self.uStar
                                self.p[ix, it] = self.pStar
                                self.e[ix, it] = self.pStar / (self.gmma - 1) / rhoStarL
                            else:  # we are inside the rarefaction wave
                                self.rho[ix, it], self.u[ix, it], self.p[ix, it], self.e[ix, it] = self.computeLeftFanQuantities(S)
                else:  # we are at the right of the discontinuity
                    if self.pStar > self.pR:  # the right wave is a shock wave
                        rhoStarR = self.rhoR * ((self.pStar / self.pR) + (self.gmma - 1) / (self.gmma + 1)) / (
                                (self.gmma - 1) / (self.gmma + 1) * self.pStar / self.pR + 1)
                        self.Sr = self.uR + self.aR * ((self.gmma + 1) / 2 / self.gmma * self.pStar / self.pR + (
                                self.gmma - 1) / 2 / self.gmma) ** 0.5
                        if S > self.Sr:  # we are in the right state zone
                            self.rho[ix, it] = self.rhoR
                            self.u[ix, it] = self.uR
                            self.p[ix, it] = self.pR
                            self.e[ix, it] = self.eR
                        else:  # we are in the right star zone
                            self.rho[ix, it] = rhoStarR
                            self.u[ix, it] = self.uStar
                            self.p[ix, it] = self.pStar
                            self.e[ix, it] = self.pStar / (self.gmma - 1) / rhoStarR
                    else:  # the right wave is a rarefaction wave
                        self.Shr = self.uR + self.aR
                        if S > self.Shr:  # we are at the right of the rarefaction
                            self.rho[ix, it] = self.rhoR
                            self.u[ix, it] = self.uR
                            self.p[ix, it] = self.pR
                            self.e[ix, it] = self.eR
                        else:
                            rhoStarR = self.rhoR * (self.pStar / self.pR) ** (1 / self.gmma)
                            aStarR = self.aR * (self.pStar / self.pR) ** ((self.gmma - 1) / 2 / self.gmma)
                            self.Str = self.uStar + aStarR
                            if S < self.Str:  # we are in the right star region, between the the contact and the rarefaction
                                self.rho[ix, it] = rhoStarR
                                self.u[ix, it] = self.uStar
                                self.p[ix, it] = self.pStar
                                self.e[ix, it] = self.pStar / (self.gmma - 1) / rhoStarR
                            else:  # we are inside the rarefaction wave
                                self.rho[ix, it], self.u[ix, it], self.p[ix, it], self.e[ix, it] = self.computeRightFanQuantities(S)

    def computeLeftFanQuantities(self, S):
        """
        Function values for left rarefaction fan
        :param S: similarity parameter (=x/t)
        :return: rho, u, p, e inside the fan
        """
        rho = self.rhoL * (2 / (self.gmma + 1) + (self.gmma - 1) / (self.gmma + 1) / self.aL * (self.uL - S)) ** (
                2 / (self.gmma - 1))
        u = 2 / (self.gmma + 1) * (self.aL + (self.gmma - 1) / 2 * self.uL + S)
        p = self.pL * (2 / (self.gmma + 1) + (self.gmma - 1) / (self.gmma + 1) / self.aL * (self.uL - S)) ** (
                2 * self.gmma / (self.gmma - 1))
        e = p / (self.gmma - 1) / rho
        return rho, u, p, e

    def computeRightFanQuantities(self, S):
        """
        Function values for right rarefaction fan
        :param S: similarity parameter (=x/t)
        :return: rho, u, p, e inside the fan
        """
        rho = self.rhoR * (2 / (self.gmma + 1) - (self.gmma - 1) / (self.gmma + 1) / self.aR * (self.uR - S)) ** (
                2 / (self.gmma - 1))
        u = 2 / (self.gmma + 1) * (-self.aR + (self.gmma - 1) / 2 * self.uR + S)
        p = self.pR * (2 / (self.gmma + 1) - (self.gmma - 1) / (self.gmma + 1) / self.aR * (self.uR - S)) ** (
                2 * self.gmma / (self.gmma - 1))
        e = p / (self.gmma - 1) / rho
        return rho, u, p, e

    def compute_fL_fR(self, p):
        """
        compute the functions needed to evaluate p*
        :param p: pressure value
        :return: value of the function (to find zero of)
        """
        if p > self.pL:
            AL = 2 / (self.gmma + 1) / self.rhoL
            BL = (self.gmma - 1) / (self.gmma + 1) * self.pL
            fL = (p - self.pL) * (AL / (p + BL)) ** 0.5
        else:
            al = np.sqrt(self.gmma * self.pL / self.rhoL)
            fL = 2 * al / (self.gmma - 1) * ((p / self.pL) ** ((self.gmma - 1) / 2 / self.gmma) - 1)

        if p > self.pR:
            AR = 2 / (self.gmma + 1) / self.rhoR
            BR = (self.gmma - 1) / (self.gmma + 1) * self.pR
            fR = (p - self.pR) * (AR / (p + BR)) ** 0.5
        else:
            ar = np.sqrt(self.gmma * self.pR / self.rhoR)
            fR = 2 * ar / (self.gmma - 1) * ((p / self.pR) ** ((self.gmma - 1) / 2 / self.gmma) - 1)
        return fL, fR

    def computeStarPressure(self):
        """
        Find p* value making use of fsolve
        :return: p*
        """
        def pFunc(p):
            dU = self.uR - self.uL
            fL, fR = self.compute_fL_fR(p)
            return fL + fR + dU

        # guess initial value of pstar. This seems like working for every test-case
        if self.pL == self.pR:
            pGuess = 0
        else:
            pGuess = 0.5 * (self.pR + self.pL)
        pSol = fsolve(pFunc, pGuess)
        return pSol

    def showAnimation(self):
        """
        Show animation of the results for all time instants
        """
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        for it in range(self.nt):
            for row in ax:
                for col in row:
                    col.cla()
            ax[0, 0].plot(self.x, self.rho[:, it], '-C0o', ms=2)
            ax[0, 0].set_ylabel(r'Density')

            ax[0, 1].plot(self.x, self.u[:, it], '-C1o', ms=2)
            ax[0, 1].set_ylabel(r'Velocity')

            ax[1, 0].plot(self.x, self.p[:, it], '-C2o', ms=2)
            ax[1, 0].set_ylabel(r'Pressure')

            ax[1, 1].plot(self.x, self.e[:, it], '-C3o', ms=2)
            ax[1, 1].set_ylabel(r'Energy')

            fig.suptitle('Time %.3f' % self.t[it])

            for row in ax:
                for col in row:
                    col.set_xlabel('x')
                    col.grid(alpha=.3)
            plt.pause(1e-3)
    
    def plotLocalSolEvolution(self):
        """
        Plot the values in time of all the points specified in self.ix_values, defined when solving the Riemann problem
        """
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        for ix in self.ix_values:
    
            ax[0, 0].plot(self.t, self.rho[ix, :], '-C0o', ms=2)
            ax[0, 0].set_ylabel(r'Density')

            ax[0, 1].plot(self.t, self.u[ix, :], '-C1o', ms=2)
            ax[0, 1].set_ylabel(r'Velocity')

            ax[1, 0].plot(self.t, self.p[ix, :], '-C0o', ms=2)
            ax[1, 0].set_ylabel(r'Pressure')

            ax[1, 1].plot(self.t, self.e[ix, :], '-C0o', ms=2)
            ax[1, 1].set_ylabel(r'Energy')

            fig.suptitle('x: %.3f' % self.x[ix])

            for row in ax:
                for col in row:
                    col.set_xlabel('t')
                    col.grid(alpha=.3)
    

    def getSolutionInTime(self):
        """
        Return the arrays storing the solutions in time
        """
        rho = self.rho[self.ix_values, :].flatten()
        u = self.u[self.ix_values, :].flatten()
        p = self.p[self.ix_values, :].flatten()
        return rho, u, p
    
    def getSolutionInSpace(self):
        """
        Return the arrays storing the solutions in space
        """
        rho = self.rho[:, self.it_values].flatten()
        u = self.u[:, self.it_values].flatten()
        p = self.p[:, self.it_values].flatten()
        return rho, u, p


    def drawSpaceTimePlot(self, folder_name = None, file_name = None):
        """
        Just draw the sketch of the space-time plot showing the characters of the waves
        """
        x0 = 0
        t0 = 0
        p = np.linspace(0, 10, 100)

        plt.figure(figsize=(4, 2.75))

        # contact wave
        alpha = np.arctan2(1 , self.uStar)
        plt.plot(x0 + np.cos(alpha) * p, t0 + np.sin(alpha) * p, '--k', label='Contact wave')

        # left wave
        if self.left_wave == 'shock':
            alpha = np.arctan2(1, self.Sl)
            plt.plot(x0+np.cos(alpha)*p, t0+np.sin(alpha)*p, 'C0', label='Shock')
        else:
            alpha1 = np.arctan2(1, self.Shl)
            alpha2 = np.arctan2(1, self.Stl)
            alphas = np.linspace(alpha1, alpha2, 15)
            for alpha in alphas:
                xt, yt = x0 + np.cos(alpha) * p, t0 + np.sin(alpha) * p
                if alpha==alphas[0]:
                    plt.plot(xt, yt, 'C0', label='Rarefaction')
                elif alpha==alphas[-1]:
                    plt.plot(xt, yt, 'C0')
                else:
                    plt.plot(xt, yt, 'C0', lw=0.75)

        # right wave
        if self.right_wave == 'shock':
            alpha = np.arctan2(1, self.Sr)
            plt.plot(x0 + np.cos(alpha) * p, t0 + np.sin(alpha) * p, 'C1', label='Shock')
        else:
            alpha1 = np.arctan2(1, self.Shr)
            alpha2 = np.arctan2(1, self.Str)
            alphas = np.linspace(alpha1, alpha2, 15)
            for alpha in alphas:
                xt, yt = x0 + np.cos(alpha) * p, t0 + np.sin(alpha) * p
                if alpha == alphas[0]:
                    plt.plot(xt, yt, 'C1', label='Rarefaction')
                elif alpha == alphas[-1]:
                    plt.plot(xt, yt, 'C1')
                else:
                    plt.plot(xt, yt, 'C1', lw=0.75)

        # plt.gca().set_aspect('equal', adjustable='box')
        # plt.legend()
        plt.xlabel(r'$x$')
        plt.ylabel(r'$t$')
        # plt.xticks([])
        # plt.yticks([])
        plt.grid(alpha=0.3)
        if file_name is not None and folder_name is not None:
            os.makedirs(folder_name, exist_ok=True)
            plt.savefig(folder_name + '/' + file_name + '_wave_struct.pdf', bbox_inches='tight')
    
    def saveSolution(self, folder_name, file_name):
        """
        Save the full object
        """
        os.makedirs(folder_name, exist_ok=True)
        full_path = folder_name+'/'+file_name+'.pik'
        with open(full_path, 'wb') as file:
            pickle.dump(self, file)
        print('Solution save to ' + full_path + ' !')



