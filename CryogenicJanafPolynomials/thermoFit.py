"""ThermoFit Class
This class calculates the polynomial coefficients for the fit of the heat
capacity (Cp), entropy (S) and enthalpy (H). See also Eqs. (19-21) in [1]

References:
[1] R.J. Kee, F.M. Rupley, and J.A. Miller, "Chemkin-II: A Fortran Chemical
    Kinetics Package for the Analysis of Gas-Phase Chemical Kinetics",
    Sandia Report SAND89-8009.UC-401, September 1989
"""

import CoolProp.CoolProp as CP
import numpy as np
import scipy.optimize as opt
from scipy.optimize import curve_fit
from .heatOfFormation import HeatOfFormation


def _CpFunc(T, coeffs):
    """Model function for the heat capacity Cp given in Eq. (19) of [1]"""
    return coeffs[0] + coeffs[1]*T + coeffs[2]*np.power(T, 2) + coeffs[3]*np.power(T, 3) + coeffs[4]*np.power(T, 4)

# Define the residuals for least squares fitting


def _CpObjective(T, Cp, coeffs):
    y_pred = _CpFunc(T, coeffs)
    return np.sum((y_pred - Cp)**2)


def _CpFuncDeriv(T, coeffs):
    """Derivative of the heat capacity function"""
    return coeffs[1] + 2 * coeffs[2] * T + 3 * coeffs[3] * np.power(T, 2) + 4 * coeffs[4] * np.power(T, 3)


def _HFunc2(T, CpCoeffs, a6):
    """Model function for the enthalpy given in Eq. (20) of [1]"""
    return CpCoeffs[0] + CpCoeffs[1]/2.0*T + CpCoeffs[2]/3.0*np.power(T, 2) + CpCoeffs[3]/4.0*np.power(T, 3) + CpCoeffs[4]/5.0*np.power(T, 4) + a6/T


def _HFunc(T, coeffs):
    """Model function for the enthalpy given in Eq. (20) of [1]"""
    return coeffs[0] + coeffs[1]/2.0*T + coeffs[2]/3.0*np.power(T, 2) + coeffs[3]/4.0*np.power(T, 3) + coeffs[4]/5.0*np.power(T, 4) + coeffs[5]/T


class ThermoFit:
    _thermoLibraryBackend = "CoolProp"
    _atm = 101325                         # One atmospheric pressure
    _fluidNameAliases = {}
    _heatOfFormation = HeatOfFormation()

    def _genThermoData(self):
        """Generate the Cp, H, and S values for the species in the temperature
        range and return the T vector, Cp, H, and S vectors
        """
        # Generate the temperature vector
        T = np.linspace(0.9*self.TMin, 1.1*self.TMax, 120)

        # Calculate the enthalpy, entropy and heat capacity values
        H = np.zeros(T.shape)
        Cp = np.zeros(T.shape)
        S = np.zeros(T.shape)

        for i in range(len(T)):
            try:
                H[i] = CP.PropsSI('H', 'P|'+self.phase, self.p, 'T', T[i], self.ThermoFluidName(self.speciesName))
                Cp[i] = CP.PropsSI('C', 'P|'+self.phase, self.p, 'T', T[i], self.ThermoFluidName(self.speciesName))
            except:
                H[i] = np.nan
                Cp[i] = np.nan

        # Check that the Cp value is physical -- non NaN and montone increasing
        # Remove all NaN entries
        ind = np.argwhere(~np.isnan(Cp)).flatten()
        Cp = Cp[ind]
        H = H[ind]
        T = T[ind]

        # Find minimum value
        ind = np.argmin(Cp)

        if ind > 0:
            # Calculate the first valid gradient after the ind position.
            # It may happen that the gradient at T[ind] is negative, which is unphysical.
            # Then higher temperature points must be found
            CpGrad = CP.PropsSI('d(C)/d(T)|P', 'P|'+self.phase, self.p, 'T',
                                T[ind], self.ThermoFluidName(self.speciesName))
            if CpGrad <= 0:
                ind += 1
                CpGrad = CP.PropsSI('d(C)/d(T)|P', 'P|'+self.phase, self.p, 'T',
                                    T[ind], self.ThermoFluidName(self.speciesName))
            for i in range(ind):
                Cp[i] = (T[i]-T[ind])*CpGrad + Cp[ind]

            # Update TMin
            self.TMin = T[ind]

        # Set the reference temperature to 298K
        H0 = CP.PropsSI('H', 'P|'+self.phase, self.p, 'T', 298.15, self.ThermoFluidName(self.speciesName))
        H = H - H0

        # Add the heat of formation
        M= CP.PropsSI('M', 'P|'+self.phase, self.p, 'T', 298, self.ThermoFluidName(self.speciesName))
        H = H + self._heatOfFormation.getStandard(self.speciesName)/M

        return T, Cp, H, S

    def __init__(self, speciesName, TMin, TMax, **kwargs):
        """
        Parameters
        ----------
        speciesName : string
            Name of the species

        TMin, TMax : float
            Minimum and maximum of the temperature range

        Keywords
        --------
        phase : string
            Phase of the species by default it is gas

        p : float
            Refernce pressure level for the calculation of thermodynamic 
            properties. By default set to 1 atm.

        coeffs : *float
            Array of floats of the 7 coefficients

        """
        self.speciesName = speciesName
        self.TMin = TMin
        self.TMax = TMax
        self.phase = 'gas'
        self.p = self._atm
        self.dataFromThermoLib = False
        self.coeffs = np.zeros(7)

        if 'phase' in kwargs:
            self.phase = kwargs['phase']

        if 'p' in kwargs:
            self.p = kwargs['p']

        if 'coeffs' in kwargs:
            self.coeffs = kwargs['coeffs']

    @classmethod
    def ContainsFluid(cls, fluidName):
        if not cls._fluidNameAliases:
            for fluid in CP.FluidsList():
                cls._fluidNameAliases[fluid] = CP.get_aliases(fluid)

        # This is a bit ridiculus but I do not know a better way now
        # Loop through all aliases
        for key, value in cls._fluidNameAliases.items():
            if fluidName in value:
                return True

        return False

    @classmethod
    def ThermoFluidName(cls, fluidName):
        """Convert the given fluid name to the accepted name for the thermo
        library
        """
        if not cls._fluidNameAliases:
            for fluid in CP.FluidsList():
                cls._fluidNameAliases[fluid] = CP.get_aliases(fluid)

        # This is a bit ridiculus but I do not know a better way now
        # Loop through all aliases
        for key, value in cls._fluidNameAliases.items():
            if fluidName in value:
                return key

        # if not found return the fluid name
        return fluidName

    def genFromThermoLib(self):
        """Create the polynomial fit to a given data set"""
        T, Cp, H, S = self._genThermoData()

        # Divide Cp by the specific gas constant to match
        # the description of Eq. (19) [1]
        RGas = (8.314/CP.PropsSI('M', 'P|'+self.phase, self.p, 'T', 298, self.ThermoFluidName(self.speciesName)))

        CpScaled = Cp/RGas

        # Create an initial guess for the Cp polynomial
        # Result is ordered with highest degree first
        Cp_polyFit = np.polyfit(T, CpScaled, deg=4)
        Cp_polyFit = np.flip(Cp_polyFit)

        # Define the constraint for monotonic increase
        constraint = {'type': 'ineq', 'fun': lambda coeffs: _CpFuncDeriv(T, coeffs)}

        CpCoeffs = opt.minimize(lambda coeffs: _CpObjective(T, CpScaled, coeffs), Cp_polyFit, constraints=[constraint])
        CpCoeffs = CpCoeffs.x

        # Fit the last coefficient of H

        HScaled = H/(RGas*T)

        # Create an initial guess for the Cp polynomial
        # Result is ordered with highest degree first
        popt, _ = curve_fit(lambda T, a6: _HFunc2(T, CpCoeffs, a6), T, HScaled)

        self.coeffs = np.zeros(7)
        self.coeffs[0:5] = CpCoeffs
        self.coeffs[5] = popt
        if np.count_nonzero(np.isnan(self.coeffs)) > 0:
            self.dataFromThermoLib = False
        else:
            self.dataFromThermoLib = True

        return self.dataFromThermoLib

    def plot(self, ax, prop, **kwargs):
        xMin = self.TMin
        xMax = self.TMax
        if 'TMin' in kwargs:
            xMin = kwargs['TMin']
        if 'TMax' in kwargs:
            xMin = kwargs['TMax']

        xx = np.linspace(xMin, xMax)
        if prop == "Cp":
            return ax.plot(xx, _CpFunc(xx, self.coeffs), **kwargs)
        elif prop == "H":
            return ax.plot(xx, _HFunc(xx, self.coeffs), **kwargs)
