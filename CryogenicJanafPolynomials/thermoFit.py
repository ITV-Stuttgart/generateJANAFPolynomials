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
from numpy.polynomial import Polynomial

class ThermoFit:
    _thermoLibraryBackend="CoolProp"
    _atm=101325                         # One atmospheric pressure
    _fluidNameFull = []
    _fluidNameAliases = {}
    
    def _genThermoData(self):
        """Generate the Cp, H, and S values for the species in the temperature
        range and return the T vector, Cp, H, and S vectors
        """
        # Generate the temperature vector
        T = np.linspace(self.TMin,self.TMax,100)

        # Calculate the enthalpy, entropy and heat capacity values
        H = np.zeros(T.shape)
        Cp = np.zeros(T.shape)
        S = np.zeros(T.shape)

        for i in range(len(T)):
            try:
                H[i] = CP.PropsSI('H','P|'+self.phase,self.p,'T',T[i],self.speciesName)
            except:
                H[i] = np.nan
            try:
                Cp[i] = CP.PropsSI('C','P|'+self.phase,self.p,'T',T[i],self.speciesName)
            except:
                Cp[i] = np.nan
        
        return T,Cp,H,S


   
    def _polyFit(self):
        """Create the polynomial fit to a given data set"""
        T,Cp,H,S = self._genThermoData()
        # Generate the polynomial fit for the heat capacity Eq. (19)
        CpFit = Polynomial.fit(T,Cp,4)
        CpCoeffs = CpFit.convert().coef

        # Divide the coefficients by the specific gas constant to match
        # the description of Eq. (19)
        CpCoeffs /=  8.314/CP.PropsSI('M','P|'+self.phase,self.p,'T',298,self.speciesName)
        return CpFit, CpCoeffs
        
        

    def __init__(self,speciesName,TMin,TMax,**kwargs):
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

        """
        self.speciesName = speciesName
        self.TMin = TMin
        self.TMax = TMax
        self.phase = 'gas'
        self.p = self._atm
        self.valid=False

        if 'phase' in kwargs:
            self.phase = kwargs['phase']
        
        if 'p' in kwargs:
            self.p = kwargs['p']

        self.CpFit, CpCoeffs = self._polyFit()
        self.coeffs = np.zeros(7)
        self.coeffs[0:len(CpCoeffs)] = CpCoeffs
        if np.count_nonzero(np.isnan(CpCoeffs)) > 0:
            self.valid = False
        else:
            self.valid = True

    @classmethod
    def ContainsFluid(cls,fluidName):
        if not cls._fluidNameFull:
            cls._fluidNameFull = CP.FluidsList()
            for fluid in cls._fluidNameFull:
                cls._fluidNameAliases[fluid] = CP.get_aliases(fluid)

        # This is a bit ridiculus but I do not know a better way now
        # Loop through all aliases
        for key, value in cls._fluidNameAliases.items():
            if fluidName in value:
                return True,key

        return False,fluidName

    def plot(self,ax, **kwargs):
        xx, yy = self.CpFit.linspace()
        ax.plot(xx,yy,**kwargs)










