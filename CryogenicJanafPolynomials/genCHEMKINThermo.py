"""Generates the CHEMKIN thermo file for a given set of species"""

import numpy as np
from .thermoFit import ThermoFit


class GenCHEMKINThermo:

    def __init__(self, TMin, TCommon, TMax, p):
        """
        Parameters
        ----------
        TMin, TCommon, TMax : float
            Minimum, maximum and common temperature for the lower and upper
            polynomial fit

        p : float
            Pressure value to evaluate thermodynamic properties


        Species Data Dictionary
        -----------------------    
        Dictionary with the species data. Contains the polynomial fits and
        species molar mass.
        Keys: 
            - lowerPolyFit   : Lower polynomial fit
            - upperPolyFit   : Upper polynomial fit
            - date           : Date format of CHEMKIN-II
            - atomicSymbol   : Atomic symbols and formula
            - phase          : Phase of species 
        Access e.g., the lower polynomial fit
          self.speciesData[speciesName][lowerPolyFit] 
        """
        self.TMin = TMin
        self.TCommon = TCommon
        self.TMax = TMax
        self.p = p
        self.species = []
        self.speciesData = {}

    def genThermoData(self,*args):
        if len(args) == 0:
            species = self.species
        else:
            species = args[0]

        for specie in species:
            if ThermoFit.ContainsFluid(specie):
                # First for the lower temperature range
                lowerFit = ThermoFit(specie, self.TMin, self.TCommon, p=self.p)
                if lowerFit.genFromThermoLib():
                    self.speciesData[specie]['lowerPolyFit'] = lowerFit

                upperFit = ThermoFit(specie, self.TCommon, self.TMax, p=self.p)
                if upperFit.genFromThermoLib():
                    self.speciesData[specie]['upperPolyFit'] = upperFit

    def writeThermoFile(self, filename):
        with open(filename, "w") as fp:
            # First write a comment block
            fp.write('! This file has been generated with the CryogenicJanafPolynmials tool\n')
            fp.write('! Polynomial coefficients are generated from the CoolProp library\n')
            fp.write('! The Min/Max temperature is set to the global values,\n')
            fp.write('! some polynomials may result in unpyhsical values.\n')
            fp.write('! However, this is required for OpenFOAM to not limit the temperature range\n')
            fp.write('! to the largest TMin or lowest TMax\n')
            fp.write("THERMO ALL\n")
            fp.write(
                f' {self.TMin: 09.3f} {self.TCommon: 09.3f} {self.TMax: 09.3f}\n')
            for specie in self.species:
                data = self.speciesData[specie]
                if 'lowerPolyFit' in data:
                    if data['lowerPolyFit'].dataFromThermoLib or data['upperPolyFit'].dataFromThermoLib:
                        fp.write(f"! {specie} generated with the polynomial fit with CoolProp data\n")
                    lowerCoeff = data['lowerPolyFit'].coeffs
                    upperCoeff = data['upperPolyFit'].coeffs
                    fp.write(
                        f"{specie: <18}{data['date']: <6}{data['atomicSymbol']: <19}{data['phase']}{self.TMin: 10.3f}{self.TMax: 10.3f}{self.TCommon: 010.3f}    1\n")
                    fp.write(
                        f"{lowerCoeff[0]: 15E}{lowerCoeff[1]: 15E}{lowerCoeff[2]: 15E}{lowerCoeff[3]: 15E}{lowerCoeff[4]: 15E}    2\n")
                    fp.write(
                        f"{lowerCoeff[5]: 15E}{lowerCoeff[6]: 15E}{upperCoeff[0]: 15E}{upperCoeff[1]: 15E}{upperCoeff[2]: 15E}    3\n")
                    fp.write(
                        f"{upperCoeff[3]: 15E}{upperCoeff[4]: 15E}{upperCoeff[5]: 15E}{upperCoeff[6]: 15E}                   4\n")
                else:
                    print(f"No thermo data is available for {specie}!")
            fp.write("END\n")

    def readChemkinFile(self, chemkinFile):
        """Read in an existing therm.dat CHEMKIN-II file to get the list of 
        species.
        """
        readTemperatureRange = False
        readSpeciesLine = -1
        TMin = 1
        TMax = 1
        TCommon = 1
        currSpeciesName = "undefined"
        with open(chemkinFile) as fp:
            for line in fp:
                # Check for comments
                if line[0] == '!':
                    continue
                if "END" in line:
                    break

                if "THERMO" in line:
                    # Read the temperature range in the next line
                    readTemperatureRange = True
                    continue

                if readTemperatureRange:
                    lineParts = line.split()
                    TMin = float(lineParts[0])
                    TCommon = float(lineParts[1])
                    TMax = float(lineParts[2])
                    readTemperatureRange = False
                    readSpeciesLine = 0
                    continue

                if readSpeciesLine == 0:
                    specieName = line[0:18]
                    specieName = specieName.strip()
                    self.species.append(specieName)
                    self.speciesData[specieName] = {}
                    self.speciesData[specieName]['date'] = line[18:24]
                    self.speciesData[specieName]['atomicSymbol'] = line[24:44]
                    self.speciesData[specieName]['phase'] = line[44]
                    
                    TRange = line[45:76].split()
                    TMin = float(TRange[0])
                    TMax = float(TRange[1])
                    TCommon = float(TRange[2])

                    # Indicate that it is read from the file
                    self.speciesData[specieName]['dataFromThermoLib'] = False
                    currSpeciesName = specieName
                    readSpeciesLine += 1
                    readSpeciesLine %= 4
                elif readSpeciesLine == 1:
                    coeff = np.zeros(7)
                    coeff[0] = float(line[0:15])
                    coeff[1] = float(line[15:30])
                    coeff[2] = float(line[30:45])
                    coeff[3] = float(line[45:60])
                    coeff[4] = float(line[60:75])
                    self.speciesData[currSpeciesName]['upperPolyFit'] = ThermoFit(
                        currSpeciesName, TCommon, TMax, phase=self.speciesData[currSpeciesName]['phase'], coeffs=coeff)
                    readSpeciesLine += 1
                    readSpeciesLine %= 4
                elif readSpeciesLine == 2:
                    upperCoeffFit = self.speciesData[currSpeciesName]['upperPolyFit']
                    lowerCoeff = np.zeros(7)
                    upperCoeffFit.coeffs[5] = float(line[0:15])
                    upperCoeffFit.coeffs[6] = float(line[15:30])
                    lowerCoeff[0] = float(line[30:45])
                    lowerCoeff[1] = float(line[45:60])
                    lowerCoeff[2] = float(line[60:75])
                    self.speciesData[currSpeciesName]['upperPolyFit'] = upperCoeffFit
                    self.speciesData[currSpeciesName]['lowerPolyFit'] = ThermoFit(
                        currSpeciesName, TMin, TCommon, phase=self.speciesData[currSpeciesName]['phase'], coeffs=lowerCoeff)
                    readSpeciesLine += 1
                    readSpeciesLine %= 4
                elif readSpeciesLine == 3:
                    lowerCoeffFit = self.speciesData[currSpeciesName]['lowerPolyFit']
                    lowerCoeffFit.coeffs[3] = float(line[0:15])
                    lowerCoeffFit.coeffs[4] = float(line[15:30])
                    lowerCoeffFit.coeffs[5] = float(line[30:45])
                    lowerCoeffFit.coeffs[6] = float(line[45:60])
                    self.speciesData[currSpeciesName]['lowerPolyFit'] = lowerCoeffFit
                    readSpeciesLine += 1
                    readSpeciesLine %= 4

    def plot(self, ax, speciesName, var, **kwargs):
        """Plot the species property in the given axis"""

        plot_kwargs = kwargs.copy()

        l1 = self.speciesData[speciesName]['lowerPolyFit'].plot(ax,var,**plot_kwargs)
        plot_kwargs.pop("label", None)
        plot_kwargs.pop("color", None)
        self.speciesData[speciesName]['upperPolyFit'].plot(ax,var,color=l1[0].get_color(),**plot_kwargs)
        return ax

    def dataFromThermoLib(self,specieName):
        if self.speciesData[specieName]['lowerPolyFit'].dataFromThermoLib or self.speciesData[specieName]['upperPolyFit'].dataFromThermoLib:
            return True
        return False
