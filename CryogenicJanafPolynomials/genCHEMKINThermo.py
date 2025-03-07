"""Generates the CHEMKIN thermo file for a given set of species"""

import numpy as np
from .thermoFit import ThermoFit


class GenCHEMKINThermo:

    def __init__(self,TMin,TCommon,TMax,p):
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
            - lowerPolyCoeff   : Lower polynomial fit
            - upperPolyCoeff   : Upper polynomial fit
            - date           : Date format of CHEMKIN-II
            - atomicSymbol   : Atomic symbols and formula
            - phase          : Phase of species 
        Access e.g., the lower polynomial fit
          self.speciesData[speciesName][lowerPolyCoeff] 
        """
        self.TMin = TMin
        self.TCommon = TCommon
        self.TMax = TMax
        self.p = p
        self.species = []
        self.speciesData = {}

    def genThermoData(self):

        for specie in self.species:
            foundSpecie,specieString = ThermoFit.ContainsFluid(specie)
            if foundSpecie:
                # First for the lower temperature range
                lowerFit = ThermoFit(specieString, self.TMin, 1.1*self.TCommon,p=self.p)
                if lowerFit.valid:
                    self.speciesData[specie]['validPolyFit'] = True
                    self.speciesData[specie]['lowerPolyCoeff'] = lowerFit.coeffs

                upperFit = ThermoFit(specieString, 0.9*self.TCommon, 1.1*self.TMax,p=self.p)
                if upperFit.valid:
                    self.speciesData[specie]['upperPolyCoeff'] = upperFit.coeffs

    def writeThermoFile(self):
        with open("thermNew.dat", "w") as fp:
            # First write a comment block
            fp.write('! This file has been generated with the CryogenicJanafPolynmials tool\n')
            fp.write('! Polynomial coefficients are generated from the CoolProp library\n')

            fp.write("THERMO ALL\n")
            fp.write(f' {self.TMin: 09.3f} {self.TCommon: 09.3f} {self.TMax: 09.3f}\n')
            for specie in self.species:
                data = self.speciesData[specie]
                fp.write(f"{specie: <18}{data['date']: <6}{data['atomicSymbol']: <19}{data['phase']}{self.TMin: 10.3f}{self.TMax: 10.3f}{self.TCommon: 010.3f}    1\n")
                
                if 'lowerPolyCoeff' in data:
                    if data['validPolyFit']:
                        fp.write("! Generated with the polynomial fit with CoolProp data\n")
                    lowerCoeff = data['lowerPolyCoeff']
                    upperCoeff = data['upperPolyCoeff']
                    fp.write(f"{lowerCoeff[0]: 15E}{lowerCoeff[1]: 15E}{lowerCoeff[2]: 15E}{lowerCoeff[3]: 15E}{lowerCoeff[4]: 15E}    2\n")
                    fp.write(f"{lowerCoeff[5]: 15E}{lowerCoeff[6]: 15E}{upperCoeff[0]: 15E}{upperCoeff[1]: 15E}{upperCoeff[2]: 15E}    3\n")
                    fp.write(f"{upperCoeff[3]: 15E}{upperCoeff[4]: 15E}{upperCoeff[5]: 15E}{upperCoeff[6]: 15E}                   4\n")
                else:
                    print(f"No thermo data is available for {specie}!")
            fp.write("END\n")
                


    def readChemkinFile(self,chemkinFile):
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

                    # Indicate that it is read from the file
                    self.speciesData[specieName]['validPolyFit'] = False
                    currSpeciesName = specieName
                    readSpeciesLine += 1
                    readSpeciesLine %= 4
                elif readSpeciesLine == 1:
                    coeff = np.zeros(7)
                    coeff[0] = float(line[ 0:15])
                    coeff[1] = float(line[15:30])
                    coeff[2] = float(line[30:45])
                    coeff[3] = float(line[45:60])
                    coeff[4] = float(line[60:75])
                    self.speciesData[currSpeciesName]['upperPolyCoeff'] = coeff
                    readSpeciesLine += 1
                    readSpeciesLine %= 4
                elif readSpeciesLine == 2:
                    upperCoeff = self.speciesData[currSpeciesName]['upperPolyCoeff']
                    lowerCoeff = np.zeros(7)
                    upperCoeff[5] = float(line[ 0:15])
                    upperCoeff[6] = float(line[15:30])
                    lowerCoeff[0] = float(line[30:45])
                    lowerCoeff[1] = float(line[45:60])
                    lowerCoeff[2] = float(line[60:75])
                    self.speciesData[currSpeciesName]['upperPolyCoeff'] = upperCoeff
                    self.speciesData[currSpeciesName]['lowerPolyCoeff'] = lowerCoeff
                    readSpeciesLine += 1
                    readSpeciesLine %= 4
                elif readSpeciesLine == 3:
                    lowerCoeff = self.speciesData[currSpeciesName]['lowerPolyCoeff']
                    lowerCoeff[3] = float(line[ 0:15])
                    lowerCoeff[4] = float(line[15:30])
                    lowerCoeff[5] = float(line[30:45])
                    lowerCoeff[6] = float(line[45:60])
                    self.speciesData[currSpeciesName]['lowerPolyCoeff'] = lowerCoeff
                    readSpeciesLine += 1
                    readSpeciesLine %= 4

    def plot(self,ax,speciesName,var,**kwargs):
        def CpPoly(T,CpCoeffs):
            val = 0
            for i in range(5):
                val = val + CpCoeffs[i]*np.power(T,i)
            return val
        
        if var == "Cp":
            TLow = np.linspace(self.TMin,self.TCommon)
            l1, = ax.plot(TLow,CpPoly(TLow,self.speciesData[speciesName]['lowerPolyCoeff']),**kwargs)
            THigh = np.linspace(self.TCommon,self.TMax)
            l1, = ax.plot(THigh,CpPoly(THigh,self.speciesData[speciesName]['upperPolyCoeff']),color=l1.get_color(),**kwargs)
        return ax


                    



                




                


