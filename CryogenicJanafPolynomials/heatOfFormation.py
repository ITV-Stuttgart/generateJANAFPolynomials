"""Dictionary with the heat of formation"""
import importlib.resources


class HeatOfFormation:

    def _getValidName(self,specie):
        # If the specie is not equal the short name, look through the long names
        if not specie in self._specieNameTable:
            for key in self._specieNameTable:
                if specie in self._specieNameTable[key]:
                    return key
        return specie


    def __init__(self):
        # Read the heatOfFormation.dat file
        self._heatOfFormationFile = 'heatOfFormation.dat'
        self._data = {}
        self._specieNameTable = {}
        self.load_data()

    def getStandard(self,specie):
        specie = self._getValidName(specie)
        if specie in self._data:
            return self._data[specie]['heatOfFormation_298K']
        else:
            print(f"Specie {specie} not found in data bank")
            return 0

    def load_data(self):
        """Reads the static heat of formation file inside the package."""
        try:
            # Open the file from the package
            with importlib.resources.open_text(__package__, self._heatOfFormationFile) as f:
                for line in f:
                    lineParts = line.split(';')
                    speciesShortName= lineParts[0]
                    prop = {'speciesFullName': lineParts[1],
                            'heatOfFormation_0K': float(lineParts[2])*1000.0,   # Convert to SI unit [J/mol]
                            'heatOfFormation_298K': float(lineParts[3])*1000.0  # Convert to SI unit [J/mol]
                            }
                    self._specieNameTable[speciesShortName]=[lineParts[1]]
                    self._data[speciesShortName] = prop

        except FileNotFoundError:
            print("heatOfFormation.dat file not found!")

