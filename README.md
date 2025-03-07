# Generate JANAF Polynoms from CoolProp

This library generates CHEMKIN files with JANAF polynomials for a user given list of species and temperature range based on a given thermodynamic library backend. The default thermodynamic library is CoolProp. However, the code can be extended to use other libraries, e.g., NIST RefProp. 


## Motivation

OpenFOAM uses the CHEMKIN file format to read the reaction, species, and thermodynamic data for reactive flows. The CHEMKIN file reader only supports the gasHThermoPhysics thermo model which is:
```
    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > gasHThermoPhysics;
```
Hence, all thermodynamic data needs to be provided as JANAF polynomials according to Eqs. (19-21) in [1]. For the simulation of reacting flashing flows of cryogenic liquids, thermophysical data down to 70K must be known. 



## References

[1] R.J. Kee, F.M. Rupley, and J.A. Miller, "Chemkin-II: A Fortran Chemical
Kinetics Package for the Analysis of Gas-Phase Chemical Kinetics",
Sandia Report SAND89-8009.UC-401, September 1989


