# CHEMKIN II File Format Description

The Chemkin thermo file format is based very closely on the NASA (1971)
file format.  The specifics are documented in:

R.J. Kee, F.M. Rupley, and J.A. Miller, "Chemkin-II: A Fortran Chemical
Kinetics Package for the Analysis of Gas-Phase Chemical Kinetics",
Sandia Report SAND89-8009.UC-401, September 1989

Thermodynamic data are provided to the Chemkin linker (interpreter) either
within a mechanism file or as a separate file.  Thermo data provided in
the mechanism file overrides any data for the same species provided in a
separate thermo file.  The separate thermo file (if any) is typically
named "therm.dat".  In either case, the thermo data is expressed in the
same format.  Mechanism files generally consist of sections containing element
data, species data, thermo data (if not taken from a separate file), and
reaction rate data.

The rest of this description of Chemkin thermodynamic data format is quoted
directly from the above reference.  Since element and species data are
prerequisite to thermo data, those sections are described before the
thermo data format.  The reaction rate data format is not described here.

			Element data
All chemical species in the reaction mechanism must be composed of chemical
elements or isotopes of chemical elements.  Each element and isotope must
be declared as a one- or two-character symbol.  The purpose of the element
data is to associate atomic weights of the elements with their character
symbol representations and to identify the order in which arrays of element
information in the Gas-Phase Subroutine Library are referenced.  For example,
a Fortran array of atomic weights for the elements is in exactly the same
order in which the elements were declared in the element data.  In other
words, if the atomic weights are stored in array AWT, then AWT(3) is the
atomic weight of the third element declared in the element data.

For the elements appearing on the periodic chart, the Interpreter has the
atomic weight (in grams per mole) stored internally.  For isotopes, a one-
or two-character symbol must be input to the Interpreter to identify each
isotope, and a symbol and an atomic weight (in grams per mole) for each must
be defined.  The same symbol must be used in the thermodynamic data to
identify the elemental composition of species involving the isotope.  Once
an isotope has been so defined, it is treated exactly as a new element.  If
an ionic species is used in the mechanism (i.e. OH+), an electron must be
declared as the element E.

Element data must start with the word ELEMENTS (or ELEM), followed by any
number of element symbols on any number of lines.  Element symbols may appear
anywhere on a line, but those on the same line must be separated by blanks.
Any line or portion of a line starting with an exclamation mark (!) is
considered a comment and will be ignored.  Blank lines are ignored.

If an element is on the periodic chart, then only the symbol identifying the
element need appear in the element data.  [Elements recognized automatically
by Chemkin are listed as H through Fm, plus E and D]  For an isotope,
the atomic weight
must follow the identifying symbol and be delimited by slashes (/).  The
atomic weight may be an integer, floating point, or E format (e.g. 2, 2.0,
0.2E1), but internally will be converted to a floating point number.  For
example, the isotope deuterium may be defined as D/2.014/.  If desired, the
atomic weight of an element in the periodic chart may be altered by including
the atomic weight as input just as though the element were an isotope.

Figure 4 shows several equivalent wasy to describe element information.  In
this example the elements are hydrogen, oxygen, nitrogen, and the isotope
deuterium.  Table I summarizes the rules for element data.

ELEMENTS   H D/2.014 O N END

ELEM           ! ELEM is equivalent to ELEMENTS
 H
 D /2.014/
 O
 N
END

ELEM H
ELEM D/2.014/
ELEM O
ELEM N
	Figure 4. Equivalent Ways to Describe Element Information.

	Table 1. Summary of the Rules for Element Data:
1. The first element line must start with the word ELEMENTS (or ELEM).
2. Element or isotope names are either one- or two-character symbols.
3. An isotope name (i.e., a name not on the periodic chart) must be followed
   by its atomic weight (in grams per mole) delimited by slashes.
4. Each element or isotope should be declared only once; however, duplicate
   element symbols are ignored.
5. An element or isotope name may appear anywhere on the line.
6. Any number of element or isotope names may appear on a line, and more than
   one line may be used.
7. Element or isotope names that appear on the same line must be separated by
   at least one blank space.
8. An element or isotope name that begins on one line may not continue to the
   next line.
9. Any blank spaces between an element or isotope name and the first slash
   are ignored and any blank spaces between slashes and an atomic weight are
   also ignored.  However, no blank spaces are allowed within an element name
   or an atomic weight.
10. There may be more than one ELEMENT statement.
11. All characters following an exclamation mark are comments.

			Species Data
Each chemical species in a problem must be identified on a species line(s).
Any set of up to 16 upper or lower case characters* can be used as a species
name.  In addition, each species must be composed of elements that have been
identified in the element data.  As for the element data, one of the primary
purposes of the species data is to identify the order in which Fortran arrays
of species information are referenced in the Gas-Phase Subroutine Library.

Species data must start with the word SPECIES (or SPEC), followed by any
number of species symbols on any number of lines.  Species symbols may appear
anywhere on a line, but those on the same line must be separated by blank
spaces.  Any line or portion of a line starting with an exclamation mark (!)
is considered a comment and will be ignored.  Blank lines are ignored.
Figure 5 shows several equivalent ways to describe species information.  The
rules for species data are summarized in Table 2.

SPECIES  H2 O2 H O OH HO2 N2 N NO END

SPEC          ! SPEC is equivalent to SPECIES
 H2 O2
 H O OH HO2 N2 N NO
END

SPEC H2
spec O2
etc.

	Figure 5. Equivalent Ways to Describe Species Information.
*Species names may not begin with a number, a +, or an =, or have imbedded
blanks; an ionic species may end with any number of +'s or -'s; an imbedded
plus sign (+) must be enclosed in parentheses.


	Table 2. Summary of the Rules for Species Data
1. Species data must start with the word SPECIES (or SPEC).
2. Species names are composed of up to 16-character upper- or lower-case
   symbols.  The names cannot begin with the characters +, =, or a number;
   an ionic species name may end with one or more +'s or -'s.
3. Each species should be declared only once; however, duplicated species
   symbols will be ignored.
4. Each species that subsequently appears in a reaction must be declared.
5. A species name may appear anywhere on a line.
6. Any number of species names may appear on a line, and more than one line
   may be used.
7. Species named on the same line must be separated by at least one blank
   space.
8. A species name that begins on one line may not continue to the next line.
9. There may be more than one SPECIES statement.
10. All characters following an exclamation mark are comments.

			Thermodynamic Data
Any chemical species that appears in a problem must have thermodynamic data
associated with it. The data may be extracted from a database (file LTHRM)
and/or read from input (file LIN).  If all the thermodynamic data are to be
extracted from the database, then no thermodynamic data input is required.
However, if the user wishes to override information in the database or to
provide data on species not in the database, then Interpreter input is needed.
In any case the format for the information is the same.

The format (see Table 3) is a minor modification of that used by Gordon and
McBride for the Thermodynamic Database in the NASA Chemical Equilibrium
code.  Our modification allows for a different midpoint temperature for
the fits to the properties of each chemical species.  We also allow a
species to be composed of a maximum of five elements, not four.  However,
the formatting is such that the Chemkin Interpreter can use the NASA
database directly without any modification.

As indicated in Table 3, the pertinent information includes the species name,
the elemental composition of the species, and the temperature ranges over
which the polynomial fits to thermodynamic data are valid.  The fits to
Cp0/R, H0/RT, and S0/R consist of seven coefficients for each of two
temperature ranges.  Further information about the fitting procedure and
data for many species can be found in a report on the Chemkin Thermodynamic
Database.

When thermodynamic data input is required, it must immediately follow
species data.  The first thermodynamic data line must start with the
word THERMO (or THER).  If all the thermodynamic data are input directly
to the Interpreter, then the first line must read THERMO ALL and the
code will not expect a Thermodynamic Database from file LTHRM; for this
option the next line must be line 2 of Table 3.  For either option, the
subsequent thermodynamic data lines must be in the format of lines 3-6
of Table 3.  (For the THERMO option the midpoint temperature is taken
from the line 2 information already in the Thermodynamic Database.)  As
many species as needed can be included as THERMO input.

Figure 6 shows some examples of thermodynamic property input.  In these
three examples for OH, OH+, and OH-, it is seen from columns 25-34 that
the elemental composition of each molecule is one O atom and one H atom.
In addition, columns 35-39 indicate that two of the species, OH+ and OH-,
are ionic since they contain -1 and +1 electrons (E), respectively.  The
G in column 45 indicates that all three species are gaseous.  The 1000.00
in columns 66-73 for OH+ indicates that the common temperature between
the high- and low-temperature fits is 1000.00 K.  If columns 66-73 are
left blank, as they are for OH+ and OH-, then the common temperature is that
given in columns 21-30 of line 2 in Table 3, which in this example is in the
Thermodynamic Database.

	Table 3. Summary of the Rules for Thermo Data
Line   Contents                                         Format     Column
Number
1      THERMO (or THERMO ALL*)                          Free       Any
2#     Temperature ranges for 2 sets of coefficients:   3F10.0     1 to 30
       lowest T, common T, and highest T
3      Species name (must start in column 1)            18A1       1 to 18
       Date (not used in the code)                      6A1        19 to 24
       Atomic symbols and formula                       4(2A1,I3)  25 to 44
       Phase of species (S, L, or G for solid,          A1         45
       liquid, or gas, respectively) 
       Low temperature                                  E10.0      46 to 55
       High temperature                                 E10.0      56 to 65
       Common temperature (if needed)                   E8.0       66 to 73
       (blank for default)
       Atomic symbols and formula (if needed)           2A1,I3     74 to 78
       (blank for default)
       The integer 1                                    I1         80
4      Coefficients a1-a5 in Eqs. (19) - (21),          5(E15.0)   1 to 75
       for upper temperature interval
       The integer 2
5      Coefficients a6, a7 for upper temperature        5(E15.0)   1 to 75
       interval, and a1, a2, and a3 for lower
       The integer 3                                    I1         80
6      Coefficients a4, a5, a6, a7                      4(E15.0)   1 to 60
       for lower temperature interval          
       The integer 4                                    I1         80
...    Repeat lines 3-6 for each species.
last   END (Optional, end of thermodynamic data)        Free       Any

* Use only when all the thermodynamic data are to be taken from Interpreter
 input.
# Include line 2 only with THERMO ALL (it is already in the Thermodynamic
 Database).

The following is a summary of the possibilities for specifying thermodynamic
data.
Case 1: All thermodynamic data from database
        1. Assign the database as file LTHRM (default Fortran unit 17)
        2. No THERMO data required as input
Case 2: Thermodynamic data from database and input
        1. Assign the database as file LTHRM (default Fortran unit 17)
        2. Include the following lines:
           THERMO
           Data in Table 3 format (lines 3-6 repeated) for species not in
           the database or to override species in database
           END
Case 3: All thermodynamic data from input
        1. Do not attach a database
        2. Include the following lines:
           THERMO ALL
           Line 2 of Table 3 format.
           Data in Table 3 format (lines 3-6 repeated) for at least all
           species named in the species data.
           END

Example (top of file):
THERMO
   300.000  1000.000  5000.000
(CH2O)3            70590C   3H   6O   3     G  0300.00   4000.00  1500.00      1
 0.01913678E+03 0.08578044E-01-0.08882060E-05-0.03574819E-08 0.06605143E-12    2
-0.06560876E+06-0.08432507E+03-0.04662286E+02 0.06091547E+00-0.04710536E-03    3
 0.01968843E-06-0.03563271E-10-0.05665404E+06 0.04525265E+03                   4

Notes:
1. The most significant difference between the Chemkin thermo file format
 and the NASA (1971) style is the use of columns 66 to 78 for a
 species-specific common (high and low) temperature and an additional element
 specification.  The NASA format does not specify the use of this area,
 although many NASA format files place the molecular weight here.
 Also, the NASA format does not allow for different common temperatures
 for the high and low temperature polynomials for each species.
2. Because of the differences (plus one or two other subtle differences),
 Chemkin thermo files are, in general, not compatible with NASA (1971)
 thermo files, although this depends on the specific usage.  In many
 cases, the two are compatible.
