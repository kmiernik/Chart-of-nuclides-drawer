#!/usr/bin/python3
"""
   Copyright Krzysztof Miernik 2012
   k.a.miernik@gmail.com 
   
   Distributed under GNU General Public Licence v3
   
   This Python script loads data from NuBase03 ascii file
   and creates xml document
   NuBase web page
   http://amdc.in2p3.fr/web/nubase_en.html

   NuBase2003 ascii format description 
   (strings are given in Python way i.e [0:3] are characters 0,1,2)
   based on G. Audi et al Nuc. Phys. A 729(2003) 3

   [0:3] mass number A
   [3] empty
   [4:7] atomic number Z
   [7] 0 for g.s or 1, 2, 3, etc. for isomers
   [8:11] empty
   [11:17] redundand - mass, element name and letter n, m, etc. for isomers
   [17] empty
   [18:35] mass defect and its uncertainity
   [35:39] empty
   [39:60] isomer excitation energy, uncertainity and code
   [60:78] half-life, unit uncertainity
   [78] empty
   [79:93] g.s. spin
   [93:105] year of apperance in nubase, reference code
   [105] empty
   [106:] decay modes

    Description of xml document structure:
    <!-- Root element -->
    <nuclear_data_table>

        <!-- Entry for nuclide, attributes:
        A - mass number,
        Z - atomic number,
        id - element name in form e.g 12C
        element - chemical element name -->
        <nuclide A="" Z="" id="">

            <!-- mass defect:
                 extrapolation is True if value is extrapolated (False otherwise)
                 uncertainity and value in keV -->
            <mass_defect extrapolation="False" uncertainity="" value=""/>

            <!-- half life:
                 extrapolation is True if value is extrapolated
                 uncertainity and value is in given units-->
            <half_life extrapolation="" uncertainity="" unit="" value=""/>

            <!-- Ground state spin 
                 extrapolation is True if value is extrapolated
                 value is total spin (J) of ground state-->
            <spin extrapolation="" value=""/>

            <!-- A list of decay modes-->
            <decay_modes>
                <!-- Each decay has a name (mode),
                     relation (=, ~, >, <, >=, <=)
                     and uncertainity -->
                <decay mode="" relation="" uncertainity="" value=""/>
                <decay .../>
            </decay_modes>

            <!-- A list of isomeric states -->
            <isomers>
                <!-- Each isomer has excitation energy in keV
                     extrapolation if energy is extrapolated
                     uncertainity in keV
                     half life, same as above
                     a list of decay modes, same as above
                     a comment including measurement method and reference -->
                <isomer energy="" extrapolation="" uncertainity="">
                    <half_life extrapolation="" uncertainity="" unit="" value=""/>
                    <decay_modes>
                    ...
                    </decay_modes>
                    <comment>...</comment>
                </isomer ...>
                ...
                </isomer>
            </isomers>

            <!-- Comment field includes year of apperance in nubase and
                 reference code -->
            <comment>...</comment>
        <nuclide>
        ...
        </nuclide>
        ...
    </nuclear_data_table>
"""

import argparse
import xml.dom.minidom
from Nuclide import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reads NuBase2003 ugly ascii file and writes xml document')
    parser.add_argument('infile', type=argparse.FileType('r'), 
                         help='Input data file (required)')
    parser.add_argument('outfile', type=argparse.FileType('w'), 
                         help='Output data file (required)')
    args = parser.parse_args()

    line_number = 0
    first = True
    data = []
    for line in args.infile:
        line_number += 1

        A = int(line[0:3])
        Z = int(line[4:7])
        isomer = False if line[7] == '0' else True
        mass_defect = line[18:38].strip()
        isomer_data = line[38:60].strip()
        half_life = line[60:78].strip()
        gs_spin = line[79:93].strip()
        comment = line[93:105].strip()
        decay_modes = line[106:-1].strip().lower()
        try:
            if not(isomer):
                if not(first):
                    data.append(isotope)
                else:
                    first = False
                isotope = NuclideNb03(Z, A, mass_defect, half_life,
                                      gs_spin, decay_modes, comment)
            else:
                isotope.nb_add_isomer(isomer_data, half_life,
                                      decay_modes, comment)

        except ParameterError as err:
            print('Line', line_number, ':', err.msg)
            if isomer:
                print('Skipping bad entry: isomer of isotope {}'.format(isotope))
            else:
                print('Skipping bad entry: isotope {}'.format(isotope))
            print()


    dom = xml.dom.minidom.getDOMImplementation()
    table = dom.createDocument(None, "nuclear_data_table", None)
    root = table.documentElement
    for isotope in data:
        isotope.add_to_xml_table(table, root)
    args.outfile.write(table.toprettyxml(indent="    ", encoding="utf-8").decode("utf-8"))

