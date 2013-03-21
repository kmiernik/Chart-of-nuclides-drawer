#!/usr/bin/python3
"""
   Copyright Krzysztof Miernik 2012
   k.a.miernik@gmail.com 
   
   Distributed under GNU General Public Licence v3
   
   This Python script loads data from Nuclear Wallet Cards ascii file
   and creates xml document
   http://www.nndc.bnl.gov/wallet/
   Nuclear Waller Cards by J.K. Tulli, National Nuclear Data Center,
   Brookhaven National Laboratory

   NWC ascii format description 
   (strings are given in Python way i.e [0:3] are characters 0,1,2)
   [0] 'F' for 235U product
   [1:4] mass number A
   [4] M for isomer
   [6:9] atomic number Z
   [10:12] element
   [16:26] spin
   [30:34] decay mode
   [34:41] % branch
   [42:49] excitation energy
   [49:56] Q value
   [63:80] T1/2
   [82:96] Abundance
   [98:105] Atomic mass
   [105:113] mass uncertainity
   [114] S for systematics
   [117:123] data and reference
   [124:132] t1/2 in seconds (WTF?)

    Description of xml document structure - see Nubase2xml.py
"""

import argparse
import xml.dom.minidom
from Nuclide import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reads Nuclear Wallet Cards ascii file and writes xml document')
    parser.add_argument('infile', type=argparse.FileType('r'), 
                         help='Input data file (required)')
    parser.add_argument('outfile', type=argparse.FileType('w'), 
                         help='Output data file (required)')
    args = parser.parse_args()

    line_number = 0
    first = True
    data = []
    previous = {'A' : 0, 'Z': 0, 'mass_defect' : None,
            'half_life': None, 'isomer': None}
    for line in args.infile:
        line_number += 1

        A = int(line[1:4])
        Z = int(line[6:9])
        isomer = True if line[4] == 'M' else False
        spin = line[16:26].strip()
        d_mode = line[30:34].strip()
        d_relation = line[34].strip()
        d_branch = line[35:41].strip()
        isomer_excitation = line[42:49].strip()
        Q_value = line[49:56].strip()
        half_life = line[63:80].strip()
        abundance = line[81:96].strip()
        mass = line[96:105].strip()
        mass_error = line[106:113].strip()
        extrapolated = True if line[114] == 'S' else False
        comment = line[117:123].strip()

        mass_defect = {'value': mass, 'uncertainity': mass_error,
                       'extrapolated' : extrapolated}

        gs_spin = {'value': spin, 'extrapolated': False}

        decay_modes = []
        if len(abundance) > 0:
            stable_data = {}
            abundance = abundance.split()
            stable_data['mode'] = 'is'
            stable_data['value'] = abundance[0].strip('%.')
            stable_data['relation'] = '='
            if len(abundance) > 1:
                stable_data['uncertainity'] = abundance[1]
            else:
                stable_data['uncertainity'] = ''
            decay_modes.append(stable_data)

        if len(d_mode) > 0: 
            decay_mode = {}
            decay_mode['value'] = d_branch
            decay_mode['uncertainity'] = '?'
            d_mode = d_mode.lower()

            if d_mode in ['ep', 'e2p', 'e3p', 'ea', 'e2a']:
                decay_mode['mode'] = 'b+' + d_mode[1:]
            else:
                decay_mode['mode'] = d_mode

            if d_relation == '':
                decay_mode['relation'] = '='
            elif d_relation in ['@', '&']:
                decay_mode['relation'] = '~'
            else:
                decay_mode['relation'] = d_relation

            decay_modes.append(decay_mode)

        if len(decay_modes) == 0:
            empty_mode = {}
            empty_mode['value'] = ''
            empty_mode['uncertainity'] = ''
            empty_mode['mode'] = '?'
            empty_mode['relation'] = ''
            decay_modes.append(empty_mode)

        current = {'A' : A, 'Z': Z, 'mass_defect' : mass_defect,
                   'half_life': half_life, 'isomer': isomer }
        same = True
        for key in current.keys():
            if current[key] != previous[key]:
                same = False
                break

        try:
            if same and not(isomer):
                # Additional decay mode of isotope
                isotope.add_decay_mode(decay_mode)
            elif same and isomer:
                # Additional decay mode of isomer
                isomer_index = len(isotope.isomers) - 1
                isotope.add_isomer_decay_mode(isomer_index, decay_mode)
            elif not(same) and isomer:
                # Add isomer to isotope
                if current['A'] != previous['A'] or current['Z'] != previous['Z']:
                    # There is unresolved dispute which state is ground and
                    # which is isomeric, table gives both as a isomeric
                    # we take first as g.s
                    # To emphasize this fact we also add the same data
                    # to isomer data
                    if not(first):
                        data.append(isotope)
                    else:
                        first = False
                    isotope = NuclideNwc11(Z, A, mass_defect, half_life,
                                           gs_spin, decay_modes, comment)
                    isomer_hl = isotope.nwc_parse_half_life(half_life)
                    isomer_data = {'energy' : isomer_excitation,
                                    'uncertainity' : '?',
                                    'extrapolated' : extrapolated,
                                    'half_life' : isomer_hl,
                                    'decay_modes': decay_modes,
                                    'comment': comment}
                    isotope.add_isomer(isomer_data)
                else:
                    isomer_hl = isotope.nwc_parse_half_life(half_life)
                    isomer_data = {'energy' : isomer_excitation,
                                    'uncertainity' : '?',
                                    'extrapolated' : extrapolated,
                                    'half_life' : isomer_hl,
                                    'decay_modes': decay_modes,
                                    'comment': comment}
                    isotope.add_isomer(isomer_data)
            else:
                # Create new isotope
                if not(first):
                    data.append(isotope)
                else:
                    first = False
                isotope = NuclideNwc11(Z, A, mass_defect, half_life,
                                       gs_spin, decay_modes, comment)

        except ParameterError as err:
            print('Line', line_number, ':', err.msg)
            if isomer:
                print('Skipping bad entry: isomer of isotope {}'.format(isotope))
            else:
                print('Skipping bad entry: isotope {}'.format(isotope))
            print()
        except IndexError:
            print('Index Error, Line', line_number, isotope, isomer)
            exit()
        else:
            previous = current

    dom = xml.dom.minidom.getDOMImplementation()
    table = dom.createDocument(None, "nuclear_data_table", None)
    root = table.documentElement
    for isotope in data:
        isotope.add_to_xml_table(table, root)
    args.outfile.write(table.toprettyxml(indent="    ", encoding="utf-8").decode("utf-8"))

