"""
Example of searching xml data file for specific characteristics of nuclides
"""

import xml.dom.minidom
from Nuclide import *


def load_xml_nuclear_table(datafile):
    try:
        dom = xml.dom.minidom.parse(datafile)
    except (EnvironmentError, xml.parsers.expat.ExpatError) as err:
        print("{0}: import error: {1}".format(datafile, err))
        return None

    data = []
    for nuclide in dom.getElementsByTagName("nuclide"):
        try:
            A = int(nuclide.getAttribute('A'))
            Z = int(nuclide.getAttribute('Z'))
            N = A - Z

            isotope = NuclideXml(Z, A, nuclide)
            data.append(isotope)
        except (ValueError, LookupError) as err:
            print("{0}: import error: {1}".format(datafile, err))
            return False
    return data


data = load_xml_nuclear_table("nubase12.xml")

for nuclide in data:
    for mode in nuclide.decay_modes:
        try:
            value = float(mode['value'])
        except ValueError:
            value = 0.0

        if (mode['mode'] == 'sf' and value >= 10.0 and
                mode['relation'] in ['=', '~', '>']):
            print('{}  {: >6} {: >3} {: >5.1f}'.format(
                    nuclide,
                    nuclide.half_life['value'],
                    nuclide.half_life['unit'], 
                    value)
                    )

    for isomer in nuclide.isomers:
        if isomer['half_life']['unit'] not in ['s', 'ms', 'us']:
            continue
        for mode in isomer['decay_modes']:
            try:
                value = float(mode['value'])
            except ValueError:
                value = 0.0

            if (mode['mode'] == 'sf' and value >= 10.0 and
                    mode['relation'] in ['=', '~', '>']):
                print('{}* {: >6} {: >3} {: >5.1f}'.format(
                       nuclide,
                       isomer['half_life']['value'],
                      isomer['half_life']['unit'], 
                      value)
                      )

