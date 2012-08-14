"""
   Copyright Krzysztof Miernik 2012
   k.a.miernik@gmail.com 
   
   Distributed under GNU General Public Licence v3

   This module provides Nuclide class used for storing and parsing data from
   NuBase2003 ascii file or xml documents
"""

import xml.dom.minidom
import re

class ParameterError(Exception):
    """Error class for all kinds of wrong parameters passed to 
    all functions and classes in this module"""
    def __init__(self, msg):
        super().__init__()
        self.msg = msg
    def __str__(self):
        return repr(self.msg)

class Nuclide(object):
    """Object containing all information about nuclide. 
       Base class for subclasses for specific format
       (like nubase, nuc. wallet cards, xml etc."""

    # Time units as used by NuBase evaluators
    # It is interesting that half-life can be as short as 1e-24 and at the same
    # time some nuclides are 'proton unstable'
    # in 1 ys the light would travel a distance of 3e-16 m (0.3 fm) 
    # the meaning of 'proton unstable' is very unclear then

    # Units are split into short time (base of 1 second) and
    # long time (base of 1 year) for the purpose of better numerical
    # calculations if neccessary

    # Base unit for short times: second
    # Also added here descriptive 'units'
    _short_time_units = {'d': 86400.0,
                        'h': 3600.0,
                        'm': 60.0,
                        's': 1.0,
                        'ms': 1e-3,
                        'us': 1e-6,
                        'ns': 1e-9,
                        'ps': 1e-12,
                        'fs': 1e-15,
                        'as': 1e-18,
                        'zs': 1e-21,
                        'ys': 1e-24,
                        'stbl' : 'stable',
                        'p-unst': 'unstable',
                        'n-unst': 'unstable',
                        '?': 'unknown'}

    # Base unit for long times: year
    _long_time_units = { 'Yy': 1e24,
                        'Zy': 1e21,
                        'Ey': 1e18,
                        'Py': 1e15,
                        'Ty': 1e12,
                        'Gy': 1e9,
                        'My': 1e6,
                        'ky': 1e3,
                        'y' : 1.0}

    # 1 year as adopted by Nubase 2003
    _long_to_short = 31556926

    # Chemical element name in order of atomic number Z,
    # note that element[0] is neutron 'n'
    _element = 'n' ,\
              'H' ,  'He',  'Li',  'Be',  'B',\
              'C' ,   'N',   'O',   'F',  'Ne',\
              'Na',  'Mg',  'Al',  'Si',   'P',\
              'S' ,  'Cl',  'Ar',   'K',  'Ca',\
              'Sc',  'Ti',   'V',  'Cr',  'Mn',\
              'Fe',  'Co',  'Ni',  'Cu',  'Zn',\
              'Ga',  'Ge',  'As',  'Se',  'Br',\
              'Kr',  'Rb',  'Sr',   'Y',  'Zr',\
              'Nb',  'Mo',  'Tc',  'Ru',  'Rh',\
              'Pd',  'Ag',  'Cd',  'In',  'Sn',\
              'Sb',  'Te',   'I',  'Xe',  'Cs',\
              'Ba',  'La',  'Ce',  'Pr',  'Nd',\
              'Pm',  'Sm',  'Eu',  'Gd',  'Tb',\
              'Dy',  'Ho',  'Er',  'Tm',  'Yb',\
              'Lu',  'Hf',  'Ta',   'W',  'Re',\
              'Os',  'Ir',  'Pt',  'Au',  'Hg',\
              'Tl',  'Pb',  'Bi',  'Po',  'At',\
              'Rn',  'Fr',  'Ra',  'Ac',  'Th',\
              'Pa',   'U',  'Np',  'Pu',  'Am',\
              'Cm',  'Bk',  'Cf',  'Es',  'Fm',\
              'Md',  'No',  'Lr',  'Rf',  'Db',\
              'Sg',  'Bh',  'Hs',  'Mt',  'Ds',\
              'Rg',  'Cn', 'Uut',  'Fl', 'Uup',\
              'Lv', 'Uus', 'Uuo', 'Uue', 'Ubn'

    def __init__(self, Z = 0, A = 0, mass_defect = None,
                 half_life = None, gs_spin = None, decay_modes = None,
                 isomers = None, comment = None):
        """Constructor for base class. All fields are properites and
        sanity check is done therein."""

        super().__init__()
        self.Z = Z
        self.A = A
        self.mass_defect = mass_defect
        self.half_life = half_life
        self.gs_spin = gs_spin
        self.decay_modes = decay_modes
        if isomers is None:
            self.isomers = []
        else:
            self.isomers = isomers
        if comment is None:
            self.comment = ""
        else:
            self.comment = comment


    def __str__(self):
        return "{}{}".format(self.A, self.element)

    @property
    def element(self):
        """Returns chemical element name (read-only)"""
        return self._element[self.Z]

    @property
    def Z(self):
        """Returns atomic number Z"""
        return self._Z

    @Z.setter
    def Z(self, Z):
        """Sets atomic number Z, Z must be integer and Z >= 0"""
        try:
            Z = int(Z)
        except ValueError:
            raise ParameterError("Expecting integer for atomic number Z but '{0}' was found".format(Z))
        if Z < 0:
            raise ParameterError("Z atomic number cannot be smaller then 0, {} was given".format(Z))
        self._Z = Z

    @property
    def A(self):
        """Returns mass number A"""
        return self._A

    @A.setter
    def A(self, A):
        """Sets mass number A. A must be integer and A > 0 with an
        exception of dummy nuclide where Z = 0, A = 0"""
        try:
            A = int(A)
        except ValueError:
            raise ParameterError("Expecting integer for mass number A but '{0}' was found".format(A))
        if A < 1 and self.Z != 0:
            raise ParameterError("A mass number cannot be smaller then 1, {} was given".format(A))
        self._A = A

    @property
    def N(self):
        """Number of neutrons (read only)"""
        return self._A - self._Z

    @property
    def mass_defect(self):
        """Returns mass defect data in format:
           dict {value, uncertainity, extrapolated, ...}
            """
        return self._mass_defect

    @mass_defect.setter
    def mass_defect(self, mass_defect):
        """Sets mass defect data, format 
        dict {mass, uncertainity, extrapolated} is expected"""
        if mass_defect is None:
            self._mass_defect = {}
        else:
            for key in ['value', 'uncertainity', 'extrapolated']:
                if mass_defect.get(key) is None:
                    raise ParameterError("Wrong format of mass defect '{}' was passed".format(mass_defect))
            self._mass_defect = mass_defect

    @property
    def half_life(self):
        """Half-life is returned as a dictionary {value, unit, uncertainity, extrapolated, relation}"""
        return self._half_life

    @half_life.setter
    def half_life(self, half_life):
        """Sets half-life data, format dict
        {value, unit, uncertainity, extrapolated, relation} is expected"""
        if half_life is None:
            self._half_life = {}
        else:
            for key in ['value', 'unit', 'uncertainity', 'relation', 'extrapolated']:
                if half_life.get(key) is None:
                    raise ParameterError("Wrong format of half life '{}' was passed".format(half_life))
            self._half_life = half_life

    @property
    def gs_spin(self):
        """Ground state spin a dict of {value, extrapolated}"""
        return self._gs_spin 

    @gs_spin.setter
    def gs_spin(self, gs_spin):
        """Sets g.s. spin, format
        {'value' : , 'extrapolated' : }
        """
        if gs_spin is None:
            self._gs_spin = {}
        else:
            for key in ['value', 'extrapolated']:
                if gs_spin.get(key) is None:
                    raise ParameterError("Wrong format of ground state spin '{}' was passed".format(gs_spin))
            self._gs_spin = gs_spin

    @property
    def decay_modes(self):
        """A list of decay modes and branching ratios, returns list of dictionaries
        [ {'mode': , 'relation': , 'value' : , 'uncertainity': }, {}, ...]
        where decay mode is 'B-', 'A', etc, 
        relation is '=', '~', '>=' ,'<=' (greater equal, and less equal coded in unicode)
        and value is given in percent or by '?' if unknown"""
        return self._decay_modes 

    @decay_modes.setter
    def decay_modes(self, decay_modes):
        """Sets decay_modes data"""
        self._decay_modes = []
        if decay_modes is not None:
            for mode in decay_modes:
                self.add_decay_mode(mode)

    def add_decay_mode(self, decay_mode):
        """Decay mode should be a dictionary:
        {'mode': , 'relation': , 'value' : , 'uncertainity': }
        """
        for key in ['mode', 'relation', 'value', 'uncertainity']:
            if decay_mode.get(key) is None:
                raise ParameterError("Wrong format of decay mode '{}' was passed".format(decay_mode))
        self._decay_modes.append(decay_mode)

    def add_isomer(self, isomer):
        """Adds isomer
        isomer data should be a dictionary:
        {energy, uncertainity, extrapolated, 
         half_life : {half_life, unit, uncertainity, extrapolated},
         decay_modes: [{mode, relation, value, uncertainity},{},...] }
        """
        for key in ['energy', 'uncertainity', 'extrapolated', 'half_life', 'decay_modes']:
            if isomer.get(key) is None:
                raise ParameterError("Wrong format of isomer data '{}' was passed; key {} error".format(isomer, key))
        for key in ['value', 'unit', 'uncertainity', 'relation', 'extrapolated']:
            if isomer['half_life'].get(key) is None:
                raise ParameterError("Wrong format of isomer data '{}' was passed; key {} error".format(isomer, key))
        for mode in isomer['decay_modes']:
            for key in ['mode', 'relation', 'uncertainity', 'value']:
                if mode.get(key) is None:
                    raise ParameterError("Wrong format of isomer data '{}' was passed; key {} error".format(isomer))
        self.isomers.append(isomer)


    def add_isomer_decay_mode(self, isomer_index, decay_mode):
        """
        Adds decay mode to isomer data selected by isomer_index
        """
        for key in ['mode', 'relation', 'uncertainity', 'value']:
            if decay_mode.get(key) is None:
                raise ParameterError("Wrong format of isomer decay data '{}' was passed".format(isomer))
        self.isomers[isomer_index]['decay_modes'].append(decay_mode)

    def add_to_xml_table(self, table, root):
        """Adds nuclide to xml document;
        table is dom Document
        root is Document root element
        e.g
            table = dom.createDocument(None, "nuclear_data_table", None)
            root = table.documentElement
        """
        nuclide = table.createElement("nuclide")
        nuclide.setAttribute("Z", str(self.Z))
        nuclide.setAttribute("A", str(self.A))
        nuclide.setAttribute("id", str(self))
        nuclide.setAttribute("element", self.element)

        mass_defect = table.createElement("mass_defect")
        for key, value in self.mass_defect.items():
            mass_defect.setAttribute(key, str(value))
        nuclide.appendChild(mass_defect)

        half_life = table.createElement("half_life")
        for key, value in self.half_life.items():
            half_life.setAttribute(key, str(value))
        nuclide.appendChild(half_life)

        spin = table.createElement("spin")
        for key, value in self.gs_spin.items():
            spin.setAttribute(key, str(value))
        nuclide.appendChild(spin)

        decay_modes = table.createElement("decay_modes")
        for mode in self.decay_modes:
            decay = table.createElement("decay")
            for key, value in mode.items():
                decay.setAttribute(key, str(value))
            decay_modes.appendChild(decay)
        nuclide.appendChild(decay_modes)

        isomers = table.createElement("isomers")
        for state in self.isomers:
            isomer = table.createElement("isomer")
            for key in ['energy', 'uncertainity', 'extrapolated']:
                isomer.setAttribute(key, str(state[key]))

            ihalf_life = table.createElement("half_life")
            for key, value in state['half_life'].items():
                ihalf_life.setAttribute(key, str(value))
            isomer.appendChild(ihalf_life)

            idecay_modes = table.createElement("decay_modes")
            for mode in state['decay_modes']:
                idecay = table.createElement("decay")
                for key, value in mode.items():
                    idecay.setAttribute(key, str(value))
                idecay_modes.appendChild(idecay)
            isomer.appendChild(idecay_modes)

            comment_text = table.createTextNode(state['comment'])
            comment = table.createElement("comment")
            comment.appendChild(comment_text)
            isomer.appendChild(comment)

            isomers.appendChild(isomer)
        
        if len(self.isomers) > 0:
            nuclide.appendChild(isomers)

        comment_text = table.createTextNode(self.comment)
        comment = table.createElement("comment")
        comment.appendChild(comment_text)
        nuclide.appendChild(comment)
        root.appendChild(nuclide)


class NuclideNb03(Nuclide):
    """Nuclide class for reading data from Nubase2003 format"""

    def __init__(self, Z, A, mass_defect, half_life,
                 gs_spin, decay_modes, comment = ""):
        """Constructor for Nubase03 version.
        Comment:
        A crucial for efficiency is to call
        super().__init__(Z, A, {}, {}, {}, [], [], "")
        putting all lists explicitely

        Otherwise objects will be linked by a common references to 
        memeber (if member is not changed e.g isomers)
        and they will not be scheduled for garbage collection!

        Now more pythonic approach
        __init__(something = None):
            if something is None:
                self.something = []
            else:
                self.something = something
        """
        super().__init__(Z, A)
        self.mass_defect = self._parse_mass_defect(mass_defect)
        self.half_life = self._parse_half_life(half_life)
        self.gs_spin = self._parse_gs_spin(gs_spin)
        self.decay_modes = self._parse_decay_modes(decay_modes)
        self.comment = comment

    def _parse_mass_defect(self, mass_defect):
        """Returns dict {value, uncertainity, extrapolated}
        parsed from format used by nubase2003"""
        result = {}
        mass_defect = mass_defect.strip()
        result['extrapolated'] = True if mass_defect.count('#') > 0 else False
        if result['extrapolated']:
            mass_defect = mass_defect.replace('#', ' ')
        try:
            mass_defect = mass_defect.split()
            for it in mass_defect:
                it = it.strip()
            result['value'] = float(mass_defect[0])
            if len(mass_defect) == 2:
                result['uncertainity'] = float(mass_defect[1])
            else:
                result['uncertainity'] = '?'
        except ValueError:
            raise ParameterError(" {} is not valid mass defect string".format(mass_defect))
        return result

    def _parse_half_life(self, half_life):
        """Half-life given as a string "value unit" white-space separated
           as in nubase
           However sometimes evaluators use "value unit error"
           or "<value unit"
           or "stbl"
           or empty string

           returns dict {half life, unit, uncertainity, extrapolated}
        """
        result = {}
        half_life = half_life.strip()
        result['extrapolated'] = True if half_life.count('#') > 0 else False
        if result['extrapolated']:
            half_life = half_life.replace('#', ' ')
        
        items = half_life.split()
        for it in items:
            it = it.strip()

        if len(items) == 0:
            result['value'] = '?'
            result['uncertainity'] = '?'
            result['unit'] = '?'
            result['relation'] = '?'
        elif items[0] in ['stbl', 'p-unst', 'n-unst']:
            result['value'] = self._short_time_units[items[0]]
            result['uncertainity'] = ''
            result['unit'] = ''
            result['relation'] = '='
        elif len(items) == 2 or len(items) == 3:
            if ( self._short_time_units.get(items[1]) is None and
                self._long_time_units.get(items[1]) is None ):
                raise ParameterError(
                      'Could not find half-life unit {}'.format(items[1]) )
            result['uncertainity'] = '?' if len(items) == 2 else items[2]
            result['value'] = items[0]
            result['unit'] = items[1] 
            result['relation'] = '='
        else:
            raise ParameterError("String {} is not a valid half life string".format(half_life))
        return result

    def _parse_gs_spin(self, gs_spin):
        """Parses nubase style spin information

        returns dictionary {value, extrapolated} """
        result = {}
        gs_spin = gs_spin.strip()
        result['extrapolated'] = True if gs_spin.count('#') > 0 else False
        if result['extrapolated']:
            gs_spin = gs_spin.replace('#', ' ')
        result['value'] = gs_spin
        return result

    def _parse_decay_modes(self, decay_modes):
        """Parses decay modes string from nubase

        returns list of dictionaries
        [ {mode, value, relation, uncertainity}, {...} ]
        """
        # NuBase evaluators have left some fortran garbage like 'LE' and 'GE'
        # We replace them with proper unicode signs 
        decay_modes = decay_modes.replace('le', '\u2264')
        decay_modes = decay_modes.replace('ge', '\u2265')
            
        decay_list = [] 
        if len(decay_modes.strip()) == 0:
            empty = {'mode' : '?', 'value' : '', 'relation' : '', 
                     'uncertainity' : ''}
            decay_list.append(empty)
            return decay_list
        try:
            for item in decay_modes.split(';'):
                # For unknown values string "mode ?" is used 
                # but sometimes it is "mode=?"
                # or "mode= ?"
                # we fix this so it always has '=' sign
                question = re.search(r" \?", item)
                if question is not None:
                    if item.count('=') == 0:
                        item = re.sub(r" \?", "=?", item)

                mode, relation, value = re.split('(=|~|>|<|\u2264|\u2265)',
                                                 item)
                error = '0'
                value = value.split()
                if len(value) > 1:
                    error = value[1]
                value = value[0]
                decay_mode = {'mode': mode, 'relation' : relation, 'value' : value,
                               'uncertainity' : error}
                decay_list.append(decay_mode)
        except ValueError:
            raise ParameterError('Error parsing decay modes string {} nuclide {}'.format(decay_modes, self))
        return decay_list


    def nb_add_isomer(self, isomer_data, half_life, decay_modes, comment):
        """Adds isomer using nubase style data"""
        half_life = self._parse_half_life(half_life)
        decay_modes = self._parse_decay_modes(decay_modes)

        isomer_data = isomer_data.strip()

        extrapolated = True if isomer_data.count('#') > 0 else False
        if extrapolated:
            isomer_data = isomer_data.replace('#', ' ')

        isomer_data = isomer_data.split()
        try:
            energy = isomer_data[0].strip()
            error = isomer_data[1].strip()
            # Default is gamma spectrometry
            code = 'Gamma spectometry'
            if len(isomer_data) >= 3 and not(isomer_data[-1].isnumeric()):
                code = isomer_data[-1].strip()
                starred = True if code.count('*') > 0 else False
                if starred:
                    code_comment = 'Uncertainity of energy is larger then energy itself'
                    code = code.replace('*', ' ')
                if code == 'MD':
                    code = 'Mass doublet'
                elif code == 'RQ':
                    code = 'Reaction energy difference'
                elif code == 'AD':
                    code = 'Alpha energy difference'
                elif code == 'BD':
                    code = 'Beta energy difference'
                elif code == 'p':
                    code = 'Proton decay'
                elif code == 'XL':
                    code = 'L X-rays'
                elif code == 'Nm':
                    code = 'Estimated value from Nilsson model'
                elif code == 'EU':
                    code = 'Existence under discussion'
                elif code == 'RN':
                    code = 'Proved not to exists'
                elif code == '&':
                    code = 'Ground state and isomer ordering reversed compared to ENSDF'
                else:
                    code = "Code '{}' is not documented'".format(code)
                    
                if starred:
                    code += " " + code_comment

        except (IndexError, ValueError):
            raise ParameterError('Error parsing isomer data string {} nuclide {}'.format(isomer_data, self))

        result = { 'energy' : energy,
                   'uncertainity': error,
                   'extrapolated': extrapolated,
                   'half_life': half_life,
                   'decay_modes' : decay_modes,
                   'comment': code + " " + comment }
        self.add_isomer(result)


class NuclideNwc11(Nuclide):
    """Nuclide class for reading data from Nuclear Wallet Cards 2011 format"""

    def __init__(self, Z, A, mass_defect, half_life, gs_spin,
                 decay_modes, comment = None):
        """ Constructor for Nuclear Wallet Cards 2011 version."""
        super().__init__(Z, A, mass_defect, half_life=None, gs_spin=gs_spin,
                         decay_modes=decay_modes, isomers=None, comment=comment)
        self.half_life = self.nwc_parse_half_life(half_life)

    def nwc_parse_half_life(self, half_life):
        """Half-life given as a string "value unit uncertainty" white-space separated
           as in nuclear waller cards
           evaluators use units in capital letters
           also units EV, KEV, MEV
           No extrapolated values (or no documentation on that?)

           returns dictionary {half life, unit, uncertainity, extrapolated}
        """
        result = {}
        half_life = half_life.lower().strip()

        # Placeholder. If documentation on extrapolation founded it 
        # might be change to some function
        result['extrapolated'] = False

        items = half_life.split()
        for it in items:
            it = it.strip()
        if len(items) == 2:
            items.append('?')

        if len(items) == 0:
            result['value'] = '?'
            result['uncertainity'] = '?'
            result['unit'] = '?'
            result['relation'] = '?'
        elif items[0] in ['stable', 'unbound']:
            result['value'] = items[0] if items[0] == 'stable' else 'unstable'
            result['uncertainity'] = ''
            result['unit'] = ''
            result['relation'] = '='
            if len(items) > 1:
                result['uncertainity'] = items[1]
        elif len(items) == 3:
            if items[1] in ['ev', 'kev', 'mev']:
                if items[1] == 'ev':
                    items[1] = 'as'
                elif items[1] == 'kev':
                    items[1] = 'zs'
                elif items[1] == 'mev':
                    items[1] = 'ys'

                try:
                    items[0] = '{0:.5f}'.format(items[0] * 0.04562)
                except TypeError:
                    pass
                try:
                    items[2] = '{0:.5f}'.format(items[2] * 0.04562)
                except TypeError:
                    pass

            if ( self._short_time_units.get(items[1]) is None and
                 self._long_time_units.get(items[1]) is None ):
                raise ParameterError(
                      'Could not find half-life unit {}'.format(items[1]) )
            result['value'] = items[0]
            result['unit'] = items[1]
            result['uncertainity'] = '?'
            if items[2] == 'ap':
                result['relation'] = '~' 
            elif items[2] == 'lt':
                result['relation'] = '<' 
            elif items[2] == 'le':
                result['relation'] = '\u2264' 
            elif items[2] == 'gt':
                result['relation'] = '>' 
            elif items[2] == 'ge':
                result['relation'] = '\u2265' 
            else:
                result['relation'] = '='
                result['uncertainity'] = items[2]
        else:
            raise ParameterError("String {} is not a valid half life string".format(half_life))
        return result

class NuclideXml(Nuclide):
    """Nuclide class for reading and writing data to xml document"""

    def __init__(self, Z, A, xml_nuclide_entry = None):
        """ Constructor
        xml_nuclide_entry should be an dom element with tag name 'nuclide'

        e.g (simplified):
        dom = xml.dom.minidom.parse(datafile)
        for xml_nuclide_entry in dom.getElementsByTagName("nuclide"):
            isotope = NuclideXml(A, Z, xml_nuclide_entry)
        """
        super().__init__(Z, A)
        if xml_nuclide_entry is not None:
            self.parse_xml_entry(xml_nuclide_entry)

    def parse_xml_entry(self, xml_nuclide_entry):
        """ This function takes entry from xml file and loads all
        data to Nuclide object (xml entry by name of nuclide)

        Potentially there's a place for imporvement in more automatic
        xml parsing
        """
        nuclide = xml_nuclide_entry
        
        self.A = nuclide.getAttribute('A')
        self.Z = nuclide.getAttribute('Z')

        mass_defect = nuclide.getElementsByTagName('mass_defect')[0]
        md_attrs = ['value', 'uncertainity', 'extrapolated']
        md_data = {}
        for attr in md_attrs:
            value = mass_defect.getAttribute(attr)
            md_data[attr] = value
        self.mass_defect = md_data

        half_life = nuclide.getElementsByTagName('half_life')[0]
        hl_attr = ['value', 'unit', 'uncertainity', 'relation', 'extrapolated']
        hl_data = {}
        for attr in hl_attr:
            value = half_life.getAttribute(attr)
            hl_data[attr] = value
        self.half_life = hl_data

        spin = nuclide.getElementsByTagName('spin')[0]
        s_attrs = ['value', 'extrapolated']
        s_data = {}
        for attr in s_attrs:
            value = spin.getAttribute(attr)
            s_data[attr] = value
        self.gs_spin = s_data

        decay_modes = nuclide.getElementsByTagName("decay_modes")[0]
        decay_attr = ["mode", "value", "relation", "uncertainity"]
        dm_data = []
        for decay in decay_modes.getElementsByTagName("decay"):
            mode_data = {}
            for attr in decay_attr:
                value = decay.getAttribute(attr)
                mode_data[attr] = value
            dm_data.append(mode_data)
        self.decay_modes = dm_data

        isomers = nuclide.getElementsByTagName("isomers")
        if len(isomers) > 0:
            for isomer in isomers[0].getElementsByTagName("isomer"):
                i_data = {}
                i_attrs = ['energy', 'extrapolated', 'uncertainity']
                for attr in i_attrs:
                    value = isomer.getAttribute(attr)
                    i_data[attr] = value

                half_life = isomer.getElementsByTagName('half_life')[0]
                hl_data = {}
                for attr in hl_attr:
                    value = half_life.getAttribute(attr)
                    hl_data[attr] = value
                i_data['half_life'] = hl_data

                decay_modes = nuclide.getElementsByTagName("decay_modes")[0]
                dm_data = []
                for decay in decay_modes.getElementsByTagName("decay"):
                    mode_data = {}
                    for attr in decay_attr:
                        value = decay.getAttribute(attr)
                        mode_data[attr] = value
                    dm_data.append(mode_data)
                i_data['decay_modes'] = dm_data
                self.add_isomer(i_data)
