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
        self.msg = msg
    def __str__(self):
        return repr(self.msg)

class Nuclide(object):
    """Object containing all information about nuclide provided by NuBase ascii file"""

    # Here are units used by NuBase evaluators
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
                        'p-unst': 'p-unstable',
                        'n-unst': 'n-unstable',
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
              'Sb',  'Te',  ' I',  'Xe',  'Cs',\
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

    def __init__(self, Z = 0, A = 0, source = "",  mass_defect = None,
                 half_life = None, gs_spin = None, decay_modes = None,
                 isomers = [], comment = ""):
        """ Constructor. Variable source potentially can be used to parse
        input data in different formats. Currently 'nubase' is supported.

        To load from xml document create dummy nuclide and use
        parse_xml_entry function"""

        super().__init__()
        if source == 'nubase':
            self.Z = Z
            self.A = A
            self.__mass_defect = self.__nb_parse_mass_defect(mass_defect)
            self.__half_life = self.__nb_parse_half_life(half_life)
            self.__gs_spin = self.__nb_parse_gs_spin(gs_spin)
            self.__decay_modes = self.__nb_parse_decay_modes(decay_modes)
        else:
            self.Z = Z
            self.A = A
            if mass_defect != None:
                self.mass_defect = mass_defect
            else:
                self.__mass_defect = None
            if half_life != None:
                self.half_life = half_life
            else:
                self.__half_life = None
            if gs_spin != None:
                self.gs_spin = gs_spin
            else:
                self.__gs_spin = None
            if decay_modes != None:
                self.decay_modes = decay_modes
            else:
                self.__decay_modes = None
        self.comment = comment
        self.isomers = isomers

    def __str__(self):
        return "{}{}".format(self.A, self.element)

    @property
    def element(self):
        """Returns chemical element name"""
        return self._element[self.Z]

    @property
    def Z(self):
        """Returns atomic number Z"""
        return self.__Z

    @Z.setter
    def Z(self, Z):
        """Sets atomic number Z, Z must be integer and Z >= 0"""
        try:
            Z = int(Z)
        except ValueError:
            raise ParameterError("Expecting integer for atomic number Z but '{0}' was found".format(Z))
        if Z < 0:
            raise ParameterError("Z atomic number cannot be smaller then 0, {} was given".format(Z))
        self.__Z = Z

    @property
    def A(self):
        """Returns mass number A"""
        return self.__A

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
        self.__A = A

    @property
    def N(self):
        """Number of neutrons (read only)"""
        return self.__A - self.__Z

    @property
    def mass_defect(self):
        """Returns mass defect data in format:
            [mass, uncertainity, extrapolated]
            extrapolated is boolean
            """
        return self.__mass_defect

    @mass_defect.setter
    def mass_defect(self, mass_defect):
        """Sets mass defect data, format 
        [mass, uncertainity, extrapolated] is expected"""
        if len(mass_defect) != 3:
            raise ParameterError('Mass defect is expected to be a list [mass, error, extrapolated]')
        self.__mass_defect = mass_defect

    @property
    def half_life(self):
        """Half-life is returned as a list [value, unit, uncertainity, extrapolated]"""
        return self.__half_life


    @half_life.setter
    def half_life(self, half_life):
        """Sets half-life data, format 
        [value, unit, uncertainity, extrapolated] is expected"""
        if len(half_life) != 4:
            raise ParameterError('Half life is expected to be a list [value, unit, uncertainity, extrapolated], {} was passed'.format(half_life))
        self.__half_life = half_life

    @property
    def gs_spin(self):
        """Ground state spin"""
        return self.__gs_spin 

    @gs_spin.setter
    def gs_spin(self, gs_spin):
        """Sets g.s. spin, format
        [value, extrapolated] is expected"""
        if len(gs_spin) != 2:
            raise ParameterError('Ground state spin is expected to be a list [value, extrapolated], {} was passed'.format(gs_spin))
        self.__gs_spin = gs_spin

    @property
    def decay_modes(self):
        """A list of decay modes and branching ratios, returns list:
        [mode, relation, value, uncertainity]
        where decay mode is 'B-', 'A', etc, 
        relation is '=', '~', '>=' ,'<=' (greater equal, and less equal coded in unicode)
        and value is given in percent or by '?' if unknown"""
        return self.__decay_modes 


    @decay_modes.setter
    def decay_modes(self, decay_modes):
        """Sets decay_modes data"""
        for mode in decay_modes:
            if len(mode) != 4:
                raise ParameterError('Decay modes is expected to be a list of lists [ [mode, relation, value, error], [mode, ...], ... ], {} was passed'.format(decay_modes))
        self.__decay_modes = decay_modes

    def __nb_parse_mass_defect(self, mass_defect):
        """Returns list [mass, uncertainity, extrapolated]
        parsed from format used by nubase2003"""
        mass_defect = mass_defect.strip()
        extrapolated = True if mass_defect.count('#') > 0 else False
        if extrapolated:
            mass_defect = mass_defect.replace('#', ' ')
        try:
            mass_defect = mass_defect.split()
            for it in mass_defect:
                it = it.strip()
            mass = float(mass_defect[0])
            if len(mass_defect) == 2:
                error = float(mass_defect[1])
            else:
                error = '?'
        except ValueError:
            raise ParameterError(" {} is not valid mass defect string".format(mass_defect))
        return [mass, error, extrapolated]

    def __nb_parse_half_life(self, half_life):
        """Half-life given as a string "value unit" white-space separated
           as in nubase
           However sometimes evaluators use "value unit error"
           or "<value unit"
           or "stbl"
           or empty string

           returns list [half life, unit, uncertainity, extrapolated]
        """
        half_life = half_life.strip()
        extrapolated = True if half_life.count('#') > 0 else False
        if extrapolated:
            half_life = half_life.replace('#', ' ')
        
        items = half_life.split()
        for it in items:
            it = it.strip()

        if len(items) == 0:
            return ['?', '?', '?', extrapolated]
        elif items[0] in ['stbl', 'p-unst', 'n-unst']:
            error = ''
            if len(items) > 1:
                error = items[1]
            return [self._short_time_units[items[0]], '', error, extrapolated]
        elif len(items) == 2 or len(items) == 3:
            if ( self._short_time_units.get(items[1]) == None and
                self._long_time_units.get(items[1]) == None ):
                raise ParameterError(
                      'Could not find half-life unit {}'.format(items[1]) )
            error = '?' if len(items) == 2 else items[2]
            return [items[0], items[1], error, extrapolated]
        else:
            raise ParameterError("String {} is not a valid half life string".format(half_life))

    def __nb_parse_gs_spin(self, gs_spin):
        """Parses nubase style spin information

        returns list [spin, extrapolated] """
        gs_spin = gs_spin.strip()
        extrapolated = True if gs_spin.count('#') > 0 else False
        if extrapolated:
            gs_spin = gs_spin.replace('#', ' ')
        return [gs_spin, extrapolated]

    def __nb_parse_decay_modes(self, decay_modes):
        """Parses decay modes string from nubase

        returns list of lists
        [ [mode, value, relation, uncertainity], [...] ]
        """
        # NuBase evaluators have left some fortran garbage like 'LE' and 'GE'
        # We replace them with proper unicode signs 
        decay_modes = decay_modes.replace('LE', '\u2264')
        decay_modes = decay_modes.replace('GE', '\u2265')
            
        decay_list = [] 
        if len(decay_modes.strip()) == 0:
            decay_list.append(['?', '=', '?', '0'])
            return decay_list
        try:
            for item in decay_modes.split(';'):
                # For unknown values string "mode ?" is used 
                # but sometimes it is "mode=?"
                # or "mode= ?"
                # we fix this so it always has '=' sign
                question = re.search(r" \?", item)
                if question != None:
                    if item.count('=') == 0:
                        item = re.sub(r" \?", "=?", item)

                mode, relation, value = re.split('(=|~|>|<|\u2264|\u2265)',
                                                 item)
                error = '0'
                value = value.split()
                if len(value) > 1:
                    error = value[1]
                value = value[0]
                decay_list.append([mode, relation, value, error])
        except ValueError:
            raise ParameterError('Error parsing decay modes string {} nuclide {}'.format(decay_modes, self))

        return decay_list


    def nb_add_isomer(self, isomer_data, half_life, decay_modes, comment):
        """Adds isomer using nubase style data"""
        half_life = self.__nb_parse_half_life(half_life)
        decay_modes = self.__nb_parse_decay_modes(decay_modes)

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

        self.isomers.append([energy, error, extrapolated, half_life, decay_modes, code + " " + comment])


    def add_isomer(self, isomer):
        """Adds isomer"""
        if len(isomer) != 5 or len(isomer[3]) != 4:
            raise ParameterError('Isomer data is expected to be a list [energy, error, extrapolated, [half_life, unit, uncertainity, extrapolated], [[decay mode, relation, value, uncerainity]]]; {} was passed'.format(isomer))
        self.isomers.append(isomer)


    def parse_xml_entry(self, nuclide):
        """ This function takes entry from xml file and loads all
        data to Nuclide object

        Potentially there's a place for imporvement in more automatic
        xml parsing
        """
        
        self.A = nuclide.getAttribute('A')
        self.Z = nuclide.getAttribute('Z')

        mass_defect = nuclide.getElementsByTagName('mass_defect')[0]
        md_attrs = ['value', 'uncertainity', 'extrapolation']
        md_data = []
        for attr in md_attrs:
            value = mass_defect.getAttribute(attr)
            md_data.append(value)
        self.mass_defect = md_data

        half_life = nuclide.getElementsByTagName('half_life')[0]
        hl_attrs = ['value', 'unit', 'uncertainity', 'extrapolation']
        hl_data = []
        for attr in hl_attrs:
            value = half_life.getAttribute(attr)
            hl_data.append(value)
        self.half_life = hl_data

        spin = nuclide.getElementsByTagName('spin')[0]
        s_attrs = ['value', 'extrapolation']
        s_data = []
        for attr in s_attrs:
            value = spin.getAttribute(attr)
            s_data.append(value)
        self.gs_spin = s_data

        decay_modes = nuclide.getElementsByTagName("decay_modes")[0]
        decay_attr = ("mode", "value",
                        "relation", "uncertainity")
        dm_data = []
        for decay in decay_modes.getElementsByTagName("decay"):
            mode_data = []
            for attr in decay_attr:
                value = decay.getAttribute(attr)
                mode_data.append(value)
            dm_data.append(mode_data)
        self.decay_modes = dm_data

        isomers = nuclide.getElementsByTagName("isomers")
        if len(isomers) > 0:
            for isomer in isomers[0].getElementsByTagName("isomer"):
                i_data = []
                i_attrs = ('energy', 'extrapolation', 'uncertainity')
                for attr in i_attrs:
                    value = isomer.getAttribute(attr)
                    i_data.append(value)

                half_life = isomer.getElementsByTagName('half_life')[0]
                hl_data = []
                hl_attr = ("value", "unit", "uncertainity",
                            "extrapolation")
                for attr in hl_attr:
                    value = half_life.getAttribute(attr)
                    hl_data.append(value)
                i_data.append(hl_data)

                decay_modes = nuclide.getElementsByTagName("decay_modes")[0]
                decay_attr = ("mode", "value",
                            "relation", "uncertainity")
                dm_data = []
                for decay in decay_modes.getElementsByTagName("decay"):
                    mode_data = []
                    for attr in decay_attr:
                        value = decay.getAttribute(attr)
                        mode_data.append(value)
                    dm_data.append(mode_data)
                i_data.append(dm_data)
                self.add_isomer(i_data)
