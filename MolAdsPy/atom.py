from ._exception import BasicException
from numpy import array,ndarray


class AtomError(BasicException):
    pass


class Atom:
    _curratomid = None  # ID available for a new Atom object
    _instances = None  # List of Atom objects created so far
    __slots__ = ["_symbol", "_x", "_fixed", "_active", "_belongs_to", "_id"]

    def __init__(
        self,
        symbol,
        x=array([0.0, 0.0, 0.0]),
        fixed=array([False, False, False]),
        active=True,
    ):
        """
        Object initialization.

        Parameters
        ----------
        symbol : string
            Symbol of the chemical element.
        x : Numpy array, optional
            Cartesian coordinates of the atom. It can also be providade as a
            Python list or tuple. The default is [0.0 0.0 0.0].
        fixed : Numpy array, optional
            Defines if a component of the atomic coordinates remains unaltered
            in the course of a geometry optimization. It can also be provided as
            a Python list or tuple. The default is [False False False].
        active : logical, optional
            Whether the atom is active or not. The default is True.
            
        """
        if isinstance(symbol, str) and len(symbol) > 0 and len(symbol) <= 4:
            self._symbol = symbol

            if not symbol in Species:
                print("WARNING: This atomic species is not in the Species dictionary!")
        else:
            raise AtomError(
                "'symbol' must be provided as a non-empty string with up to 4 characters!"
            )

        if isinstance(x, (list, tuple)):
            x = array(x)

        if isinstance(x, ndarray) and x.shape[0] == 3:
            self._x = x.astype(float)
        else:
            raise AtomError("'x' must be an array with three components!")

        if isinstance(fixed, (list, tuple)):
            fixed = array(fixed)

        if isinstance(fixed, ndarray) and fixed.shape[0] == 3:
            self._fixed = fixed.astype(bool)
        else:
            raise AtomError("'fixed' must be an array with three components!")

        if isinstance(active, bool):
            self._active = active
        else:
            raise AtomError("'active' must be a boolean!")

        self._belongs_to = None

        type(self)._append(self)

    def displace(self, disp):
        """
        Displaces the atom.

        Parameters
        ----------
        disp : Numpy array
             Displacement vector. It can also be provided as a Python list or
             tuple.
             
        """
        if isinstance(disp, (list, tuple)):
            disp = array(disp)

        if isinstance(disp, ndarray) and disp.shape[0] == 3:
            self._x += disp.astype(float)

            if self._belongs_to is not None:
                self._belongs_to._update()
        else:
            raise AtomError("'disp' must be an array with three components!")

    def copy(self):
        """
        Returns a copy of the current Atom object.

        Returns
        -------
        Atom object
            Copy of the atom.
            
        """
        return Atom(self._symbol, self._x, self._fixed, self._active)

    def __str__(self):
        """
        Returns the symbol, the element and the ID of the Atom object.

        Returns
        -------
        String
            A string containing the symbol, the element and the ID of the Atom object.
            
        """
        return "<Atom object> Symbol: %s; Element: %s; ID: %d" % (
            self._symbol,
            self.element,
            self._id,
        )

    @classmethod
    def _append(cls, atom):
        """
        Appends an Atom object to the list of Atom objects created so far.

        Parameters
        ----------
        atom : Atom object
            Atom to be added to the list of Atom objects created so far.

        """
        if cls._curratomid is None:
            cls._curratomid = 0
            cls._instances = []

        atom._id = cls._curratomid
        cls._instances.append(atom)
        cls._curratomid += 1

    """
    Properties
    ----------
    symbol : string, readonly
        Symbol of the chemical element.
    element : string, readonly
        Name of the chemical element.
    atomic_number : integer, readonly
        Atomic number of the element.
    atomic_mass : float, readonly
        Atomic mass of the element, in atomic units.
    atomic_radius : float, readonly
        Atomic radius of the element, in Angstroms.
    valence_electrons : integer, readonly
        Number of electrons in the element's outermost shell.
    valency : Python tuple, readonly
        Number of electrons that the atom can lose or gain when forming bonds.
    electronegativity : float, readonly
        Tendency of the element to attract electrons, in Pauling scale.
    oxidation_states : Python tuple, readonly
        Theoretical charges of the atom if its bonds were fully ionic.    
    coords : Numpy array
        Cartesian coordinates of the atom. It can also be provided as a Python 
        list or tuple.
    isfixed : Numpy array
        Defines if an atom component should remain unaltered in the course
        of a geometry optimization. It can also be provided as a Python list or 
        tuple.
    active : logical
        Determine whether the atom is active within the structure it is part of.
    belongs_to : AtomCollection-derived object, readonly
        Structure to which the AtomCollection object belongs.
    ID : integer, readonly
        Unique atom identifier.
        
    """

    @property
    def symbol(self):
        return self._symbol

    @property
    def element(self):
        if self._symbol in Species and "element" in Species[self._symbol]:
            return Species[self._symbol]["element"]
        else:
            return None

    @property
    def atomic_number(self):
        if self._symbol in Species and "atomic number" in Species[self._symbol]:
            return Species[self._symbol]["atomic number"]
        else:
            return None

    @property
    def atomic_mass(self):
        if self._symbol in Species and "atomic mass" in Species[self._symbol]:
            return Species[self._symbol]["atomic mass"]
        else:
            return None

    @property
    def atomic_radius(self):
        if self._symbol in Species and "atomic radius" in Species[self._symbol]:
            return Species[self._symbol]["atomic radius"]
        else:
            return None

    @property
    def valence_electrons(self):
        if self._symbol in Species and "valence electrons" in Species[self._symbol]:
            return Species[self._symbol]["valence electrons"]
        else:
            return None

    @property
    def valency(self):
        if self._symbol in Species and "valency" in Species[self._symbol]:
            return Species[self._symbol]["valency"]
        else:
            return None

    @property
    def electronegativity(self):
        if self._symbol in Species and "electronegativity" in Species[self._symbol]:
            return Species[self._symbol]["electronegativity"]
        else:
            return None

    @property
    def oxidation_states(self):
        if self._symbol in Species and "oxidation states" in Species[self._symbol]:
            return Species[self._symbol]["oxidation states"]
        else:
            return None

    @property
    def coords(self):
        return self._x

    @coords.setter
    def coords(self, val):
        if isinstance(val, (list, tuple)):
            val = array(val)

        if isinstance(val, ndarray) and val.shape[0] == 3:
            self._x = val.astype(float)

            if self._belongs_to is not None:
                self._belongs_to._update()
        else:
            raise AtomError("Atom coordinates must be an array with three components!")

    @property
    def isfixed(self):
        return self._fixed

    @isfixed.setter
    def isfixed(self, val):
        if isinstance(val, (list, tuple)):
            val = array(val)

        if isinstance(val, ndarray) and val.shape[0] == 3:
            self._fixed = val.astype(bool)
        else:
            raise AtomError("'isfixed' must be an array with three components!")

    @property
    def active(self):
        return self._active

    @active.setter
    def active(self, val):
        if isinstance(val, bool):
            self._active = val

            if self._belongs_to is not None:
                self._belongs_to._update()
        else:
            raise AtomError("'active' must be boolean!")

    @property
    def belongs_to(self):
        return self._belongs_to

    @property
    def ID(self):
        return self._id


Species = {
    "H": {
        "element": "hydrogen",
        "atomic number": 1,
        "atomic mass": 1.0078,
        "atomic radius": 0.53,
        "valence electrons": 1,
        "valency": (1,),
        "electronegativity": 2.2,
        "oxidation states": (-1, 1),
    },
    "He": {
        "element": "helium",
        "atomic number": 2,
        "atomic mass": 4.0026,
        "atomic radius": 0.31,
        "valence electrons": 2,
        "valency": (0,),
        "electronegativity": 0.0,
        "oxidation states": (0,),
    },
    "Li": {
        "element": "lithium",
        "atomic number": 3,
        "atomic mass": 6.941,
        "atomic radius": 1.67,
        "valence electrons": 1,
        "valency": (1,),
        "electronegativity": 0.98,
        "oxidation states": (1,),
    },
    "Be": {
        "element": "berylium",
        "atomic number": 4,
        "atomic mass": 9.0122,
        "atomic radius": 1.12,
        "valence electrons": 2,
        "valency": (2,),
        "electronegativity": 1.57,
        "oxidation states": (0, 1, 2),
    },
    "B": {
        "element": "boron",
        "atomic number": 5,
        "atomic mass": 10.811,
        "atomic radius": 0.87,
        "valence electrons": 3,
        "valency": (3,),
        "electronegativity": 2.04,
        "oxidation states": (-5, -1, 0, 1, 2, 3),
    },
    "C": {
        "element": "carbon",
        "atomic number": 6,
        "atomic mass": 12.011,
        "atomic radius": 0.67,
        "valence electrons": 4,
        "valency": (4,),
        "electronegativity": 2.55,
        "oxidation states": (-4, -3, -2, -1, 0, 1, 2, 3, 4),
    },
    "N": {
        "element": "nitrogen",
        "atomic number": 7,
        "atomic mass": 14.007,
        "atomic radius": 0.56,
        "valence electrons": 5,
        "valency": (3,),
        "electronegativity": 3.04,
        "oxidation states": (-3, -2, -1, 1, 2, 3, 4, 5),
    },
    "O": {
        "element": "oxygen",
        "atomic number": 8,
        "atomic mass": 15.999,
        "atomic radius": 0.48,
        "valence electrons": 6,
        "valency": (2,),
        "electronegativity": 3.44,
        "oxidation states": (-2, -1, 0, 1, 2),
    },
    "F": {
        "element": "fluorine",
        "atomic number": 9,
        "atomic mass": 18.998,
        "atomic radius": 0.42,
        "valence electrons": 7,
        "valency": (1,),
        "electronegativity": 3.98,
        "oxidation states": (-1, 0),
    },
    "Ne": {
        "element": "neon",
        "atomic number": 10,
        "atomic mass": 20.18,
        "atomic radius": 0.38,
        "valence electrons": 8,
        "valency": (0,),
        "electronegativity": 0.0,
        "oxidation states": (0,),
    },
    "Na": {
        "element": "sodium",
        "atomic number": 11,
        "atomic mass": 22.99,
        "atomic radius": 1.90,
        "valence electrons": 1,
        "valency": (1,),
        "electronegativity": 0.93,
        "oxidation states": (-1, 1),
    },
    "Mg": {
        "element": "magnesium",
        "atomic number": 12,
        "atomic mass": 24.305,
        "atomic radius": 1.45,
        "valence electrons": 2,
        "valency": (2,),
        "electronegativity": 1.31,
        "oxidation states": (1, 2),
    },
    "Al": {
        "element": "aluminium",
        "atomic number": 13,
        "atomic mass": 26.982,
        "atomic radius": 1.18,
        "valence electrons": 3,
        "valency": (3,),
        "electronegativity": 1.61,
        "oxidation states": (-2, -1, 1, 2, 3),
    },
    "Si": {
        "element": "silicon",
        "atomic number": 14,
        "atomic mass": 28.086,
        "atomic radius": 1.11,
        "valence electrons": 4,
        "valency": (4,),
        "electronegativity": 1.9,
        "oxidation states": (-4, -3, -2, -1, 0, 1, 2, 3, 4),
    },
    "P": {
        "element": "phosphorus",
        "atomic number": 15,
        "atomic mass": 30.974,
        "atomic radius": 0.98,
        "valence electrons": 5,
        "valency": (3, 5),
        "electronegativity": 2.19,
        "oxidation states": (-3, -2, -1, 0, 1, 2, 3, 4, 5),
    },
    "S": {
        "element": "sulfur",
        "atomic number": 16,
        "atomic mass": 32.065,
        "atomic radius": 0.87,
        "valence electrons": 6,
        "valency": (2,),
        "electronegativity": 2.58,
        "oxidation states": (-2, -1, 0, 1, 2, 3, 4, 5, 6),
    },
    "Cl": {
        "element": "chlorine",
        "atomic number": 17,
        "atomic mass": 35.453,
        "atomic radius": 0.79,
        "valence electrons": 7,
        "valency": (1,),
        "electronegativity": 3.16,
        "oxidation states": (-1, 1, 2, 3, 4, 5, 6, 7),
    },
    "Ar": {
        "element": "argon",
        "atomic number": 18,
        "atomic mass": 39.948,
        "atomic radius": 0.71,
        "valence electrons": 8,
        "valency": (0,),
        "electronegativity": 0.0,
        "oxidation states": (0,),
    },
    "K": {
        "element": "potassium",
        "atomic number": 19,
        "atomic mass": 39.098,
        "atomic radius": 2.43,
        "valence electrons": 1,
        "valency": (1,),
        "electronegativity": 0.82,
        "oxidation states": (-1, 1),
    },
    "Ca": {
        "element": "calcium",
        "atomic number": 20,
        "atomic mass": 40.078,
        "atomic radius": 1.94,
        "valence electrons": 2,
        "valency": (2,),
        "electronegativity": 1.0,
        "oxidation states": (1, 2),
    },
    "Sc": {
        "element": "scandium",
        "atomic number": 21,
        "atomic mass": 44.956,
        "atomic radius": 1.84,
        "valence electrons": 3,
        "valency": (3,),
        "electronegativity": 1.36,
        "oxidation states": (0, 1, 2, 3),
    },
    "Ti": {
        "element": "titanium",
        "atomic number": 22,
        "atomic mass": 47.867,
        "atomic radius": 1.76,
        "valence electrons": 4,
        "valency": (2, 3, 4),
        "electronegativity": 1.54,
        "oxidation states": (-2, -1, 0, 1, 2, 3, 4),
    },
    "V": {
        "element": "vanadium",
        "atomic number": 23,
        "atomic mass": 50.942,
        "atomic radius": 1.71,
        "valence electrons": 5,
        "valency": (2, 3, 4, 5),
        "electronegativity": 1.63,
        "oxidation states": (-3, -1, 0, 1, 2, 3, 4, 5),
    },
    "Cr": {
        "element": "chromium",
        "atomic number": 24,
        "atomic mass": 51.996,
        "atomic radius": 1.66,
        "valence electrons": 6,
        "valency": (2,),
        "electronegativity": 1.66,
        "oxidation states": (-4, -2, -1, 0, 1, 2, 3, 4, 5, 6),
    },
    "Mn": {
        "element": "manganese",
        "atomic number": 25,
        "atomic mass": 54.938,
        "atomic radius": 1.61,
        "valence electrons": 7,
        "valency": (2, 4, 7),
        "electronegativity": 1.55,
        "oxidation states": (-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7),
    },
    "Fe": {
        "element": "iron",
        "atomic number": 26,
        "atomic mass": 55.845,
        "atomic radius": 1.56,
        "valence electrons": 8,
        "valency": (2, 3),
        "electronegativity": 1.83,
        "oxidation states": (-4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7),
    },
    "Co": {
        "element": "cobalt",
        "atomic number": 27,
        "atomic mass": 58.933,
        "atomic radius": 1.52,
        "valence electrons": 9,
        "valency": (2, 3),
        "electronegativity": 1.88,
        "oxidation states": (-3, -1, 0, 1, 2, 3, 4, 5),
    },
    "Ni": {
        "element": "nickel",
        "atomic number": 28,
        "atomic mass": 58.693,
        "atomic radius": 1.49,
        "valence electrons": 10,
        "valency": (3, 2),
        "electronegativity": 1.91,
        "oxidation states": (-2, -1, 0, 1, 2, 3, 4),
    },
    "Cu": {
        "element": "copper",
        "atomic number": 29,
        "atomic mass": 63.546,
        "atomic radius": 1.45,
        "valence electrons": 1,
        "valency": (1, 2),
        "electronegativity": 1.9,
        "oxidation states": (-2, 0, 1, 2, 3, 4),
    },
    "Zn": {
        "element": "zinc",
        "atomic number": 30,
        "atomic mass": 65.380,
        "atomic radius": 1.42,
        "valence electrons": 2,
        "valency": (2,),
        "electronegativity": 1.65,
        "oxidation states": (-2, 0, 1, 2),
    },
    "Ga": {
        "element": "gallium",
        "atomic number": 31,
        "atomic mass": 69.723,
        "atomic radius": 1.36,
        "valence electrons": 3,
        "valency": (3, 1),
        "electronegativity": 1.81,
        "oxidation states": (-5, -4, -3, -2, -1, 1, 2, 3),
    },
    "Ge": {
        "element": "germanium",
        "atomic number": 32,
        "atomic mass": 72.64,
        "atomic radius": 1.25,
        "valence electrons": 4,
        "valency": (4, 2),
        "electronegativity": 2.01,
        "oxidation states": (-4, -3, -2, -1, 0, 1, 2, 3, 4),
    },
    "As": {
        "element": "arsenic",
        "atomic number": 33,
        "atomic mass": 74.922,
        "atomic radius": 1.14,
        "valence electrons": 5,
        "valency": (5, 3),
        "electronegativity": 2.18,
        "oxidation states": (-3, -2, -1, 0, 1, 2, 3, 4, 5),
    },
    "Se": {
        "element": "selenium",
        "atomic number": 34,
        "atomic mass": 78.96,
        "atomic radius": 1.03,
        "valence electrons": 6,
        "valency": (2, 4, 6),
        "electronegativity": 2.55,
        "oxidation states": (-2, -1, 1, 2, 3, 4, 5, 6),
    },
    "Br": {
        "element": "bromine",
        "atomic number": 35,
        "atomic mass": 79.904,
        "atomic radius": 0.94,
        "valence electrons": 7,
        "valency": (1, 3, 5),
        "electronegativity": 2.96,
        "oxidation states": (-1, 1, 3, 4, 5, 7),
    },
    "Kr": {
        "element": "krypton",
        "atomic number": 36,
        "atomic mass": 83.709,
        "atomic radius": 0.87,
        "valence electrons": 8,
        "valency": (0,),
        "electronegativity": 3.0,
        "oxidation states": (0, 1, 2),
    },
    "Rb": {
        "element": "rubidium",
        "atomic number": 37,
        "atomic mass": 85.468,
        "atomic radius": 2.65,
        "valence electrons": 1,
        "valency": (1,),
        "electronegativity": 0.82,
        "oxidation states": (-1, 1),
    },
    "Sr": {
        "element": "strontium",
        "atomic number": 38,
        "atomic mass": 87.620,
        "atomic radius": 2.19,
        "valence electrons": 2,
        "valency": (2,),
        "electronegativity": 0.95,
        "oxidation states": (1, 2),
    },
    "Y": {
        "element": "yttrium",
        "atomic number": 39,
        "atomic mass": 88.906,
        "atomic radius": 2.12,
        "valence electrons": 3,
        "valency": (3,),
        "electronegativity": 1.22,
        "oxidation states": (0, 1, 2, 3),
    },
    "Zr": {
        "element": "zirconium",
        "atomic number": 40,
        "atomic mass": 91.224,
        "atomic radius": 2.06,
        "valence electrons": 4,
        "valency": (4,),
        "electronegativity": 1.33,
        "oxidation states": (-2, 0, 1, 2, 3, 4),
    },
    "Nb": {
        "element": "niobium",
        "atomic number": 41,
        "atomic mass": 92.906,
        "atomic radius": 1.98,
        "valence electrons": 5,
        "valency": (3, 5),
        "electronegativity": 1.6,
        "oxidation states": (-3, -1, 0, 1, 2, 3, 4, 5),
    },
    "Mo": {
        "element": "molybdenum",
        "atomic number": 42,
        "atomic mass": 95.95,
        "atomic radius": 1.90,
        "valence electrons": 6,
        "valency": (2, 3, 4, 5, 6),
        "electronegativity": 2.16,
        "oxidation states": (-4, -2, -1, 0, 1, 2, 3, 4, 5, 6),
    },
    "Tc": {
        "element": "technetium",
        "atomic number": 43,
        "atomic mass": 98,
        "atomic radius": 1.83,
        "valence electrons": 7,
        "valency": (4,),
        "electronegativity": 1.9,
        "oxidation states": (-3, -1, 0, 1, 2, 3, 4, 5, 6, 7),
    },
    "Ru": {
        "element": "ruthenium",
        "atomic number": 44,
        "atomic mass": 101.07,
        "atomic radius": 1.78,
        "valence electrons": 8,
        "valency": (3,),
        "electronegativity": 2.2,
        "oxidation states": (-4, -2, 0, 1, 2, 3, 4, 5, 6, 7, 8),
    },
    "Rh": {
        "element": "rhodium",
        "atomic number": 45,
        "atomic mass": 102.91,
        "atomic radius": 1.73,
        "valence electrons": 9,
        "valency": (3,),
        "electronegativity": 2.28,
        "oxidation states": (-3, -1, 0, 1, 2, 3, 4, 5, 6),
    },
    "Pd": {
        "element": "palladium",
        "atomic number": 46,
        "atomic mass": 106.42,
        "atomic radius": 1.69,
        "valence electrons": 10,
        "valency": (4,),
        "electronegativity": 2.2,
        "oxidation states": (0, 1, 2, 3, 4),
    },
    "Ag": {
        "element": "silver",
        "atomic mass": 47,
        "atomic number": 107.87,
        "atomic radius": 1.65,
        "valence electrons": 1,
        "valency": (1,),
        "electronegativity": 1.93,
        "oxidation states": (-2, -1, 1, 2, 3),
    },
    "Cd": {
        "element": "cadmium",
        "atomic number": 48,
        "atomic mass": 112.41,
        "atomic radius": 1.61,
        "valence electrons": 2,
        "valency": (2,),
        "electronegativity": 1.69,
        "oxidation states": (-2, 1, 2),
    },
    "In": {
        "element": "indium",
        "atomic number": 49,
        "atomic mass": 114.82,
        "atomic radius": 1.56,
        "valence electrons": 3,
        "valency": (3,),
        "electronegativity": 1.78,
        "oxidation states": (-5, -2, -1, 1, 2, 3),
    },
    "Sn": {
        "element": "tin",
        "atomic number": 50,
        "atomic mass": 118.71,
        "atomic radius": 1.45,
        "valence electrons": 4,
        "valency": (2, 4),
        "electronegativity": 1.96,
        "oxidation states": (-4, -3, -2, -1, 0, 1, 2, 3, 4),
    },
    "Sb": {
        "element": "antimony",
        "atomic number": 51,
        "atomic mass": 121.76,
        "atomic radius": 1.33,
        "valence electrons": 5,
        "valency": (3, 5),
        "electronegativity": 2.05,
        "oxidation states": (-3, -2, -1, 0, 1, 2, 3, 4, 5),
    },
    "Te": {
        "element": "tellurium",
        "atomic number": 52,
        "atomic mass": 127.6,
        "atomic radius": 1.23,
        "valence electrons": 6,
        "valency": (2, 4, 6),
        "electronegativity": 2.1,
        "oxidation states": (-2, -1, 1, 2, 3, 4, 5, 6),
    },
    "I": {
        "element": "iodine",
        "atomic number": 53,
        "atomic mass": 126.9,
        "atomic radius": 1.15,
        "valence electrons": 7,
        "valence": (1, 3, 5, 7),
        "electronegativity": 2.66,
        "oxidation states": (-1, 1, 3, 4, 5, 6, 7),
    },
    "Xe": {
        "element": "xenon",
        "atomic number": 54,
        "atomic mass": 131.29,
        "atomic radius": 1.08,
        "valence electrons": 8,
        "valency": (0,),
        "electronegativity": 2.6,
        "oxidation states": (0, 1, 2, 4, 6, 8),
    },
    "Cs": {
        "element": "cesium",
        "atomic number": 55,
        "atomic mass": 132.91,
        "atomic radius": 2.98,
        "valence electrons": 1,
        "valency": (1,),
        "electronegativity": 0.79,
        "oxidation states": (-1, 1),
    },
    "Ba": {
        "element": "barium",
        "atomic number": 56,
        "atomic mass": 137.33,
        "atomic radius": 2.53,
        "valence electrons": 2,
        "valency": (2,),
        "electronegativity": 0.89,
        "oxidation states": (1, 2),
    },
    "La": {
        "element": "lanthanum",
        "atomic number": 57,
        "atomic mass": 138.91,
        "atomic radius": 2.50,
        "valence electrons": 3,
        "valency": (3,),
        "electronegativity": 1.1,
        "oxidation states": (0, 1, 2, 3),
    },
    "Ce": {
        "element": "cerium",
        "atomic number": 58,
        "atomic mass": 140.12,
        "atomic radius": 2.48,
        "valence electrons": 4,
        "valency": (4,),
        "electronegativity": 1.12,
        "oxidation states": (1, 2, 3, 4),
    },
    "Pr": {
        "element": "praseodymium",
        "atomic number": 59,
        "atomic mass": 140.91,
        "atomic radius": 2.47,
        "valence electrons": 5,
        "valency": (3,),
        "electronegativity": 1.13,
        "oxidation states": (0, 1, 2, 3, 4, 5),
    },
    "Nd": {
        "element": "neodymium",
        "atomic number": 60,
        "atomic mass": 144.24,
        "atomic radius": 2.06,
        "valence electrons": 6,
        "valency": (3,),
        "electronegativity": 1.14,
        "oxidation states": (0, 2, 3, 4),
    },
    "Pm": {
        "element": "promethium",
        "atomic number": 61,
        "atomic mass": 145,
        "atomic radius": 2.05,
        "valence electrons": 7,
        "valency": (3,),
        "electronegativity": 1.13,
        "oxidation states": (2, 3),
    },
    "Sm": {
        "element": "samarium",
        "atomic number": 62,
        "atomic mass": 150.36,
        "atomic radius": 2.38,
        "valence electrons": 8,
        "valency": (3,),
        "electronegativity": 1.17,
        "oxidation states": (0, 2, 3),
    },
    "Eu": {
        "element": "europium",
        "atomic number": 63,
        "atomic mass": 151.96,
        "atomic radius": 2.31,
        "valence electrons": 9,
        "valency": (3,),
        "electronegativity": 1.2,
        "oxidation states": (0, 2, 3),
    },
    "Gd": {
        "element": "gadolinium",
        "atomic number": 64,
        "atomic mass": 157.25,
        "atomic radius": 2.33,
        "valence electrons": 10,
        "valency": (3,),
        "electronegativity": 1.2,
        "oxidation states": (0, 1, 2, 3),
    },
    "Tb": {
        "element": "terbium",
        "atomic number": 65,
        "atomic mass": 158.93,
        "atomic radius": 2.25,
        "valence electrons": 11,
        "valency": (3,),
        "electronegativity": 1.1,
        "oxidation states": (0, 1, 2, 3, 4),
    },
    "Dy": {
        "element": "dysprosium",
        "atomic number": 66,
        "atomic mass": 162.5,
        "atomic radius": 2.28,
        "valence electrons": 12,
        "valency": (3,),
        "electronegativity": 1.22,
        "oxidation states": (0, 1, 2, 3, 4),
    },
    "Ho": {
        "element": "holmium",
        "atomic number": 67,
        "atomic mass": 164.93,
        "atomic radius": 2.26,
        "valence electrons": 13,
        "valency": (3,),
        "electronegativity": 1.23,
        "oxidation states": (0, 1, 2, 3),
    },
    "Er": {
        "element": "erbium",
        "atomic number": 68,
        "atomic mass": 167.26,
        "atomic radius": 2.26,
        "valence electrons": 14,
        "valency": (3,),
        "electronegativity": 1.24,
        "oxidation states": (0, 1, 2, 3),
    },
    "Tm": {
        "element": "thulium",
        "atomic number": 69,
        "atomic mass": 168.93,
        "atomic radius": 2.22,
        "valence electrons": 15,
        "valency": (3,),
        "electronegativity": 1.25,
        "oxidation states": (0, 2, 3),
    },
    "Yb": {
        "element": "ytterbium",
        "atomic number": 70,
        "atomic mass": 173.04,
        "atomic radius": 2.22,
        "valence electrons": 16,
        "valency": (3,),
        "electronegativity": 1.1,
        "oxidation states": (0, 1, 2, 3),
    },
    "Lu": {
        "element": "lutetium",
        "atomic number": 71,
        "atomic mass": 174.97,
        "atomic radius": 2.17,
        "valence electrons": 3,
        "valency": (3,),
        "electronegativity": 1.27,
        "oxidation states": (0, 1, 2, 3),
    },
    "Hf": {
        "element": "hafnium",
        "atomic number": 72,
        "atomic mass": 178.49,
        "atomic radius": 2.08,
        "valence electrons": 4,
        "valency": (4,),
        "electronegativity": 1.3,
        "oxidation states": (-2, 0, 1, 2, 3, 4),
    },
    "Ta": {
        "element": "tantalum",
        "atomic number": 73,
        "atomic mass": 180.95,
        "atomic radius": 2.00,
        "valence electrons": 5,
        "valency": (5,),
        "electronegativity": 1.5,
        "oxidation states": (-3, -1, 0, 1, 2, 3, 4, 5),
    },
    "W": {
        "element": "tungsten",
        "atomic number": 74,
        "atomic mass": 183.84,
        "atomic radius": 1.93,
        "valence electrons": 6,
        "valency": (6,),
        "electronegativity": 2.36,
        "oxidation states": (-4, -2, -1, 0, 1, 2, 3, 4, 5, 6),
    },
    "Re": {
        "element": "rhenium",
        "atomic number": 75,
        "atomic mass": 186.21,
        "atomic radius": 1.88,
        "valence electrons": 7,
        "valency": (7,),
        "electronegativity": 1.9,
        "oxidation states": (-3, -1, 0, 1, 2, 3, 4, 5, 6, 7),
    },
    "Os": {
        "element": "osmium",
        "atomic number": 76,
        "atomic mass": 190.23,
        "atomic radius": 1.85,
        "valence electrons": 8,
        "valency": (2, 3, 4),
        "electronegativity": 2.2,
        "oxidation states": (-4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8),
    },
    "Ir": {
        "element": "iridium",
        "atomic number": 77,
        "atomic mass": 192.22,
        "atomic radius": 1.80,
        "valence electrons": 9,
        "valency": (2, 3, 4),
        "electronegativity": 2.2,
        "oxidation states": (-3, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
    },
    "Pt": {
        "element": "platinum",
        "atomic number": 78,
        "atomic mass": 195.08,
        "atomic radius": 1.77,
        "valence electrons": 10,
        "valency": (2, 4),
        "electronegativity": 2.28,
        "oxidation states": (-3, -2, -1, 0, 1, 2, 3, 4, 5, 6),
    },
    "Au": {
        "element": "gold",
        "atomic number": 79,
        "atomic mass": 196.97,
        "atomic radius": 1.74,
        "valence electrons": 1,
        "valency": (1,),
        "electronegativity": 2.54,
        "oxidation states": (-3, -2, -1, 0, 1, 2, 3, 5),
    },
    "Hg": {
        "element": "mercury",
        "atomic number": 80,
        "atomic mass": 200.59,
        "atomic radius": 1.71,
        "valence electrons": 2,
        "valency": (1, 2),
        "electronegativity": 2.0,
        "oxidation states": (-2, 1, 2),
    },
    "Tl": {
        "element": "thallium",
        "atomic number": 81,
        "atomic mass": 204.38,
        "atomic radius": 1.56,
        "valence electrons": 3,
        "valency": (1, 3),
        "electronegativity": 1.62,
        "oxidation states": (-5, -2, -1, 1, 2, 3),
    },
    "Pb": {
        "element": "lead",
        "atomic number": 82,
        "atomic mass": 207.20,
        "atomic radius": 1.54,
        "valence electrons": 4,
        "valency": (2, 4),
        "electronegativity": 1.87,
        "oxidation states": (-4, -2, -1, 1, 2, 3, 4),
    },
    "Bi": {
        "element": "bismuth",
        "atomic number": 83,
        "atomic mass": 208.98,
        "atomic radius": 1.43,
        "valence electrons": 5,
        "valency": (3, 5),
        "electronegativity": 2.02,
        "oxidation states": (-3, -2, -1, 1, 2, 3, 4, 5),
    },
}
