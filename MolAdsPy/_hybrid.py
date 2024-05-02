from __future__ import print_function
from ._atomcollection import AtomCollection
from ._exception import BasicException
from abc import ABC, abstractmethod
from numpy import array


class HybridError(BasicException):
    pass


class Hybrid(AtomCollection, ABC):
    __slots__ = [
        "_components",
        "_a0",
        "_latvec",
        "_n",
        "_m",
        "_l",
        "_maxn",
        "_maxm",
        "_maxl",
        "_origin",
    ]

    def __init__(self, **kwargs):
        """
        Object initialization.

        """
        self._components = (
            {}
        )  # Dictionary of objects that are put together to form the hybrid system
        self._a0 = 1.0  # Lattice constant of the hybrid object
        self._latvec = array(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        )  # Lattice vectors of the hybrid object
        self._n = 0  # Number of translations along the first lattice vector
        self._m = 0  # Number of translations along the second lattice vector
        self._l = 0  # Number of translations along the third lattice vector
        self._maxn = 0  # Maximum number of translations already perfomed along the first lattice vector
        self._maxm = 0  # Maximum number of translations already perfomed along the second lattice vector
        self._maxl = 0  # Maximum number of translations already perfomed along the third lattice vector
        self._origin = array(
            [0.0, 0.0, 0.0]
        )  # Position with minimum values of X, Y, and Z coordinates of the structure.

        self._verbose = kwargs.get("verbose", False)
        type(self)._append(self)

    @abstractmethod
    def add_component(self, *args, **kwargs):
        """
        To be implemented in a derived hybrid structure, this method is expected
        to add atom collection objects to the hybrid system.

        """
        pass

    @abstractmethod
    def remove_component(self, *args, **kwargs):
        """
        To be implemented in a derived hybrid structure, this method is expected
        to remove an atom collection object from the hybrid system.

        """
        pass

    def write_xyz(self, file_name="coords.xyz", **kwargs):
        """
        Saves the atomic coordinates of the hybrid system into an XYZ file.

        Parameters
        ----------
        file_name : string, optional
            Name of the XYZ file. The default is "coords.xyz".

        """
        super().write_xyz(file_name, ucell=False, **kwargs)

    def write_pw_input(
        self, file_name="pw.in", pseudopotentials={}, pwargs={}, **kwargs
    ):
        """
        Creates a basic input file for geometry relaxation of the structure using
        the pw.x code found in the Quantum Espresso package.

        Parameters
        ----------
        file_name : string, optional
            Name of the input file. The default is "pw.in".
        ucell : logical, optional
            Write only the coordinates of atoms in the unit cell. The default is
            False.
        pseudopotentials : Python dictionary, optional
            Specify for each element provided as a key the name of a pseudopotential
            file given as the corresponding value. The default is {}.
        pwargs : Python dictionary.
            Dictionary containing key-value pairs that allow some customization
            of the input file. For additional information, check out 'pw.x'
            documentation. These are currently the accepted keys:
                calculation : string, optional
                    Type of calculation to be carried out. It must be either
                    'relax' or 'vc-relax'. The default is 'relax'.
                ibrav : integer, optional
                    Bravais lattice. The default is 0.
                celldm(2)-celldm(6) : float, optional (depending on 'ibrav')
                    Crystallographic constants that can be provided depending on
                    'ibrav' value. 'celldm(1)', the lattice parameter, is calculated
                    from the 'a0' property of the atom collection object.
                ecutwfc : float, optional
                    Plane-wave energy cutoff. The default is 32 Ry.
                ecutrho : float, optional
                    Charge density cutoff. The default is 128 Ry.
                nspin : integer, optional
                    Spin polarization. It can be either 1 (non-polarized) or
                    2 (polarized, magnetization along z axis). The default is 1.
                occupations : string, optional
                    Occupation function. The default is 'fixed'.
                smearing : string, optional
                    Smearing for metals, if the occupation function is 'smearing'.
                    The defautl is 'gaussian'.
                degauss : float, optional
                    Gaussian spreading, in Ry, if the occupation function is
                    'smearing'. The default is 0.02 Ry.
                starting_magnetization : Python list, optional
                    List of starting magnetic moments for the different atom
                    types in the slab. The default is 1.0 Bohr magneton
                    for the first atom type and 0.0 for the others, if any.
                input_dft : string, optional
                    It allows to activate non-local vdW interaction by setting
                    it to 'vdw-df'. The default is ''.
                mixing_mode : string ,optional
                    Mixing charge mode. It can be 'plain' (Broyden mixing),
                    'TF' (Thomas-Fermi screening for homogeneous systems) or
                    'local-TF' (local density Thomas-Fermi screening). The
                    default is 'plain'.
                mixing_beta : float, optional
                    Beta parameter in charge mixing. The default is 0.7.
                kvec : Python list, optional
                    Grid (first three elements) and shift (last three
                    components) used to generate k-points according to the
                    Monhorst-Pack scheme. The default is [1,1,1,0,0,0].
        """
        super().write_pw_input(
            file_name, pseudopotentials, pwargs, ucell=False, **kwargs
        )

    @abstractmethod
    def _begin_handler(self, obj, method, method_args, **kwargs):
        """
        To be implemented in a derived hybrid structure, this method will typically
        be used to check whether a specific operation performed by an object belonging
        to the hybrid object is permitted or not.

        """
        pass

    @abstractmethod
    def _end_handler(self, obj, method, method_args, **kwargs):
        """
        To be implemented in a derived hybrid structure, this method will typically
        be used to adjust the result of some operation performed by an object that
        belongs to the hybrid object.

        """
        pass

    @abstractmethod
    def _update(self, obj, **kwargs):
        """
        To be implemented in a derived hybrid structure, this method updates
        attributes of the hybrid object when an object that belongs to it undergoes
        changes in its own attributes.

        """
        pass

    def copy(self):
        """
        This method is disabled for hybrid objects and issues an error message.

        """
        raise HybridError("Not implemented for hybrid objects!")

    def add_atom(self, *args, **kwargs):
        """
        This method is disabled for hybrid objects and issues an error message.

        """
        raise HybridError("Not implemented for hybrid objects!")

    def remove_atom(self, *args, **kwargs):
        """
        This method is disabled for hybrid objects and issues an error message.

        """
        raise HybridError("Not implemented for hybrid objects!")

    def location(self, *args, **kwargs):
        """
        This method is disabled for hybrid objects and issues an error message.

        """
        raise HybridError("Not implemented for hybrid objects!")

    def _read_xyz(self, *args, **kwargs):
        """
        This method is disabled for hybrid objects and issues an error message.

        """
        raise HybridError("Not implemented for hybrid objects!")

    def __setitem__(self, *args, **kwargs):
        """
        This method is disabled for hybrid objects and issues an error message.

        """
        raise HybridError("Not implemented for hybrid objects!")

    """
    Properties
    ----------
    _atoms : Python tuple, readonly
        Tuple of atoms in the Hybrid object.
    _species : Python tuple, readonly
        Tuple of atomic species in the Hybrid object.
    _loc : Python dictionary, readonly
        Dictionary of three indices, indicating the position of every atom in 
        the hybrid system supercell.
    """

    @property
    def _atoms(self):
        atom_list = []

        for key in self._components:
            if isinstance(self._components[key], (list, tuple)):
                for component in self._components[key]:
                    atom_list += component._atoms
            else:
                atom_list += self._components[key]._atoms

        return tuple(atom_list)

    @property
    def _species(self):
        spec_list = []

        for key in self._components:
            if isinstance(self._components[key], (list, tuple)):
                for component in self._components[key]:
                    spec_list += component._species
            else:
                spec_list += self._components[key]._species

        return tuple(set(spec_list))

    @property
    def _loc(self):
        loc_dict = {}

        for key in self._components:
            if isinstance(self._components[key], (list, tuple)):
                for component in self._components[key]:
                    loc_dict.update(component._loc)
            else:
                loc_dict.update(self._components[key]._loc)

        return loc_dict
