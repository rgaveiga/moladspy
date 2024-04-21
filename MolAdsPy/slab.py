from __future__ import print_function
from ._atomcollection import AtomCollection
from ._exception import BasicException
from .atom import Atom
from numpy import array,ndarray,min,max
from copy import deepcopy
from multipledispatch import dispatch


class SlabError(BasicException):
    pass


class Slab(AtomCollection):
    __slots__ = ["_vaccuum", "_ads_sites", "_top", "_bottom"]

    @dispatch(str, vaccuum=(int, float))
    def __init__(self, label, vaccuum=10.0, **kwargs):
        """
        Object initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the slab, e.g., "graphene".
        vaccuum : float, optional
            Vaccuum separating the slab top from the slab bottom in a periodic
            boundary conditions scheme. The default is 10.0 Angstroms.

        """
        #TODO: Check in Molecule, same question
        super().__init__(**kwargs)

        if len(label) > 0:
            self._label = label
        else:
            raise SlabError("'label' must be a non-empty string!")

        if vaccuum > 0.0:
            self._vaccuum = vaccuum
        else:
            raise SlabError("'vaccuum' must be a number greater than zero!")

        self._ads_sites = {}  # Dictionary of adsorption sites
        self._top = None  # Maximum Z coordinate of the slab
        self._bottom = None  # Miniumum Z coordinate of the slab

    @dispatch(str, str, file_type=str, vaccuum=(int, float))
    def __init__(self, label, file_name, file_type="XYZ", vaccuum=10.0):
        """
        Object initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the slab, e.g., "graphene".
        file_name : string
            Name of the file containing the slab structure.
        file_type : string
            Type of file containing the slab coordinates. The default is "XYZ",
            which is the only type currently supported.
        vaccuum : float, optional
            Vaccuum separating the slab top from the slab bottom in a periodic
            boundary conditions scheme. The default is 10.0 Angstroms.

        """
        super().__init__()

        if len(label) > 0:
            self._label = label
        else:
            raise SlabError("'label' must be a non-empty string!")

        if vaccuum > 0.0:
            self._vaccuum = vaccuum
        else:
            raise SlabError("'vaccuum' must be a number greater than zero!")

        if file_type.lower() == "xyz":
            if len(file_name) > 0:
                self._read_xyz(file_name)
            else:
                raise SlabError("'file_name' must be a valid file name!")
        else:
            raise SlabError("'file_type' must be 'XYZ'!")

        self._ads_sites = {}  # Dictionary of adsorption sites

    @dispatch(str, list, float, (ndarray, list, tuple), vaccuum=(int, float))
    def __init__(self, label, atom_list, a0, lattice_vectors, vaccuum=10.0):
        """
        Object  initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the slab, e.g., "graphene".
        atom_list : Python list
            List of Atom objects to be added to the slab's unit cell.
        a0 : float
            Slab's lattice constant.
        lattice_vectors : Python array
            Slab's lattice vectors. Only the first and second lattice vectors
            are considered, although the three vectors must be provided. It can
            also be provided as a Python list or tuple.
        vaccuum : float, optional
            Vaccuum separating the slab top from the slab bottom in a periodic
            boundary conditions scheme. The default is 10.0 Angstroms.

        """
        super().__init__()

        if len(label) > 0:
            self._label = label
        else:
            raise SlabError("'label' must be a non-empty string!")

        if vaccuum > 0.0:
            self._vaccuum = vaccuum
        else:
            raise SlabError("'vaccuum' must be a number greater than zero!")

        if a0 > 0.0:
            self._a0 = a0
        else:
            raise SlabError("'a0' must be a positive number greater than zero!")

        if isinstance(lattice_vectors, (list, tuple)):
            lattice_vectors = array(lattice_vectors)

        if lattice_vectors.shape[0] == 3 and lattice_vectors.size == 9:
            self._latvec = lattice_vectors.astype(float)
        else:
            raise SlabError(
                "'lattice_vectors' must be a Numpy array with three vectors!"
            )

        if len(atom_list) > 0:
            for atom in atom_list:
                if isinstance(atom, (Atom, int)):
                    self.add_atom(atom, loc=(0, 0, 0), update=False)
                else:
                    if self._verbose:
                        print(
                            "WARNING! An element in the atom list must be either an Atom object or an atom ID!"
                        )
        else:
            raise SlabError("'atom_list' must be a non-empyt list!")

        self._update()

        self._ads_sites = {}  # Dictionary of adsorption sites

    def add_adsorption_site(self, ads_site, pos):
        """
        Adds an adsorption site on the slab.

        Parameters
        ----------
        ads_site : string
            Label of the adsorption site.
        pos : Numpy array
            Coordinates of the adsorption site on the XY plane. It can also be
            provided as a Python list or tuple.

        """
        if isinstance(pos, (list, tuple)):
            pos = array(pos)

        if isinstance(pos, ndarray) and pos.shape[0] == 2:
            self._ads_sites[ads_site] = pos.astype(float)
        else:
            raise SlabError("'pos' must be an array with two elements!")

    def remove_adsorption_site(self, ads_site):
        """
        Removes an adsorption site.

        Parameters
        ----------
        ads_site : string
            Label of the adsorption site.

        """
        if ads_site in self._ads_sites.keys():
            self._ads_sites.pop(ads_site)
        else:
            raise SlabError("'ads_site' is not a valid label for an adsorption site!")

    def displace(self, disp):
        """
        displace(disp) -> rigidly displaces the slab.

        Parameters
        ----------
        disp : Numpy array
            Displacement vector. It can also be provided as a Python list or tuple.
            
        """
        super().displace(disp)

        for key in self._ads_sites.keys():
            self._ads_sites[key] += disp[:2]

    def resize(self, n, m):
        """
        Resizes the slab in the XY plane.

        Parameters
        ----------
        n : integer
            Number of repetitions of the slab's unit cell along the first lattice
            vector.
        m : integer
            Number of repetitions of the slab's unit cell along the second lattice
            vector.

        """
        l = 0

        super().resize(n, m, l)

    def write_xyz(self, file_name="slab.xyz", ucell=False):
        """
        write_xyz(file_name,ucell) -> saves the atomic coordinates of the slab into an
            XYZ file.

        Parameters
        ----------
        file_name : string, optional
            Name of the XYZ file. The default is "slab.xyz".
        ucell : logical, optional
            Write only the coordinates of atoms in the unit cell. The default is
            False.

        """
        super().write_xyz(file_name, ucell)

    def write_pw_input(
        self, file_name="slab.in", ucell=False, pseudopotentials={}, pwargs={}
    ):
        """
        Creates a basic input file for geometry relaxation of the slab using the
        pw.x code found in the Quantum Espresso package.

        Parameters
        ----------
        file_name : string, optional
            Name of the input file. The default is "slab.in".
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
        super().write_pw_input(file_name, ucell, pseudopotentials, pwargs)

    def copy(self):
        """
        Returns a copy of the Slab object.

        Returns
        -------
        Slab object
            Copy of the Slab object.
            
        """
        newslab = super().copy()
        newslab._ads_sites = deepcopy(self._ads_sites)

        return newslab

    def _update(self):
        """
        _update() -> updates the slab's top and bottom, as well as its origin,
            if an event (e.g., adding a new atom to the slab) changes the atomic
            coordinates.

        """
        valid_coords = array(
            [
                atom._x
                for atom in self._atoms
                if atom._active
                and self._loc[atom._id][0] <= self._n
                and self._loc[atom._id][1] <= self._m
                and self._loc[atom._id][2] <= self._l
            ]
        )

        if valid_coords.shape[0] > 0:
            self._origin = min(valid_coords, axis=0)
            self._top = max(valid_coords[:, 2])
            self._bottom = min(valid_coords[:, 2])
            self._latvec[2] = array(
                [0.0, 0.0, (self._top - self._bottom + self._vaccuum) / self._a0]
            )
        else:
            self._origin = array([0.0, 0.0, 0.0])
            self._top = self._bottom = None
            self._latvec[2] = array([0.0, 0.0, 1.0])

    def __str__(self):
        """
        Returns the name and ID of the Slab object.

        Returns
        -------
        String.
            A string containing the name and ID of the Slab object.
            
        """
        return "<Slab object> Name: %s; ID: %d" % (self._label, self._id)

    """
    Properties
    ----------
    adsorption_sites : Python dictionary, readonly
        Dictionary of adsorption sites on the slab surface.
    top : float, readonly
        Top slab surface.
    bottom : float, readonly
        Bottom slab surface.
    vaccuum : float
        Empty space separating the top from the bottom of the slab, considering 
        periodic boundary conditions along the direction perpendicular to the 
        slab surfaces.
    """

    @property
    def adsorption_sites(self):
        return self._ads_sites

    @property
    def top(self):
        return self._top

    @property
    def bottom(self):
        return self._bottom

    @property
    def vaccuum(self):
        return self._vaccuum

    @vaccuum.setter
    def vaccuum(self, val):
        if isinstance(val, (float, int)) and val > 0.0:
            self._vaccuum = val
            self._latvec[2] = array(
                [0.0, 0.0, (self._top - self._bottom + val) / self._a0]
            )
        else:
            raise SlabError("Vaccuum must be provided as a number greater than zero!")
