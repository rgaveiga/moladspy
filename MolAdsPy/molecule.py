from __future__ import print_function
from ._atomcollection import AtomCollection
from ._utils import apply_owner_handlers, update_owner_attributes
from ._exception import BasicException
from numpy import array, ndarray, min, max, sum, cos, sin, radians, dot
from copy import deepcopy
from multipledispatch import dispatch


class MoleculeError(BasicException):
    pass


class Molecule(AtomCollection):
    __slots__ = [
        "_vaccuum",
        "_maxx",
        "_maxy",
        "_maxz",
        "_minx",
        "_miny",
        "_minz",
        "_anchors",
    ]

    @dispatch(str, vaccuum=(int, float))
    def __init__(self, label, vaccuum=10.0, **kwargs):
        """
        Object initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the type of molecule, for instance,
            "benzene" or "C6H12".
        vaccuum : float, optional
            Vaccuum separating the molecule from its images in a periodic
            boundary conditions scheme. The default is 10.0 Angstroms.

        """
        super().__init__(**kwargs)

        if len(label) > 0:
            self._label = label
        else:
            raise MoleculeError("'label' must be a non-empty string!")

        if vaccuum > 0.0:
            self._vaccuum = vaccuum
        else:
            raise MoleculeError("'vaccuum' must be a number greater than zero!")

        self._maxx = self._maxy = self._maxz = None  # Maximum Cartesian coordinates
        self._minx = self._miny = self._minz = None  # Minimum Cartesian coordinates
        self._anchors = {
            "com": array([0.0, 0.0, 0.0])
        }  # Anchor points for translations and rotations

    @dispatch(str, str, file_type=str, vaccuum=(int, float))
    def __init__(self, label, file_name, file_type="XYZ", vaccuum=10.0, **kwargs):
        """
        Object initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the type of molecule, for instance,
            "benzene" or "C6H12".
        file_name : string
            Name of the file containing the molecule structure.
        file_type : string
            Type of file containing the molecule coordinates. The default is
            "XYZ", which is the only type currently supported.
        vaccuum : float, optional
            Vaccuum separating the molecule from its images in a periodic
            boundary conditions scheme. The default is 10.0 Angstroms.

        """
        super().__init__(**kwargs)

        if len(label) > 0:
            self._label = label
        else:
            raise MoleculeError("'label' must be a non-empty string!")

        self._anchors = {
            "com": array([0.0, 0.0, 0.0])
        }  # Anchor points for translations and rotations

        if vaccuum > 0.0:
            self._vaccuum = vaccuum
        else:
            raise MoleculeError("'vaccuum' must be a number greater than zero!")

        if len(file_name) > 0:
            if file_type.lower() == "xyz":
                self._read_xyz(file_name)
            else:
                raise MoleculeError("'file_type' must be 'XYZ'!")
        else:
            raise MoleculeError("'file_name' must be a valid file name!")
            
        self._update(**kwargs)

    @dispatch(str, list, vaccuum=(int, float))
    def __init__(self, label, atom_list, vaccuum=10.0, **kwargs):
        """
        Object initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the type of molecule, for instance,
            "benzene" or "C6H12".
        atom_list : Python list
            List of Atom objects or atom IDs to be added to the molecule.
        vaccuum : float, optional
            Vaccuum separating the molecule from its images in a periodic
            boundary conditions scheme. The default is 10.0 Angstroms.

        """
        super().__init__(**kwargs)

        if len(label) > 0:
            self._label = label
        else:
            raise MoleculeError("'label' must be a non-empty string!")

        self._anchors = {
            "com": array([0.0, 0.0, 0.0])
        }  # Anchor points for translations and rotations
        self._origin = None

        if vaccuum > 0.0:
            self._vaccuum = vaccuum
        else:
            raise MoleculeError("'vaccuum' must be a number greater than zero!")

        self._get_from_atom_list(atom_list, **kwargs)

        self._update(**kwargs)

    def add_anchor(self, anchor, pos):
        """
        add_anchor(anchor,pos) -> adds a new anchor point that can be used as
            the reference point of the molecule for translations and rotations.

        Parameters
        ----------
        anchor : string
            Name of the anchor point to be added to the molecule.
        pos : Numpy array
            Cartesian coordinates of the anchor point. It can also be provided
            as a Python list or tuple.

        """
        if isinstance(pos, (list, tuple)):
            pos = array(pos)

        if isinstance(anchor, str) and len(anchor) > 0:
            if isinstance(pos, ndarray) and pos.shape[0] == 3:
                self._anchors[anchor] = pos.astype(float)
            else:
                raise MoleculeError("'pos' must be an array with three components!")
        else:
            raise MoleculeError("'anchor' must be a non-empty string!")

    def remove_anchor(self, anchor):
        """
        remove_anchor(anchor) -> removes an anchor point.

        Parameters
        ----------
        anchor : string
            Name of the anchor point to be removed from the molecule.

        """
        if isinstance(anchor, str) and anchor in self._anchors.keys():
            if anchor == "com":
                raise MoleculeError("The center of mass cannot be removed!")
            else:
                self._anchors.pop(anchor)
        else:
            raise MoleculeError("'anchor' is not a valid anchor point!")

    def write_xyz(self, file_name="molecule.xyz", **kwargs):
        """
        Saves the atomic coordinates of the molecule into an XYZ file.

         Parameters
         ----------
         file_name : string, optional
             Name of the XYZ file. The default is "molecule.xyz".

        """       
        super().write_xyz(file_name, ucell = True, **kwargs)

    def write_pw_input(self, file_name="molecule.in", pseudopotentials={}, pwargs={}, **kwargs):
        """
        Creates a basic input file for geometry relaxation of the molecule using
        the pw.x code found in the Quantum Espresso package.

        Parameters
        ----------
        file_name : string, optional
            Name of the input file. The default is "molecule.in".
        pseudopotentials : Python dictionary, optional
            Specify for each element provided as a key the name of a pseudopotential
            file given as the corresponding value. The default is {}.
        pwargs : Python dictionary.
            Dictionary containing key-value pairs that allow some customization
            of the input file. For additional information, check out 'pw.x'
            documentation. These are currently the accepted keys:
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

        """
        pwargs["calculation"] = "relax"
        pwargs["ibrav"] = 0
        pwargs["kvec"] = [1, 1, 1, 0, 0, 0]

        super().write_pw_input(file_name, pseudopotentials, pwargs, ucell = True, **kwargs)

    def copy(self):
        """
        Returns a copy of the Molecule object.

        Returns
        -------
        Molecule object
            Copy of the Molecule object.
        """
        newmol = super().copy()
        newmol._anchors = deepcopy(self._anchors)

        return newmol

    def displace(self, disp, **kwargs):
        """
        displace(disp) -> rigidly displaces the molecule.

        Parameters
        ----------
        disp : Numpy array
            Displacement vector. It can also be provided as a Python list or tuple.
        """        
        super().displace(disp, **kwargs)

        for key in self._anchors.keys():
            if key != "com":
                self._anchors[key] += disp

    def move_to(self, x, anchor="com", **kwargs):
        """
        move_to(x,anchor) -> rigidly moves the molecule such that the anchor
        point is located at 'x'.

        Parameters
        ----------
        x : Numpy array
            New position of the molecule's anchor point. It can also be provided
            as a Python list or tuple.
        anchor : string, optional
            Anchor point used as reference for translations. The default is
            'com', which means the molecule's center of mass.

        """       
        if isinstance(x, (list, tuple)):
            x = array(x)

        if not (isinstance(anchor, str) and anchor in self._anchors.keys()):
            raise MoleculeError("'anchor' is not a valid anchor point!")

        if isinstance(x, ndarray) and x.shape[0] == 3:
            disp = x.astype(float) - self._anchors[anchor]

            self.displace(disp, **kwargs)
        else:
            raise MoleculeError("'x' must be an array with three components!")

    @apply_owner_handlers
    def rotate(self, theta, phi, psi, anchor="com", **kwargs):
        """
        rotate(theta,phi,psi,anchor) -> rotates the molecule around an anchor point.

        Parameters
        ----------
        theta : float, optional
            First rotation angle (around Y axis), in degrees.
        phi : float, optional
            Second rotation angle (around X axis), in degrees.
        psi : float, optional
            Third rotation angle (around Z axis), in degrees.
        anchor : string, optional
            Anchor point around which the molecule will be rotated. The default
            is 'com', which means the molecule's center of mass.

        """
        if (
            isinstance(theta, (float, int))
            and isinstance(phi, (float, int))
            and isinstance(psi, (float, int))
        ):
            theta = radians(theta)
            phi = radians(phi)
            psi = radians(psi)
        else:
            raise MoleculeError("Rotation angles must be provided as numbers!")

        if not (isinstance(anchor, str) and anchor in self._anchors.keys()):
            raise MoleculeError("'anchor' is not a valid anchor point!")

        a11 = cos(theta) * cos(psi)
        a12 = -cos(phi) * sin(psi) + sin(phi) * sin(theta) * cos(psi)
        a13 = sin(phi) * sin(psi) + cos(phi) * sin(theta) * cos(psi)
        a21 = cos(theta) * sin(psi)
        a22 = cos(phi) * cos(psi) + sin(phi) * sin(theta) * sin(psi)
        a23 = -sin(phi) * cos(psi) + cos(phi) * sin(theta) * sin(psi)
        a31 = -sin(theta)
        a32 = sin(phi) * cos(theta)
        a33 = cos(phi) * cos(theta)
        rotmatrix = array([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
        rotpoint = self._anchors[anchor]

        for atom in self._atoms:
            atom._x = dot((atom._x - rotpoint), rotmatrix) + rotpoint

        for key in self._anchors.keys():
            if key != anchor:
                self._anchors[key] = (
                    dot((self._anchors[key] - rotpoint), rotmatrix) + rotpoint
                )

        self._update(**kwargs)

    def resize(self):
        """
        Does nothing, because this method is not implemented for Molecule objects.

        """
        raise MoleculeError("Method not available for Molecule objects!")

    @update_owner_attributes
    def _update(self, **kwargs):
        """
        _update() -> simultaneously updates the molecule's center of mass and
        the values of its extremities.

        """
        valid_coords = array([atom._x for atom in self.active_atoms])
        atomic_mass = array([atom.atomic_mass for atom in self.active_atoms])

        if valid_coords.shape[0] > 0:
            self._maxx, self._maxy, self._maxz = max(valid_coords, axis=0)
            self._minx, self._miny, self._minz = min(valid_coords, axis=0)
            totmass = sum(atomic_mass)
            self._anchors["com"] = (
                sum(atomic_mass * valid_coords.transpose(), axis=1) / totmass
            )
            self._latvec = array(
                [
                    [self._maxx - self._minx + self._vaccuum, 0.0, 0.0],
                    [0.0, self._maxy - self._miny + self._vaccuum, 0.0],
                    [0.0, 0.0, self._maxz - self._minz + self._vaccuum],
                ]
            )
            self._origin = array([self._minx, self._miny, self._minz])
        else:
            self._maxx = self._maxy = self._maxz = None
            self._minx = self._miny = self._minz = None
            self._anchors["com"] = array([0.0, 0.0, 0.0])
            self._latvec = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            self._origin = array([0.0, 0.0, 0.0])

    def __str__(self):
        """
        __str__() -> returns the type and the ID of the Molecule object.

        Returns
        -------
        String.
            A string containing the type and ID of the Molecule object.

        """
        return "<Molecule object> Type: %s; ID: %d" % (self._label, self._id)

    """
    Properties
    ----------
    a0 : float, readonly
        Lattice parameter. Meaningless for molecules, it is still necessary  
        when writing extended XYZ files or PW input files. It is set to 1.
    lattice_vectors : Numpy array, readonly
        Lattice vectors. Meaningless for molecules, it is still necessary  
        when writing extended XYZ files or PW input files. The vectors are 
        determined by adding the prescribed length of vaccuum to the extremities 
        of the molecule.
    anchors : Python dictionary, readonly
        Dictionary of anchor points.
    center_of_mass : Numpy array, readonly
        Center of mass of the molecule.
    maxx,minx,maxy,miny,maxz,minz : floats, readonly
        Molecule extremities.
    origin : Numpy array
        Position with minimum values of X, Y, and Z coordinates of molecule.
    """

    @property
    def a0(self):
        return self._a0

    @property
    def lattice_vectors(self):
        return self._latvec

    @property
    def anchors(self):
        return self._anchors

    @property
    def center_of_mass(self):
        return self._anchors["com"]

    @property
    def maxx(self):
        return self._maxx

    @property
    def maxy(self):
        return self._maxy

    @property
    def maxz(self):
        return self._maxz

    @property
    def minx(self):
        return self._minx

    @property
    def miny(self):
        return self._miny

    @property
    def minz(self):
        return self._minz

    @property
    def origin(self):
        return self._origin
