from __future__ import print_function
from ._hybrid import Hybrid
from ._exception import BasicException
from .molecule import Molecule
from .slab import Slab
from numpy import array, dot, max, min, ndarray
from multipledispatch import dispatch


class MoleculeAdsorptionError(BasicException):
    pass


class MoleculeAdsorption(Hybrid):
    __slots__ = ["_substrate", "_minsep", "_vaccuum", "_side", "_top", "_bottom"]

    def __init__(self, substrate, minsep=1.0, vaccuum=10.0, **kwargs):
        """
        Object initialization.

        Parameters
        ----------
        substrate : Slab object
            Substrate where one or more molecules can be adsorbed.
        minsep : float, optional
            Minimum separation between the molecule and the substrate.
        vaccuum : float, optional
            Vaccuum separating the atom with the highest Z coordinate from the
            atom with the lowest Z coordinate in a periodic boundary conditions
            scheme. The default is 10.0 Angstroms.

        """
        super().__init__(**kwargs)

        if isinstance(substrate, Slab):
            if substrate._belongs_to is not None:
                raise MoleculeAdsorptionError(
                    "The substrate already belongs to another structure!"
                )

            self._components["substrate"] = substrate
            self._origin = self._components["substrate"]._origin
            self._latvec = self._components["substrate"]._latvec
            self._components["molecules@top"] = []
            self._components["molecules@bottom"] = []
            self._components["polymer@top"] = []
            self._components["polymer@bottom"] = []

            if isinstance(minsep, (int, float)) and minsep > 0.0:
                self._minsep = minsep
            else:
                raise MoleculeAdsorptionError("'minsep' must be a number greater than zero!")

            if isinstance(vaccuum, (int, float)) and vaccuum > 0.0:
                self._vaccuum = vaccuum
            else:
                raise MoleculeAdsorptionError("'vaccuum' must be a number greater than zero!")
        else:
            raise MoleculeAdsorptionError("Substrate must be a Slab object!")

        substrate._belongs_to = self

        self._top = substrate._top
        self._bottom = substrate._bottom

        self._update(substrate, **kwargs)


    @dispatch(Molecule, str, n=int, m=int, anchor=str, vertsep=(int, float), side=str)
    def add_component(
        self,
        molecule,
        ads_site,
        n=0,
        m=0,
        anchor="com",
        vertsep=1.0,
        side="top",
        **kwargs
    ):
        """
        Adsorbs a molecule onto the substrate.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        ads_site : string
            Label of the adsorption site.
        n : integer, optional
            Number of translations applied to the position of the adsorption site
            along the first lattice vector. The default is zero.
        m : integer, optional
            Number of translations applied to the position of the adsorption site
            along the second lattice vector. The default is zero.
        anchor : string, optional
            Anchor point in the molecule that will be vertically aligned with
            the adsorption site on the substrate. The default is "com", which
            is the positon of the molecule's center of mass.
        vertsep : float, optional
            Nearest vertical distance separating the molecule from the
            substrate. The default is 1.0 Angstrom.
        side : string, optional
            Side of the slab where the molecule will be adsorbed. It must be either
            "top" or "bottom" (self-explanatory). The default is "top".
        """
        if molecule._belongs_to is not None:
            raise MoleculeAdsorptionError("The molecule already belongs to another structure!")
        elif self._components["substrate"] is None:
            raise MoleculeAdsorptionError("No substrate has been defined!")
        elif len(anchor) == 0 or anchor not in molecule.anchors.keys():
            raise MoleculeAdsorptionError("'anchor' must be a valid anchor point!")
        elif n < 0:
            raise MoleculeAdsorptionError(
                "'n' must be an integer greater than or equal to zero!"
            )
        elif m < 0:
            raise MoleculeAdsorptionError(
                "'m' must be an integer greater than or equal to zero!"
            )

        if len(ads_site) > 0 and ads_site in self._components["substrate"]._ads_sites:
            sep = max([vertsep, self._minsep])

            if side == "top":
                dz = self._components["substrate"]._top + sep - molecule._minz
            elif side == "bottom":
                dz = self._components["substrate"]._bottom - sep - molecule._maxz
            else:
                raise MoleculeAdsorptionError("'side' must be either 'top' or 'bottom'!")

            self._components["molecules@" + side].append(molecule)

            molecule._belongs_to = self
            x0 = self._components["substrate"]._ads_sites[ads_site] + dot(
                self._a0 * array([n, m]), self._latvec[:2, :2]
            )
            x = array([x0[0], x0[1], molecule.anchors[anchor][2] + dz])

            molecule.move_to(x, anchor, **kwargs)
        else:
            raise MoleculeAdsorptionError("A valid adsorption site must be provided!")


    @dispatch(
        Molecule, (ndarray, list, tuple), anchor=str, vertsep=(int, float), side=str
    )
    def add_component(
        self, molecule, pos, anchor="com", vertsep=1.0, side="top", **kwargs
    ):
        """
        Adsorbs a molecule onto the substrate.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        pos : Numpy array
            XY coordinates in the substrate above which the molecule will be
            placed.
        anchor : string, optional
            Anchor point in the molecule that will be vertically aligned with
            the adsorption site on the substrate. The default is "com", which
            is the positon of the molecule's center of mass.
        vertsep : float, optional
            Nearest vertical distance separating the molecule from the
            substrate. The default is 1.0 Angstrom.
        side : string, optional
            Side of the slab where the molecule will be adsorbed. It must be either
            "top" or "bottom" (self-explanatory). The default is "top".
        """
        if molecule._belongs_to is not None:
            raise MoleculeAdsorptionError("The polymer already belongs to another structure!")
        elif self._components["substrate"] is None:
            raise MoleculeAdsorptionError("No substrate has been defined!")
        elif len(anchor) == 0 or anchor not in molecule.anchors.keys():
            raise MoleculeAdsorptionError("'anchor' must be a valid anchor point!")

        if isinstance(pos, (list, tuple)):
            pos = array(pos)

        if isinstance(pos, ndarray) and pos.shape[0] == 2:
            sep = max([vertsep, self._minsep])

            if side == "top":
                dz = self._components["substrate"]._top + sep - molecule._minz
            elif side == "bottom":
                dz = self._components["substrate"]._bottom - sep - molecule._maxz
            else:
                raise MoleculeAdsorptionError("'side' must be either 'top' or 'bottom'!")

            x = array([pos[0], pos[1], molecule.anchors[anchor][2] + dz])

            self._components["molecules@" + side].append(molecule)

            molecule._belongs_to = self

            molecule.move_to(x, anchor, **kwargs)
        else:
            raise MoleculeAdsorptionError(
                "'pos' must be provided as a list with two components!"
            )


    @dispatch(Molecule)
    def remove_component(self, molecule, **kwargs):
        """
        Removes an adsorbed molecule.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule to be removed.
        """
        if molecule in self._components["molecules@top"]:
            self._components["molecules@top"].remove(molecule)
        elif molecule in self._components["molecules@bottom"]:
            self._components["molecules@bottom"].remove(molecule)
        else:
            raise MoleculeAdsorptionError("Molecule not found!")

        molecule._belongs_to = None

        self._update(molecule, **kwargs)


    @dispatch(int)
    def remove_component(self, molid, **kwargs):
        """
        Removes an adsorbed molecule.

        Parameters
        ----------
        molid : integer
            ID of the adsorbed molecule.
        """
        for molecule in self.adsorbed_molecules:
            if molecule.ID == molid:
                if molecule in self._components["molecules@top"]:
                    self._components["molecules@top"].remove(molecule)
                elif molecule in self._components["molecules@bottom"]:
                    self._components["molecules@bottom"].remove(molecule)

                break
        else:
            raise MoleculeAdsorptionError("Molecule ID not found!")

        molecule._belongs_to = None

        self._update(molecule, **kwargs)


    @dispatch(Molecule, (int, float))
    def set_separation(self, molecule, sep, **kwargs):
        """
        Rigidly displaces the molecule along the direction perpendicular to the
        substrate, such that the separation between the top/bottom of the substrate
        and the nearest atom in the molecule be equal to 'sep'.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.

        """
        if sep < self._minsep:
            raise MoleculeAdsorptionError("Molecule will be too close to the substrate!")

        dz = 0.0

        if molecule in self._components["molecules@top"]:
            dz = self._components["substrate"]._top + sep - molecule._minz
        elif molecule in self._components["molecules@bottom"]:
            dz = self._components["substrate"]._bottom - sep - molecule._maxz
        else:
            raise MoleculeAdsorptionError("Molecule not found!")

        if abs(dz) > 0.0:
            molecule.displace(array([0.0, 0.0, dz]), call_handler=False, **kwargs)


    @dispatch(int, (int, float))
    def set_separation(self, molid, sep, **kwargs):
        """
        Rigidly displaces the molecule along the direction perpendicular to the
        substrate, such that the separation between the top/bottom of the substrate
        and the nearest atom in the molecule be equal to 'sep'.

        Parameters
        ----------
        molid : integer
            ID of the adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.
        """
        if sep < self._minsep:
            raise MoleculeAdsorptionError("Molecule will be too close to the substrate!")

        for molecule in self.adsorbed_molecules:
            if molecule.ID == molid:
                dz = 0.0

                if molecule in self._components["molecules@top"]:
                    dz = self._components["substrate"]._top + sep - molecule._minz
                elif molecule in self._components["molecules@bottom"]:
                    dz = self._components["substrate"]._bottom - sep - molecule._maxz
                else:
                    raise MoleculeAdsorptionError("Molecule not found!")

                if abs(dz) > 0.0:
                    molecule.displace(
                        array([0.0, 0.0, dz]), call_handler=False, **kwargs
                    )

                break
        else:
            raise MoleculeAdsorptionError("Molecule ID not found!")


    @dispatch(Slab, object, list)
    def _begin_handler(self, slab, method, method_args, **kwargs):
        pass

    @dispatch(Slab, object, list)
    def _end_handler(self, slab, method, method_args, **kwargs):
        pass

    @dispatch(Molecule, object, list)
    def _begin_handler(self, molecule, method, method_args, **kwargs):
        pass

    @dispatch(Molecule, object, list)
    def _end_handler(self, molecule, method, method_args, **kwargs):
        dz = 0.0

        if (
            molecule in self._components["molecules@top"]
            and (molecule.minz - self._components["substrate"]._top) < self._minsep
        ):
            dz = self._components["substrate"]._top + self._minsep - molecule._minz
        elif (
            molecule in self._components["molecules@bottom"]
            and (self._components["substrate"]._bottom - molecule._maxz) < self._minsep
        ):
            dz = self._components["substrate"]._bottom - self._minsep - molecule._maxz

        if abs(dz) > 0.0:
            molecule.displace(array([0.0, 0.0, dz]), call_handler=False)


    def _update(self, obj=None, **kwargs):
        """
        Updates attributes of the MoleculeAdsorption object when atomic coordinates or
        cell shape or cell size change.

        Parameters
        ----------
        obj : AtomCollection object, optional
            An atom collection object representing either the substrate or an
            adsorbed molecule. The default is None.

        """
        if isinstance(obj, Slab):
            self._a0 = obj._a0
            self._m = obj._m
            self._n = obj._n
            self._maxm = obj._maxm
            self._maxn = obj._maxn

        valid_z = array([atom._x[2] for atom in self._atoms if atom._active])
        self._top = max(valid_z, axis=0)
        self._bottom = min(valid_z, axis=0)
        self._latvec[2] = array(
            [0.0, 0.0, (self._top - self._bottom + self._vaccuum) / self._a0]
        )

    def __str__(self):
        """
        Returns the ID of the MoleculeAdsorption object.

        Returns
        -------
        String
            A string containing the ID of the MoleculeAdsorption object.
        """
        return "<MoleculeAdsorption object> ID: %d" % self._id

    """
    Properties
    ----------
    substrate : Slab, readonly
        Slab object that represents the substrate where molecules are adsorbed.
    adsorbed_molecules : Python tuple, readonly
        Tuple of adsorbed molecule(s).
    minimum_separation : float
        Minimum separation between the atom(s) in the molecule(s) that are the 
        closest to the substrate.
    vaccuum : float
        Empty space separating the highest and lowest Z coordinates in the adsorbed, 
        system considering periodic boundary conditions along the Z direction.
    top : float, readonly
        Maximum Z atomic coordinate.
    bottom : float, readonly
        Minimum Z atomic coordinate.
    a0 : float
        Lattice constant.
    lattice vectors : Numpy array, readonly
        Lattice vectors.
    origin : Numpy array
        Position with minimum values of X, Y, and Z coordinates of the substrate.
    """

    @property
    def substrate(self):
        return self._components["substrate"]

    @property
    def adsorbed_molecules(self):
        return tuple(
            self._components["molecules@top"] + self._components["molecules@bottom"]
        )

    @property
    def adsorbed_polymers(self):
        return tuple(
            self._components["polymer@top"] + self._components["polymer@bottom"]
        )

    @property
    def minimum_separation(self):
        return self._minsep

    @minimum_separation.setter
    def minimum_separation(self, val):
        if isinstance(val, (float, int)) and val > 0.0:
            self._minsep = val

            for molecule in self.adsorbed_molecules:
                dz = 0.0

                if (
                    molecule in self._components["molecules@top"]
                    and (molecule._minz - self._components["substrate"]._top)
                    < self._minsep
                ):
                    dz = (
                        self._components["substrate"]._top
                        + self._minsep
                        - molecule._minz
                    )
                elif (
                    molecule in self._components["molecules@bottom"]
                    and (self._components["substrate"]._bottom - molecule._maxz)
                    < self._minsep
                ):
                    dz = (
                        self._components["substrate"]._bottom
                        - self._minsep
                        - molecule._maxz
                    )

                if abs(dz) > 0.0:
                    molecule.displace(
                        array([0.0, 0.0, dz]), call_handler=False, call_update=False
                    )

            self._update()
        else:
            raise MoleculeAdsorptionError(
                "Minimum separation must be a number greater than zero!"
            )

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
            raise MoleculeAdsorptionError(
                "Vaccuum must be provided as a number greater than zero!"
            )

    @property
    def top(self):
        return self._top

    @property
    def bottom(self):
        return self._bottom

    @property
    def a0(self):
        return self.components["substrate"].a0

    @a0.setter
    def a0(self, val):
        self._components["substrate"].a0 = val

    @property
    def lattice_vectors(self):
        return self._components["substrate"]._latvec

    @property
    def origin(self):
        return self._components["substrate"]._origin

    @origin.setter
    def origin(self, val):
        if isinstance(val, (list, tuple)):
            val = array(val)

        if isinstance(val, ndarray) and val.shape[0] == 3:
            disp = val - self._components["substrate"]._origin

            self._components["substrate"].displace(disp, call_handler=False)

            for molecule in self.adsorbed_molecules:
                molecule.displace(disp, call_handler=False, call_update=False)

            self._update()
        else:
            raise MoleculeAdsorptionError(
                "The new origin must be a Numpy array with three components!"
            )
