# MolAdsPy

This is a package written entirely in Python to place and manipulate molecules on a substrate and then generate coordinates that can be used as inputs to first-principles calculation codes like [Quantum Espresso](https://www.quantum-espresso.org/) or [VASP](https://www.vasp.at/).

In more specific terms, the package allows the user to create objects such as molecules and slabs, which are made of atomic objects, and combine them into a system in which one or more molecular objects are adsorbed onto a slab object that acts as a substrate. An adsorbed molecule can be moved to a specific location on the substrate (e.g., a predefined adsorption site), can be shifted a certain distance from its current position, and can be rotated in three dimensions. Last but not least, an adsorbed molecule can be desorbed (removed) from the substrate.

## Object model

The package consists of five classes:

* **Species** defines the properties of the chemical elements of atoms.
* **Atom** defines the attributes of atoms and the methods for manipulating this type of object.
* **Molecule** defines the attributes of molecules and the methods for manipulating this type of object.
* **Slab** defines the attributes of substrates and the methods for handling this type of object.
* **Adsorption** defines the attributes of a system that contains a substrate (Slab object) and one or more molecules (Molecule object(s)) adsorbed on it and the methods for manipulating these molecules.

## Basic usage

The main use of the package is to generate coordinates for first principle calculations of systems in which one or more molecules are adsorbed onto a substrate.

There are different ways to do this and the Jupyter notebooks in the *tests* directory offer examples of usage. In general, the process starts by instantiating a Slab object. This object will play the role of substrate where the molecules are adsorbed. Then one or more Molecule objects must be created. Substrate and molecule(s) are brought together by means of an Adsorption object.

### Creation of an Atom object

Molecule and Slab objects are both made of Atom objects. To create an Atom object, first a Species object must be instantiated:

```
Species(symbol,element,atomicnumber,atomicmass,pseudopotential)
        Parameters
        ----------
        symbol : string
            Symbol of the chemical element, e.g., "C" for carbon.
        element : string, optional
            Name of the atomic element, e.g., "carbon". The default is None.
        atomicnumber : integer, optional
            Atomic number of the element. The default is None.
        atomicmass : float, optional
            Atomic mass of the element, usually in atomic units. The default is 
            None.
        pseudopotential : string, optional
            Name of the file containing the pseudopotential for the element,
            to be employed in DFT calculations. The default is None.
```

Then an Atom object can be instantiated:

```
Atom(species,x,fixed)
        Parameters
        ----------
        species : Species object
            Atomic species of the atom.
        x : Python list, optional
            Cartesian coordinates of the atom. The default is [0.0,0.0,0.0].
        fixed : Python list, optional
            Defines if a component of the atomic coordinates remains unaltered 
            in the course of a geometry optimization. The default is 
            [False,False,False].
```

> Most of the time, Molecule and Slab objects will be created by reading atomic coordinates from files, such as XYZ format files, so Species and Atom objects will not need to be created explicitly.

### Creation of a Slab object

The simplest (but not the only) way to instantiate a Slab object is to get the atomic positions and lattice vectors from an XYZ file:

```
Slab(slabname,filename,filetype)
        Parameters
        ----------
        slabname : string
            A label allowing you to identify the slab, for instance, "graphene".
        filename : string
            Name of the file containing the slab structure.
        filetype : string
            Type of the file containing the slab structure. For now, only 'XYZ' 
            is accepted.
```

It is usually a good idea to define adsorption sites on the Slab object where the molecule(s) will be adsorbed. The adsorption sites are grouped in a Python dictionary where each key identifies a site type. To add an adsorption site to the Slab object, the **add_adsorption_site** method should be used:

```
add_adsorption_site(ads_site,pos)
        Parameters
        ----------
        ads_site : string
            Label of the adsorption site.
        pos : Python list
            Coordinates of the adsorption site on the XY plane.
```

Sometimes, it is convenient to increase the size of the slab using the **replicate** method:

```
replicate(n,m,replicate_ads_sites)
        Parameters
        ----------
        n : integer
            Number of replications along the first lattice vector.
        m : integer
            Number of replications along the second lattice vector.
        replicate_ads_sites : logical, optional
            Whether also replicate the adsorption sites or not. The default is 
            True.
```

The atomic coordinates and lattice vectors of the Slab object can be saved in a file in XYZ format by calling the **write_xyz** method:

```
write_xyz(filename)
        Parameters
        ----------
        filename : string, optional
            Name of the XYZ file. The default is "slab.xyz".
```

Once created, there are several attributes of the Slab object that can be accessed by the user:

```
    label : string, readonly
        Label that identifies the slab.
    atoms : Python list, readonly
        List of atoms in the slab.
    atomtypes : Python list, readonly
        List of atom types in the slab.
    adsorptionsites : Python dictionary, readonly
        Dictionary of adsorption sites on the slab surface.
    top : float, readonly
        Top of the substrate surface.
    bottom : float, readonly
        Bottom of the substrate surface.
    origin : Python list
        Position with minimum values of X, Y, and Z coordinates of the slab.
    latticevectors : Python list
        Repetition vectors in the XY plane.
    ID : integer, readonly
        Unique slab identifier.
```

> Information about other Slab object methods and how to use them can be found in the comments in the *\__slab.py\__* module.

### Creation of a Molecule object

A Molecule object, particularly if it represents a molecule made up of more than a handful of atoms, can be created from a file in XYZ format:

```
Molecule(moltype,filename,filetype)
        Parameters
        ----------
        moltype : string
            A label allowing you to identify the type of the molecule, for 
            instance, "benzene" or "C6H12".
        filename : string
            Name of the file containing the molecule structure.
        filetype : string
            Type of the file containing the molecule structure. For now, only 
            'XYZ' is accepted.
```

Anchors are reference points on the Molecule object that are used when the molecule is moved or rotated. Every Molecule object has at least one default anchor point, whose label is "com", corresponding to the center of mass of the molecule. Anchors points are grouped in a Python dictionary and can be created using the **add_anchor** method:

```
add_anchor(anchor,pos)
        Parameters
        ----------
        anchor : string
            Name of the anchor point.
        pos : Python list
            Cartesian coordinates of the anchor point.
```

The atomic coordinates of the Molecule object can be saved in a file in XYZ format by calling the **write_xyz** method:

```
write_xyz(filename)
        Parameters
        ----------
        filename : string, optional
            Name of the XYZ file. The default is "molecule.xyz".
```

Once created, there is a number of attributes of the Molecule object that can be accessed by the user:

    moleculetype : string, readonly
        Label that identifies the type of the molecule.
    atoms : Python list, readonly
        List of atoms in the molecule.
    anchors : Python dictionary, readonly
        Dictionary of anchor points.
    atomtypes : Python list, readonly
        List of atomic species in the molecule.
    centerofmass : Python list
        Center of mass of the molecule.
    maxx,minx,maxy,miny,maxz,minz : floats, readonly
        Molecule extremities.
    ID : int, readonly
        Unique molecule identifier.

Additionally, a Molecule object can be created as a copy of another existing Molecule object using the **copy** method:

```
copy()
        Returns
        -------
        Molecule object
            Copy of the molecule.
```

> Information about other Molecule object methods and how to use them can be found in the comments in the *\__molecule.py\__* module.

### Adsorbing molecule(s) onto a substrate

The adsorption of one or more molecules onto a substrate is done by creating an Adsorption object taking as its argument a Slab object that will serve as the substrate:

```
Adsorption(substrate)   
        Parameters
        ----------
        substrate : Slab object
            Substrate where one or more molecules can be adsorbed.
```

After the Adsorption object is created, a molecule is adsorbed using the **add_molecule** method which has two versions:

```
add_molecule(molecule,ads_site,ads_site_index,anchor,vert_sep)
        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        ads_site : string
            Label of the adsorption site.
        ads_site_index : int
            Index of an adsorption site on the substrate.
        anchor : string
            Anchor point in the molecule that will be vertically aligned with 
            the adsorption site on the substrate.
        vert_sep : float
            Nearest vertical distance separating the molecule from the 
            substrate.
            
 add_molecule(molecule,pos,anchor,vert_sep)
        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        pos : Python list
            XY coordinates in the substrate above which the molecule will be 
            placed.
        anchor : string
            Anchor point in the molecule that will be vertically aligned with 
            the adsorption site on the substrate.
        vert_sep : float
            Nearest vertical distance separating the molecule from the 
            substrate.
```

An adsorbed molecule can be removed using one of two versions of the **remove_molecule** method:

```
remove_molecule(molecule)
        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule to be removed.
            
remove_molecule(molid) -> removes an adsorbed molecule.
        Parameters
        ----------
        molid : int
            ID of the adsorbed molecule.
```

A molecule can be moved to a specific location on the substrate using one of the four versions of the **move_molecule__to** method:

```
move_molecule_to(molid,x,anchor)
        Parameters
        ----------
        molid : int
            ID of the adsorbed molecule.
        x : Python list
            New position of the molecule's anchor point.
        anchor : string
            Anchor point used as reference for translations.
            
move_molecule_to(molecule,x,anchor)
        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        x : Python list
            New position of the molecule's anchor point.
        anchor : string
            Anchor point used as reference for translations.
            
move_molecule_to(molid,ads_site,ads_site_index,anchor,vert_sep)
        Parameters
        ----------
        molid : int
            ID of the adsorbed molecule.
        ads_site : string
            Label of the adsorption site.
        ads_site_index : int
            Index of an adsorption site on the substrate.
        anchor : string
            Anchor point in the molecule that will be vertically aligned with 
            the adsorption site on the substrate.
        vert_sep : float
            Nearest vertical distance separating the molecule from the 
            substrate.
            
move_molecule_to(molecule,ads_site,ads_site_index,anchor,vert_sep)
        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        ads_site : string
            Label of the adsorption site.
        ads_site_index : int
            Index of an adsorption site on the substrate.
        anchor : string
            Anchor point in the molecule that will be vertically aligned with 
            the adsorption site on the substrate.
        vert_sep : float
            Nearest vertical distance separating the molecule from the 
            substrate.
```

A molecule can be rotated around a previously defined anchor point using one of two versions of the **rotate_molecule** method:

```
rotate_molecule(molid,theta,phi,psi,anchor)
        Parameters
        ----------
        molid : int
            ID of the molecule that will be rotated around the anchor point.
        theta : float
            First Euler angle, in degrees.
        phi : float
            Second Euler angle, in degrees.
        psi : float
            Third Euler angle, in degrees.
        anchor : string
            Anchor point around which the molecule will be rotated.

rotate_molecule(molecule,theta,phi,psi,anchor)
        Parameters
        ----------
        molecule : Molecule object
            Molecule that will be rotated around the anchor point.
        theta : float
            First Euler angle, in degrees.
        phi : float
            Second Euler angle, in degrees.
        psi : float
            Third Euler angle, in degrees.
        anchor : string
            Anchor point around which the molecule will be rotated.
```

The separation between a molecule and the substrate (i.e., the difference between the lowest z-coordinate of the molecule and the top of the substrate) can be set by calling one of the versions of the **set_separation** method:

```
set_separation(molid,sep)
        Parameters
        ----------
        molid : int
            ID of the adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.
            
set_separation(molecule,sep)
        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.
```

> It is worth noting that when using the methods of the Adsorption object to manipulate a molecule adsorbed onto the substrate, a user-defined minimum separation distance is enforced in order to prevent the molecule and substrate from overlapping.

The atomic coordinates of the adsorbed system can be saved in a file in XYZ format by calling the **write_xyz** method:

```
write_xyz(filename,vacuum)
        Parameters
        ----------
        filename : string, optional
            Name of the XYZ file. The default is "adsorbed.xyz".
        vacuum : float, optional
            Height of the vacuum layer to be added to the top of the adsorbed 
            system. The default is 10.0 (Angstroms).
```

Once created, there is a number of attributes of the Adsorption object that can be accessed by the user:

```
substrate : Slab, readonly
    Slab object that represents the substrate where molecules are adsorbed.
adsorbedmolecules : Python list, readonly
    List of adsorbed molecule(s).
atomtypes : Python list, readonly
    List of atom types in the slab+molecule(s) system.
minimumseparation : float
    Minimum separation between the atom(s) with lowest Z value in the 
    adsorbed molecule(s) and the atom(s) with highest Z value in the 
    substrate. The default is 1.0 Angstroms.
```

> Information about other Adsorption object methods and how to use them can be found in the comments in the *\__adsorption.py\__* module.

## Installation

To install the latest version of the package, the *pip* command can be used as follows:

```
pip install MolAdsPy
```

The *multipledispatch* package is a requirement and will also be installed if not found on the machine.

## Citation

If you use *MolAdsPy* in your own projects and if you don't mind, when publishing the results, please cite the following papers of mine:

* [*Structural, electronic, and magnetic properties of pristine and oxygen-adsorbed graphene nanoribbons*](https://www.sciencedirect.com/science/article/abs/pii/S0169433210003855)
* [*Interfacial properties of polyethylene/Ti3C2Tx mxene nanocomposites investigated by first-principles calculations*](https://www.sciencedirect.com/science/article/abs/pii/S0169433222028720)
* [*Self-assembly of NiTPP on Cu(111): a transition from disordered 1D wires to 2D chiral domains*](https://pubs.rsc.org/en/content/articlelanding/2015/cp/c5cp01288k)

## Final disclaimer

The package documentation is at a very early stage. It will be updated over time, but the user can always refer to comments on the methods and properties of the modules that make up the package for help on how to use them.

This code is provided as is. The author makes no guarantee that its results are accurate and is not responsible for any losses caused by the use of the code.

---

If you have any comments, suggestions or corrections, please feel free to [drop me a message](mailto:roberto.veiga@ufabc.edu.br).
