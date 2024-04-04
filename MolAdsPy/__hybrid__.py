from __future__ import print_function
from __atomcollection__ import AtomCollection
from __exception__ import BasicException
from abc import ABC,abstractmethod
from numpy import array

class HybridError(BasicException):
    pass

class Hybrid(AtomCollection,ABC):
    __slots__=["_components","_a0","_latvec","_n","_m","_l","_maxn","_maxm",
               "_maxl","_origin"]
    
    def __init__(self,**kwargs):
        '''
        Object initialization.
        '''
        self._components={}                     # Dictionary of objects that are put together to form the hybrid system
        self._a0=1.0                            # Lattice constant
        self._latvec=array([[1.0,0.0,0.0],
                            [0.0,1.0,0.0],
                            [0.0,0.0,1.0]])     # Lattice vectors
        self._n=self._m=self._l=0               # Number of translations along the first, second and third lattice vectors
        self._maxn=self._maxm=self._maxl=0      # Maximum number of translations for resizing purposes
        self._origin=array([0.0,0.0,0.0])       # Position with minimum values of X, Y, and Z coordinates of the structure.
        
        self._verbose = kwargs.get("verbose",False)
        type(self)._append(self)
        
    @abstractmethod
    def add_component(self,*args):
        '''
        To be implemented in a derived hybrid structure, this method is expected 
        to add atom collection objects to the hybrid system.
        '''        
        pass
    
    @abstractmethod
    def remove_component(self,*args):
        '''
        To be implemented in a derived hybrid structure, this method is expected 
        to remove an atom collection object from the hybrid system.
        '''        
        pass    

    def force_update(self):
        '''
        Forces the hybrid system's attributes to be updated if there is a change 
        in the attributes of any of its components.
        '''
        self._update()
        
    def write_xyz(self,file_name="coords.xyz"):
        '''
        Saves the atomic coordinates of the hybrid system into an XYZ file.

        Parameters
        ----------
        file_name : string, optional
            Name of the XYZ file. The default is "coords.xyz".
        '''
        super().write_xyz(file_name,ucell=False)        
        
    def write_pw_input(self,file_name="pw.in",pseudopotentials={},pwargs={}):
        '''
        Creates a basic input file for geometry relaxation of the hybrid system 
        using the pw.x code found in the Quantum Espresso package.

        Parameters
        ----------
        file_name : string, optional
            Name of the input file. The default is "pw.in".
        pseudopotentials : Python dictionary, optional
            Specify for each element provided as a key the name of a pseudopotential 
            file given as the corresponding value. The default is {}.
        pwargs : Python dictionary
            Dictionary containing key-value pairs that allow some customization 
            of the input file. At the moment, those are the keys accepted:
            Dictionary containing key-value pairs that allow some customization 
            of the input file. At the moment, those are the keys accepted:
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
        '''
        ucell=False
        
        super().write_pw_input(file_name,ucell,pseudopotentials,pwargs)        
    
    def copy(self):
        '''
        This method is disabled for hybrid objects and issues an error message.
        '''
        raise HybridError("Not implemented for hybrid objects!")
        
    def add_atom(self,*args):
        '''
        This method is disabled for hybrid objects and issues an error message.
        '''
        raise HybridError("Not implemented for hybrid objects!")
        
    def remove_atom(self,*args):
        '''
        This method is disabled for hybrid objects and issues an error message.
        '''
        raise HybridError("Not implemented for hybrid objects!")
        
    def location(self,*args):
        '''
        This method is disabled for hybrid objects and issues an error message.
        '''
        raise HybridError("Not implemented for hybrid objects!")
           
    def _read_xyz(self,*args):
        '''
        This method is disabled for hybrid objects and issues an error message.
        '''
        raise HybridError("Not implemented for hybrid objects!")
        
    def __setitem__(self,*args):
        '''
        This method is disabled for hybrid objects and issues an error message.
        '''
        raise HybridError("Not implemented for hybrid objects!")

    '''
    Properties
    ----------
    _atoms : Python tuple, readonly
        Tuple of atoms in the Hybrid object.
    _species : Python tuple, readonly
        Tuple of atomic species in the Hybrid object.
    _loc : Python dictionary, readonly
        Dictionary of three indices, indicating the position of every atom in 
        the hybrid system supercell.
    '''
    @property
    def _atoms(self):
        atom_list=[]
        
        for key in self._components:
            if isinstance(self._components[key],(list,tuple)):
                for component in self._components[key]:
                    atom_list+=component._atoms
            else:
                atom_list+=self._components[key]._atoms
            
        return tuple(atom_list)
        
    @property
    def _species(self):
        spec_list=[]
        
        for key in self._components:
            if isinstance(self._components[key],(list,tuple)):
                for component in self._components[key]:
                    spec_list+=component._species
            else:
                spec_list+=self._components[key]._species
            
        return tuple(set(spec_list))
        
    @property
    def _loc(self):
        loc_dict={}
        
        for key in self._components:
            if isinstance(self._components[key],(list,tuple)):
                for component in self._components[key]:
                    loc_dict.update(component._loc)
            else:
                loc_dict.update(self._components[key]._loc)        
            
        return loc_dict