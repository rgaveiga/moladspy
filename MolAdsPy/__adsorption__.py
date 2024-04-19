from __future__ import print_function
from __molecule__ import Molecule
from __slab__ import Slab
from __hybrid__ import Hybrid
from __polymer__ import Polymer
from __exception__ import BasicException
from numpy import array,dot,max,min,ndarray
from multipledispatch import dispatch
from inspect import signature

class AdsorptionError(BasicException):
    pass

class Adsorption(Hybrid):
    __slots__=["_substrate","_minsep","_vaccuum","_side","_top","_bottom"]
    __MAXREPMISMATCH__ = 1000
    __MISMATCHTOLSTD__ = 5.0
    def __init__(self,substrate,minsep=1.0,vaccuum=10.0,**kwargs):        
        r'''
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
        
        '''
        super().__init__(**kwargs)
        
        if isinstance(substrate,Slab):
            if substrate._belongs_to is not None:
                raise AdsorptionError("The substrate already belongs to another structure!")
            
            self._components["substrate"]=substrate
            self._origin=self._components["substrate"]._origin
            self._latvec=self._components["substrate"]._latvec
            self._components["molecules@top"]=[]
            self._components["molecules@bottom"]=[]  
            self._components["polymers@top"]=[]
            self._components["polymers@bottom"]=[] 
                
            if isinstance(minsep,(int,float)) and minsep>0.0:
                self._minsep=minsep
            else:
                raise AdsorptionError("'minsep' must be a number greater than zero!")

            if isinstance(vaccuum,(int,float)) and vaccuum>0.0:
                self._vaccuum=vaccuum
            else:
                raise AdsorptionError("'vaccuum' must be a number greater than zero!")
        else:
            raise AdsorptionError("Substrate must be a Slab object!")
            
        substrate._belongs_to=self
            
        self._update()
    
    #region METHODS ADD AND REMOVE IMPLEMENTATIONS
    @dispatch(Molecule,str,n=int,m=int,anchor=str,vertsep=(int,float),side=str)    
    def add_component(self,molecule,ads_site,n=0,m=0,anchor="com",vertsep=1.0,
                      side="top"):
        '''
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
        '''
        if molecule._belongs_to is not None:
            raise AdsorptionError("The molecule already belongs to another structure!")
        elif self._components["substrate"] is None:
            raise AdsorptionError("No substrate has been defined!")
        elif len(anchor)==0 or anchor not in molecule.anchors.keys():
            raise AdsorptionError("'anchor' must be a valid anchor point!")
        elif n<0:
            raise AdsorptionError("'n' must be an integer greater than or equal to zero!")
        elif m<0:
            raise AdsorptionError("'m' must be an integer greater than or equal to zero!")
        
        if len(ads_site)>0 and ads_site in self._components["substrate"]._ads_sites:
            sep=max([vertsep,self._minsep])
            
            if side=="top":
                dz=self._components["substrate"]._top+sep-molecule._minz
            elif side=="bottom":
                dz=self._components["substrate"]._bottom-sep-molecule._maxz
            else:
                raise AdsorptionError("'side' must be either 'top' or 'bottom'!")
                
            self._components["molecules@"+side].append(molecule)

            molecule._belongs_to=self
            x0=self._components["substrate"]._ads_sites[ads_site]+\
                dot(self._a0*array([n,m]),self._latvec[:2,:2])
            x=array([x0[0],x0[1],molecule.anchors[anchor][2]+dz])
            
            molecule.move_to(x,anchor)            
            self._update()
        else:
            raise AdsorptionError("A valid adsorption site must be provided!")
    
    # TODO: Verify if Align if not "Z". Throw error if it is.
    # TODO: Allow just one polymer. Make a separate dictionary.
    # TODO: Verify _lacvec same for slab and polymer because of mismatch to make sense. x and y by default.
    @dispatch(Polymer,anchor=str,vertsep=(int,float),side=str)    
    def add_component(self,polymer,anchor="com",vertsep=1.0, side="top", **kwargs):
        '''
        Adsorbs a polymer molecule onto the substrate.

        Parameters
        ----------
        polymer : Polymer molecule object
            Adsorbed polymer molecule.
        vertsep : float, optional
            Nearest vertical distance separating the polymer molecule from the 
            substrate. The default is 1.0 Angstrom.
        side : string, optional
            Side of the slab where the molecule will be adsorbed. It must be either 
            "top" or "bottom" (self-explanatory). The default is "top".
        **kwargs : See below

        Key arguments
        -----------

        skip_mismatch_adjust : boolean
            False by standard. Just used when the user wants to skip this adjustment inside this function.
        
        '''

        # Get key arguments
        skip_mismatch_adjust = kwargs.get("skip_mismatch_adjust",False)

        if polymer._orientation == 2:
            raise AssertionError("The polymer molecule is perpendicular to slab!")
        if polymer._belongs_to is not None:
            raise AdsorptionError("The polymer already belongs to another structure!")
        elif self._components["substrate"] is None:
            raise AdsorptionError("No substrate has been defined!")
        
        sep=max([vertsep,self._minsep])
        
        if side=="top":
            dz=self._components["substrate"]._top+sep-polymer._minz
        elif side=="bottom":
            dz=self._components["substrate"]._bottom-sep-polymer._maxz
        else:
            raise AdsorptionError("'side' must be either 'top' or 'bottom'!")
        
        # This object just allow one polymer because of mismatch problems with the slab.
        if(len(self._components["polymers@top"]) > 0 or len(self._components["polymers@bottom"]) > 0):
            raise AdsorptionError("The adsorbed already has a polymer!")
        
        self._components["polymers@"+side].append(polymer)

        self._adjustMismatchPolymerSlab(polymer)
        
        polymer._belongs_to=self
        #x0=self._components["substrate"]._ads_sites[ads_site]+\
        #    dot(self._a0*array([n,m]),self._latvec[:2,:2])
        #x=array([x0[0],x0[1],polymer.anchors[anchor][2]+dz])
        
        #polymer.move_to(x,anchor)            
        self._update()

    @dispatch(Molecule,(ndarray,list,tuple),anchor=str,vertsep=(int,float),side=str)
    def add_component(self,molecule,pos,anchor="com",vertsep=1.0,side="top"):
        '''
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
        '''
        if molecule._belongs_to is not None:
            raise AdsorptionError("The polymer already belongs to another structure!")
        elif self._components["substrate"] is None:
            raise AdsorptionError("No substrate has been defined!")
        elif len(anchor)==0 or anchor not in molecule.anchors.keys():
            raise AdsorptionError("'anchor' must be a valid anchor point!")
            
        if isinstance(pos,(list,tuple)):
            pos=array(pos)
        
        if isinstance(pos,ndarray) and pos.shape[0]==2:
            sep=max([vertsep,self._minsep])
            
            if side=="top":
                dz=self._components["substrate"]._top+sep-molecule._minz
            elif side=="bottom":
                dz=self._components["substrate"]._bottom-sep-molecule._maxz
            else:
                raise AdsorptionError("'side' must be either 'top' or 'bottom'!")
            
            x=array([pos[0],pos[1],molecule.anchors[anchor][2]+dz])
                        
            self._components["molecules@"+side].append(molecule)
            
            molecule._belongs_to=self
            
            molecule.move_to(x,anchor)
            self._update()
        else:
            raise AdsorptionError("'pos' must be provided as a list with two components!")
    
    # TODO: Implement
    @dispatch(Polymer,(ndarray,list,tuple),anchor=str,vertsep=(int,float),side=str)
    def add_component(self,polymer,pos,anchor="com",vertsep=1.0,side="top"):
        # TODO: Allow just one polymer because of mismatch problem
        '''
        Adsorbs a Polymer molecule onto the substrate.

        Parameters
        ----------
        polymer : Polymer object
            Adsorbed polymer.
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
        '''
        if polymer._belongs_to is not None:
            raise AdsorptionError("The polymer already belongs to another structure!")
        elif self._components["substrate"] is None:
            raise AdsorptionError("No substrate has been defined!")
        elif len(anchor)==0 or anchor not in polymer.anchors.keys():
            raise AdsorptionError("'anchor' must be a valid anchor point!")
            
        if isinstance(pos,(list,tuple)):
            pos=array(pos)
        
        if isinstance(pos,ndarray) and pos.shape[0]==2:
            sep=max([vertsep,self._minsep])
            
            if side=="top":
                dz=self._components["substrate"]._top+sep-polymer._minz
            elif side=="bottom":
                dz=self._components["substrate"]._bottom-sep-polymer._maxz
            else:
                raise AdsorptionError("'side' must be either 'top' or 'bottom'!")
            
            x=array([pos[0],pos[1],polymer.anchors[anchor][2]+dz])
                        
            self._components["molecules@"+side].append(polymer)
            
            polymer._belongs_to=self
            
            polymer.move_to(x,anchor)
            self._update()
        else:
            raise AdsorptionError("'pos' must be provided as a list with two components!")
    
    @dispatch(Molecule)    
    def remove_component(self,molecule):
        '''
        Removes an adsorbed molecule.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule to be removed.
        '''
        if molecule in self._components["molecules@top"]:
            self._components["molecules@top"].remove(molecule)
        elif molecule in self._components["molecules@bottom"]:
            self._components["molecules@bottom"].remove(molecule)
        else:
            raise AdsorptionError("Molecule not found!")
            
        molecule._belongs_to=None
        
        self._update()
    
    # TODO: Implement
    @dispatch(Polymer)    
    def remove_component(self,polymer):
        '''
        Removes an adsorbed polymer.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule to be removed.
        '''
        if polymer in self._components["molecules@top"]:
            self._components["molecules@top"].remove(polymer)
        elif polymer in self._components["molecules@bottom"]:
            self._components["molecules@bottom"].remove(polymer)
        else:
            raise AdsorptionError("Polymer not found!")
            
        polymer._belongs_to=None
        
        self._update()
    
    @dispatch(int)
    def remove_component(self,molid):
        '''
        Removes an adsorbed molecule.    

        Parameters
        ----------
        molid : integer
            ID of the adsorbed molecule.
        '''
        for molecule in self.adsorbed_molecules:
            if molecule.ID==molid:
                if molecule in self._components["molecules@top"]:
                    self._components["molecules@top"].remove(molecule)
                elif molecule in self._components["molecules@bottom"]:
                    self._components["molecules@bottom"].remove(molecule)
                    
                break
        else:
            raise AdsorptionError("Molecule ID not found!")
            
        molecule._belongs_to=None
        
        self._update()

    #endregion
    
    #region AUXILIARY FUNCTIONS FOR MISMATCH ADJUSTMENT
    def _adjustMismatchPolymerSlab(self,polymer):
        r'''
        Adjust Resize of Slab and Polymer to adjust mismatch in defined tolerance   

        Parameters
        ----------
        polymer : Polymer object
            Object of the polymer that will be adjusted over the slab
        '''
        slab = self._components["substrate"]

        # Reset resize of Polymer and Slab on the polymer orientation
        
        resize_vector = [slab._n,slab._m,slab._l]
        resize_vector[polymer._orientation] = 0
        slab.resize(resize_vector[0],resize_vector[1])

        # Resize necessary to minimize mismatch between polymer and slab
        rezize_slab,resize_polymer = self._polymer_slab_match_box_length(polymer)
        resize_vector = [slab._n,slab._m,slab._l]
        resize_vector[polymer._orientation] = rezize_slab 

        print()
        print("Previous Resize Slab: " + str(slab._n) + ", " + str(slab._m) + ", " + str(slab._l))
        print("Previous Resize Polymer: " + str(polymer._n) + ", " + str(polymer._m) + ", " + str(polymer._l))
        print("Slab will be resized by : " + str(resize_vector))
        print("Polymer will be resized: " + str(resize_polymer))
        # Resizing objects
        slab.resize(resize_vector[0],resize_vector[1])
        polymer.resize(resize_polymer)
    
    def _polymer_slab_match_box_length(self,polymer,mismatch_tol_percent=__MISMATCHTOLSTD__,maxrep=__MAXREPMISMATCH__):
        '''
        Returns the resize number i,j respective to Slab and Polymer necessary for mismatch to be under given tolerance (%)

        Parameters
        ----------
        n : integer
            Number of repetitions of the substrate's unit cell along the first lattice 
            vector.
        m : integer
            Number of repetitions of the substrate's unit cell along the second lattice 
            vector.
        
        Output
        -----------
        i : Resize number for Slab
        j : Resize number for Polymer
        '''        
        mismatch_tol = mismatch_tol_percent/100.0 # Parameter is given in percentage

        pol_ori = polymer._orientation
        slab = self._components["substrate"]
        
        # This exists to get original size of cells independent of previous resizing
        polymer_resize_vector = [polymer._n,polymer._m,polymer._l]
        polymer_original_resize = polymer_resize_vector[pol_ori]
        slab_original_resize_vector = [slab._n,slab._m,slab._l]
        slab_original_resize = slab_original_resize_vector[pol_ori]

        a0_1 = slab._a0*slab._latvec[pol_ori][pol_ori]
        a0_2 = polymer._a0*polymer._latvec[pol_ori][pol_ori]

        print()
        print("Polymer Original Resize: " + str(polymer_original_resize))
        print("Slab Original Resize: " + str(slab_original_resize))
        print("a0 polymer: " + str(a0_2))
        print("a0 slab: " + str(a0_1))
        print("latvec polymer: " + str(polymer._latvec))
        print("latvec slab: " + str(slab._latvec))
        
        i = 1
        while(i < maxrep):
            j = int(a0_1*i/a0_2)
            mismatch = abs(j*a0_2 - i*a0_1)/(i*a0_1)

            if(mismatch <= mismatch_tol):
                return int(i-1),int(j-1)
            i+=1
    #endregion

    #region METHODS SET_SEPARATION
    @dispatch(Molecule,(int,float))
    def set_separation(self,molecule,sep,**kwargs):
        '''
        Rigidly displaces the molecule along the direction perpendicular to the 
        substrate, such that the separation between the top/bottom of the substrate 
        and the nearest atom in the molecule be equal to 'sep'.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.
        **kwargs : See below
            
        Keywork Arguments
        ----------
        sep_violation : boolean
            The user gets the freedom to violate the set minsep using this flag

        '''

        sep_violation = kwargs.get('sep_violation', False)

        if(not sep_violation):
            if sep<self._minsep:
                raise AdsorptionError("Molecule will be too close to the substrate!")
        
        dz=0.0    
        
        if molecule in self._components["molecules@top"]:
            dz=self._components["substrate"]._top+sep-molecule._minz
        elif molecule in self._components["molecules@bottom"]:
            dz=self._components["substrate"]._bottom-sep-molecule._maxz
        else:
            raise AdsorptionError("Molecule not found!")
            
        if abs(dz)>0.0:
            molecule.displace(array([0.0,0.0,dz]))
            self._update()
    
    # TODO: Avoid recursion. Make flag to update function not to call set_separation again
    @dispatch(Polymer,(int,float))
    def set_separation(self,polymer,sep,update=True):
        '''
        Rigidly displaces the polymer along the direction perpendicular to the 
        substrate, such that the separation between the top/bottom of the substrate 
        and the nearest atom in the molecule be equal to 'sep'.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.
        **kwargs : See below

        Keywork Arguments
        ----------
        sep_violation : boolean
            The user gets the freedom to violate the set minsep using this flag

        '''
        
        if sep<self._minsep:
            raise AdsorptionError("Polymer will be too close to the substrate!")
        
        dz=0.0

        if polymer in self._components["molecules@top"]:
            dz=self._components["substrate"]._top+sep-polymer._minz
        elif polymer in self._components["molecules@bottom"]:
            dz=self._components["substrate"]._bottom-sep-polymer._maxz
        else:
            raise AdsorptionError("Molecule not found!")
            
        if abs(dz)>0.0:
            polymer.displace(array([0.0,0.0,dz]))
            if(update):
                self._update()
    
    @dispatch(int,(int,float))
    def set_separation(self,molid,sep,**kwargs):
        '''
        Rigidly displaces the molecule along the direction perpendicular to the 
        substrate, such that the separation between the top/bottom of the substrate 
        and the nearest atom in the molecule be equal to 'sep'.

        Parameters
        ----------
        molid : integer
            ID of the adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.
        '''
        sep_violation = kwargs.get('sep_violation', False)

        if(not sep_violation):
            if sep<self._minsep:
                raise AdsorptionError("Molecule will be too close to the substrate!")
        
        for molecule in self.adsorbed_molecules:
            if molecule.ID==molid:
                dz=0.0
                
                if molecule in self._components["molecules@top"]:
                    dz=self._components["substrate"]._top+sep-molecule._minz
                elif molecule in self._components["molecules@bottom"]:
                    dz=self._components["substrate"]._bottom-sep-molecule._maxz
                else:
                    raise AdsorptionError("Molecule not found!")
                    
                if abs(dz)>0.0:
                    molecule.displace(array([0.0,0.0,dz]))
                    self._update()
                    
                break
        else:
            raise AdsorptionError("Molecule ID not found!")
    #endregion


    #region Begin and End Hanlders for children methods

    def beginHandler(self,obj,method,method_locals,**kwargs):
        '''
        This Handler is used on the begining of AtomCollection objects
        to communicate to the Hybrid arguments and the method and check if the operation
        if valis when the object is bounded by this Hybrid

        Parameters
        ----------
        obj : AtomCollection
            ID of the adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.
        '''
        if isinstance(obj,Polymer):
            if(method.__name__ == 'align'):
                if(method_locals['axis'] == 'z'):
                    raise AdsorptionError("Cant align to z!")
                obj.resize(0)


    def endHandler(self,obj,method,method_locals,**kwargs):
        if isinstance(obj,Polymer):
            if(method.__name__ == 'align'):
                self._adjustMismatchPolymerSlab(obj)
    #endregion

    def resize(self,n,m):
        '''
        Resizes the substrate in the XY plane.

        Parameters
        ----------
        n : integer
            Number of repetitions of the substrate's unit cell along the first lattice 
            vector.
        m : integer
            Number of repetitions of the substrate's unit cell along the second lattice 
            vector.
        '''        
        self._components["substrate"].resize(n,m)
        self._update()
            
    def _update(self):
        '''
        Updates attributes of the Adsorption object when atomic coordinates or 
        cell shape or cell size change.
        '''

        self._a0=self._components["substrate"]._a0
        self._m=self._components["substrate"]._m
        self._n=self._components["substrate"]._n
        self._maxm=self._components["substrate"]._maxm
        self._maxn=self._components["substrate"]._maxn

        for molecule in self.adsorbed_molecules:
            dz=0.0
            
            #self.set_separation(update=False)

            if molecule in self._components["molecules@top"] and \
                (molecule.minz-self._components["substrate"]._top)<self._minsep:
                dz=self._components["substrate"]._top+self._minsep-molecule._minz
            elif molecule in self._components["molecules@bottom"] and \
                (self._components["substrate"]._bottom-molecule._maxz)<self._minsep:
                dz=self._components["substrate"]._bottom-self._minsep-molecule._maxz
                
            if abs(dz)>0.0:
                molecule.displace(array([0.0,0.0,dz]))
        for polymer in self.adsorbed_polymers:
            dz=0.0
            
            #self.set_separation(update=False)

            if polymer in self._components["polymers@top"] and \
                (polymer.minz-self._components["substrate"]._top)<self._minsep:
                dz=self._components["substrate"]._top+self._minsep-polymer._minz
            elif polymer in self._components["polymers@bottom"] and \
                (self._components["substrate"]._bottom-polymer._maxz)<self._minsep:
                dz=self._components["substrate"]._bottom-self._minsep-polymer._maxz
                
            if abs(dz)>0.0:
                polymer.displace(array([0.0,0.0,dz]))
        
        valid_z=array([atom._x[2] for atom in self._atoms if atom._active])
        self._top=max(valid_z,axis=0)
        self._bottom=min(valid_z,axis=0)
        self._latvec[2]=array([0.0,0.0,(self._top-self._bottom+self._vaccuum)
                               /self._a0])

    def __str__(self):
        '''
        Returns the ID of the Adsorption object.

        Returns
        -------
        String
            A string containing the ID of the Adsorption object.
        '''
        return "<Adsorption object> ID: %d" % self._id


    '''
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
    '''     
    @property
    def substrate(self):
        return self._components["substrate"]
    
    @property
    def adsorbed_molecules(self):
        return tuple(self._components["molecules@top"]+
                     self._components["molecules@bottom"])
    
    @property
    def adsorbed_polymers(self):
        return tuple(self._components["polymers@top"]+
                     self._components["polymers@bottom"])
    @property
    def minimum_separation(self):
        return self._minsep
    
    @minimum_separation.setter
    def minimum_separation(self,val):
        if isinstance(val,(float,int)) and val>0.0:
            self._minsep=val
            
            self._update()
        else:
            raise AdsorptionError("Minimum separation must be a number greater than zero!")
    
    @property
    def vaccuum(self):
        return self._vaccuum
    
    @vaccuum.setter
    def vaccuum(self,val):
        if isinstance(val,(float,int)) and val>0.0:
            self._vaccuum=val            
            self._latvec[2]=array([0.0,0.0,(self._top-self._bottom+val)/self._a0])
        else:
            raise AdsorptionError("Vaccuum must be provided as a number greater than zero!")
            
    @property
    def top(self):
        return self._top
    
    @property
    def bottom(self):
        return self._bottom
    
    @property
    def a0(self):
        return self.components["substrate"]._a0
    
    @a0.setter
    def a0(self,val):
        self._components["substrate"].a0=val
        
        self._update()
        
    @property
    def lattice_vectors(self):
        return self._components["substrate"]._latvec
    
    @property
    def origin(self):
        return self._components["substrate"]._origin
    
    @origin.setter
    def origin(self,val):
        if isinstance(val,(list,tuple)):
            val=array(val)
            
        if isinstance(val,ndarray) and val.shape[0]==3:
            disp=val-self._components["substrate"]._origin  
            
            self._components["substrate"].displace(disp)
        
            for molecule in self.adsorbed_molecules:
                molecule.displace(disp)
                
            self._update()
        else:
            raise AdsorptionError("The new origin must be a Numpy array with three components!")