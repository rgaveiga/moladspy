from __future__ import print_function
from __atom__ import Atom
from __atomcollection__ import AtomCollection
from __exception__ import BasicException
from numpy import array,ndarray,min,max,sum,cos,sin,radians,dot,zeros
from copy import deepcopy
from multipledispatch import dispatch

class PolymerError(BasicException):
    pass

class Polymer(AtomCollection):
    __slots__=["_vaccuum","_maxx","_maxy","_maxz","_minx","_miny","_minz","_anchors"]
    
    @dispatch(str,vaccuum=(int,float))
    def __init__(self,label,vaccuum=10.0,repetition_orientation=0,metric_method='min_max'):
        '''
        Object initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the type of polymer, for instance, 
            "Polyethylene" or "(C2H4)n".
        vaccuum : float, optional
            Vaccuum separating the polymer from its images in a periodic 
            boundary conditions scheme. The default is 10.0 Angstroms.
            
        '''
        super().__init__()
        
        self._orientation = self._old_orientation = repetition_orientation # x by standard
        self._metric_method = metric_method
        self._dic_metric_method = None

        if len(label)>0:
            self._label=label
        else:
            raise PolymerError("'label' must be a non-empty string!")
            
        if vaccuum>0.0:
            self._vaccuum=vaccuum
        else:
            raise PolymerError("'vaccuum' must be a number greater than zero!")
                    
        self._maxx=self._maxy=self._maxz=None       # Maximum Cartesian coordinates
        self._minx=self._miny=self._minz=None       # Minimum Cartesian coordinates
        self._anchors={"com":array([0.0,0.0,0.0])}  # Anchor points for translations and rotations

    @dispatch(str,str,file_type=str,vaccuum=(int,float)) 
    def __init__(self,label,file_name,file_type="XYZ",vaccuum=10.0,repetition_orientation=0, metric_method='min_max'):
        '''
        Object initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the type of polymer, for instance, 
            "benzene" or "C6H12".
        file_name : string
            Name of the file containing the polymer structure.
        file_type : string
            Type of file containing the polymer coordinates. The default is 
            "XYZ", which is the only type currently supported.
        vaccuum : float, optional
            Vaccuum separating the polymer from its images in a periodic 
            boundary conditions scheme. The default is 10.0 Angstroms.

        '''
        super().__init__()
        
        
        self._metric_method = metric_method
        self._dic_metric_method = None
        self._orientation = self._old_orientation = repetition_orientation # x by standard

        if len(label)>0:
            self._label=label
        else:
            raise PolymerError("'label' must be a non-empty string!")
            
        self._anchors={"com":array([0.0,0.0,0.0])}   # Anchor points for translations and rotations
            
        if vaccuum>0.0:
            self._vaccuum=vaccuum
        else:
            raise PolymerError("'vaccuum' must be a number greater than zero!")

        if len(file_name)>0:       
            if file_type.lower()=='xyz':
                self._read_xyz(file_name)
            else:
                raise PolymerError("'file_type' must be 'XYZ'!")
        else:
            raise PolymerError("'file_name' must be a valid file name!")
        
        
    @dispatch(str,list,vaccuum=(int,float))    
    def __init__(self,label,atom_list,vaccuum=10.0,metric_method='min_max',repetition_orientation=0):
        '''
        Object initialization.

        Parameters
        ----------
        label : string
            A label allowing you to identify the type of polymer, for instance, 
            "benzene" or "C6H12".
        atom_list : Python list
            List of Atom objects or atom IDs to be added to the polymer.
        vaccuum : float, optional
            Vaccuum separating the polymer from its images in a periodic 
            boundary conditions scheme. The default is 10.0 Angstroms.

        '''
        super().__init__()
        self._metric_method = metric_method
        self._dic_metric_method = None
        self._orientation = self._old_orientation = repetition_orientation # x by standard

        if len(label)>0:
            self._label=label
        else:
            raise PolymerError("'label' must be a non-empty string!")
            
        self._anchors={"com":array([0.0,0.0,0.0])}   # Anchor points for translations and rotations
        self._origin=None
            
        if vaccuum>0.0:
            self._vaccuum=vaccuum
        else:
            raise PolymerError("'vaccuum' must be a number greater than zero!")        
                        
        if len(atom_list)>0:
            for atom in atom_list:
                if isinstance(atom,(Atom,int)):
                    self.add_atom(atom,loc=(0,0,0),update=False)
                else:
                    print("WARNING! An element in the atom list must be either an Atom object or an atom ID!")
        else:
            raise PolymerError("'atom_list' must be a non-empyt list!")
            
        self._update()

    def add_anchor(self,anchor,pos):
        '''
        add_anchor(anchor,pos) -> adds a new anchor point that can be used as 
            the reference point of the polymer for translations and rotations.

        Parameters
        ----------
        anchor : string
            Name of the anchor point to be added to the polymer.
        pos : Numpy array
            Cartesian coordinates of the anchor point. It can also be provided 
            as a Python list or tuple.

        '''
        if isinstance(pos,(list,tuple)):
            pos=array(pos)
        
        if isinstance(anchor,str) and len(anchor)>0:
            if isinstance(pos,ndarray) and pos.shape[0]==3:
                self._anchors[anchor]=pos.astype(float)
            else:
                raise PolymerError("'pos' must be an array with three components!")
        else:
            raise PolymerError("'anchor' must be a non-empty string!")
            
    def remove_anchor(self,anchor):
        '''
        remove_anchor(anchor) -> removes an anchor point.

        Parameters
        ----------
        anchor : string
            Name of the anchor point to be removed from the polymer.

        '''
        if isinstance(anchor,str) and anchor in self._anchors.keys():
            if anchor=="com":
                raise PolymerError("The center of mass cannot be removed!")                
            else:
                self._anchors.pop(anchor)
        else:
            raise PolymerError("'anchor' is not a valid anchor point!")
        
    def write_xyz(self,file_name="polymer.xyz"):
        '''
       Saves the atomic coordinates of the polymer into an XYZ file.

        Parameters
        ----------
        file_name : string, optional
            Name of the XYZ file. The default is "polymer.xyz".

        '''
        super().write_xyz(file_name,ucell=True)          
                
    def write_pw_input(self,file_name="polymer.in",pseudopotentials={},pwargs={}):
        '''
        Creates a basic input file for geometry relaxation of the polymer using 
        the pw.x code found in the Quantum Espresso package.

        Parameters
        ----------
        file_name : string, optional
            Name of the input file. The default is "polymer.in".        
        pseudopotentials : Python dictionary, optional
            Specify for each element provided as a key the name of a pseudopotential 
            file given as the corresponding value. The default is {}.        
        pwargs : Python dictionary.
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

        '''
        pwargs["kvec"]=[1,1,1,0,0,0]
        ucell=True
        
        super().write_pw_input(file_name,ucell,pseudopotentials,pwargs)

    def copy(self):
        '''
        Returns a copy of the Molecule object.

        Returns
        -------
        Molecule object
            Copy of the Molecule object.
        '''
        newmol=super().copy()
        newmol._anchors=deepcopy(self._anchors)
        
        return newmol
        
    def displace(self,disp):
        '''
        displace(disp) -> rigidly displaces the polymer.

        Parameters
        ----------
        disp : Numpy array
            Displacement vector. It can also be provided as a Python list or tuple.
        '''
        super().displace(disp)
        
        for key in self._anchors.keys():
            if key!="com":
                self._anchors[key]+=disp
            
    def move_to(self,x,anchor="com"):
        '''
        move_to(x,anchor) -> rigidly moves the polymer such that the anchor 
        point is located at 'x'.

        Parameters
        ----------
        x : Numpy array
            New position of the polymer's anchor point. It can also be provided 
            as a Python list or tuple.
        anchor : string, optional
            Anchor point used as reference for translations. The default is 
            'com', which means the polymer's center of mass.

        '''
        if isinstance(x,(list,tuple)):
            x=array(x)
        
        if not (isinstance(anchor,str) and anchor in self._anchors.keys()):
            raise PolymerError("'anchor' is not a valid anchor point!")
        
        if isinstance(x,ndarray) and x.shape[0]==3:
            disp=x.astype(float)-self._anchors[anchor]
            
            self.displace(disp)
        else:
            raise PolymerError("'x' must be an array with three components!")

    def _rotate(self,theta,phi,psi,anchor="com"):
        '''
        _rotate(theta,phi,psi,anchor) -> rotates the polymer around an anchor point. 
        This method is private because it allows to rotate the polymer object with all degrees of freedom.

        Parameters
        ----------
        theta : float, optional
            First rotation angle (around Y axis), in degrees.
        phi : float, optional
            Second rotation angle (around X axis), in degrees.
        psi : float, optional
            Third rotation angle (around Z axis), in degrees.
        anchor : string, optional
            Anchor point around which the polymer will be rotated. The default 
            is 'com', which means the polymer's center of mass.

        '''
        if isinstance(theta,(float,int)) and \
            isinstance(phi,(float,int)) and \
            isinstance(psi,(float,int)):       
            theta=radians(theta)
            phi=radians(phi)
            psi=radians(psi)
        else:
            raise PolymerError("Rotation angles must be provided as numbers!")
            
        if not (isinstance(anchor,str) and anchor in self._anchors.keys()):
            raise PolymerError("'anchor' is not a valid anchor point!")
            
        a11=cos(theta)*cos(psi)
        a12=-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)
        a13=sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)
        a21=cos(theta)*sin(psi)
        a22=cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)
        a23=-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)
        a31=-sin(theta)
        a32=sin(phi)*cos(theta)
        a33=cos(phi)*cos(theta)
        rotmatrix=array([[a11,a12,a13],
                         [a21,a22,a23],
                         [a31,a32,a33]])
        rotpoint=self._anchors[anchor]
        
        for atom in self._atoms:
            atom._x=dot((atom._x-rotpoint),rotmatrix)+rotpoint
                
        for key in self._anchors.keys():
            if key!=anchor:
                self._anchors[key]=dot((self._anchors[key]-rotpoint),
                                       rotmatrix)+rotpoint                
        
        self._update()
    
    def align(self, axis):
        """
        Aligns the polymer along the specified axis.

        Parameters
        ----------
        axis : str
            The axis along which to align the polymer ('x', 'y', or 'z').
        """
        _axis_dic = {
            'x':0,
            'y':1,
            'z':2
        }
        # Check if axis is correctly given
        try:
            new_orientation = _axis_dic[axis]
        except:
            raise PolymerError(f"Invalid axis '{axis}'. Please choose 'x', 'y', or 'z'.")

        # Just align if orientation changes
        if(new_orientation != self._orientation):
            # Hold previous orientation
            self.old_orientation = self._orientation
            self._orientation = new_orientation
            
            # Dictionary to map the axis of the rotation angles (theta, phi, psi)
            # This values are examples and must be adjusted by necessity
            # Theta: Axis Y Rotation
            # Phi: Axis X Rotation
            # Psi: Axis Z Rotation
            
            
            '''
            theta : float, optional
                First rotation angle (around Y axis), in degrees.
            phi : float, optional
                Second rotation angle (around X axis), in degrees.
            psi : float, optional
                Third rotation angle (around Z axis), in degrees.
            anchor : string, optional
                Anchor point around which the polymer will be rotated. The default 
                is 'com', which means the polymer's center of mass.

            '''
            # Identify which angle should be rotate by the sum of the old orientation and the new one
            choose_angle = {
                1:2,
                3:0,
                2:1
            }
            right_angle = choose_angle[self._orientation + self._old_orientation]
            rotation_angles = array([0.0,0.0,0.0])
            di = self._orientation - self._old_orientation
            
            # Just trust this works
            angle_sign = (-1)**abs(di)*di/abs(di)
            rotation_angles[right_angle] = angle_sign*90.0
            phi, theta, psi = rotation_angles
            print(phi,theta,psi)
            self._rotate(theta,phi,psi,anchor="origin")

            self._update()

            
            # TODO: Permutate new align, change max and min of old axis by new one
            # Obtain polymer orientation from the maximum diff max - min. TODO: Reimplement
            
    def rotate(self,angle_degrees):
        """
        Rotate the polymer in its orientation axis center,

        Parameters
        ----------
        angle_degrees : str
            The angle in degrees following right-hand rule.
        """
        self._calculate_polymer_rotation_axis_center()
        rotation_angles = array([angle_degrees if i==self._orientation else 0.0 for i in range(3)])
        phi, theta, psi = rotation_angles
        self._rotate(theta,phi,psi,anchor='polymer_axis')

    def _calculate_polymer_rotation_axis_center(self,metric_method='min_max'):
        '''
        Calculates the center position in the polymer perpendicular plane through a given metric
        
        Parameters
        ----------
        metric : str
            The metric that will be used as standard to calculate the position of the rotation axis in the orientation perpedicular plane.
        """
        Metrics accepted:
            - 'min_max': Uses the center point between the min and max of a coordinate
            - 'geometric_center': Makes the average of atom centroids
            - 'mass_center': Uses the center of mass coordinates

        '''
        point = array([0.0,0.0,0.0])
        
        if(self._metric_method == None):
            self._metric_method = {
                'max_min':lambda: self._anchors["middle_point"],
                'mass_center': lambda: self._anchors["com"],
                'geometric_center': lambda: self._anchors["geometric_center"]
            }
        point = self._metric_method[metric_method]()
        
        point[self._orientation] = 0.0
        self._anchors['polymer_axis'] = point

    def resize(self,n_axis):
        '''
        Resizes the polymer in the orientation plane by self._orientation index
        '''
        
        n, m, l = [0,0,0]
        if(self._orientation == 0):
            n = n_axis
        elif(self._orientation == 1):
            m = n_axis
        elif(self._orientation == 2):
            l = n_axis
        
        super().resize(n,m,l)
                       
    def _update(self):
        '''
        _update() -> simultaneously updates the polymer's center of mass and 
        the values of its extremities.
        '''
        

        valid_coords=array([atom._x for atom in self.active_atoms])
        atomic_mass=array([atom.atomic_mass for atom in self.active_atoms])

        if valid_coords.shape[0]>0:
            self._maxx,self._maxy,self._maxz=max(valid_coords,axis=0)
            self._minx,self._miny,self._minz=min(valid_coords,axis=0)
            totmass=sum(atomic_mass)
            self._anchors["com"]=sum(atomic_mass*valid_coords.transpose(),axis=1)/totmass
            self._anchors["middle_point"]=array([(self.maxx + self.minx)/2.0,(self.maxy + self.miny)/2.0,(self.maxz + self.minz)/2.0])
            self._anchors["geometric_center"]=sum(valid_coords.transpose(),axis=1)/len(self)
            self._anchors["origin"] = array([0.0,0.0,0.0])

            # TODO: Calculate only vacuum on perpendicular orientations
            '''
            if(self._old_orientation != self._orientation):
                tmpvec = self._latvec[self._orientation]
                self._latvec[self._orientation] = self._latvec[self._old_orientation]
                self._latvec[self._old_orientation] = tmpvec
            '''

            
            vaccuum_vec =array([0.0,0.0,0.0]),array([0.0,0.0,0.0]),array([0.0,0.0,0.0])
            # Iterate in every axis
            for i in range(3):
                if(i != self._orientation):
                    vaccuum_vec[i][i] = self._vaccuum

            self._latvec=array([[self._maxx-self._minx,0.0,0.0],
                                [0.0,self._maxy-self._miny,0.0],
                                [0.0,0.0,self._maxz-self._minz]])
            self._latvec = self._latvec + vaccuum_vec
            print("Vaccum Vector")
            print(vaccuum_vec)
            self._origin=array([self._minx,self._miny,self._minz])


        else:
            self._maxx=self._maxy=self._maxz=None
            self._minx=self._miny=self._minz=None
            self._anchors["com"]=array([0.0,0.0,0.0])
            self._latvec=array([[1.0,0.0,0.0],
                                [0.0,1.0,0.0],
                                [0.0,0.0,1.0]])
            self._origin=array([0.0,0.0,0.0])
            
    def __str__(self):
        '''
        __str__() -> returns the type and the ID of the Polymer object.

        Returns
        -------
        String.
            A string containing the type and ID of the Polymer object.
        '''
        return "<Polymer object> Type: %s; ID: %d" % (self._label,self._id)
    
    '''
    Properties
    ----------
    a0 : float, readonly
        Lattice parameter. Meaningless for polymers, it is still necessary  
        when writing extended XYZ files or PW input files. It is set to 1.
    lattice_vectors : Numpy array, readonly
        Lattice vectors. Meaningless for polymers, it is still necessary  
        when writing extended XYZ files or PW input files. The vectors are 
        determined by adding the prescribed length of vaccuum to the extremities 
        of the polymer.
    anchors : Python dictionary, readonly
        Dictionary of anchor points.
    center_of_mass : Numpy array, readonly
        Center of mass of the polymer.
    maxx,minx,maxy,miny,maxz,minz : floats, readonly
        Polymer extremities.
    origin : Numpy array
        Position with minimum values of X, Y, and Z coordinates of polymer.
    '''
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
    
    @property
    def min(self):
        return array([self._minx,self._miny, self._minz])
    
    @property
    def max(self):
        return array([self._maxx,self._maxy, self._maxz])