from MolAdsPy.__atom__ import Atom,Species
from MolAdsPy.__exception__ import BasicException
from multipledispatch import dispatch
from numpy import array,ndarray,dot,sqrt
from copy import copy
from abc import ABC,abstractmethod
import os.path

class AtomCollectionError(BasicException):
    pass

class AtomCollection(ABC):
    _currid=None        # ID available for a new atom collection object
    _instances=None     # List of atom collection objects created so far
    __slots__=["_ucell","_atoms","_loc","_species","_a0","_latvec","_n","_m","_l",
               "_maxn","_maxm","_maxl","_id","_origin","_belongs_to","_label"]
    
    def __init__(self):
        '''
        Object initialization.
        '''
        self._label=""          # Label used as a convenient name of the atom collection object
        self._belongs_to=None   # Hybrid structure to which an atom collection object belongs
        self._ucell=[]          # List of Atom objects in the unit cell
        self._atoms=[]          # List of Atom objects in the atomic structure
        self._loc={}            # Location of the atoms in the supercell
        self._species=[]        # List of atomic species in the atom collection object
        self._a0=1.0            # Lattice constant
        self._latvec=array([[1.0,0.0,0.0],
                            [0.0,1.0,0.0],
                            [0.0,0.0,1.0]])     # Lattice vectors
        self._n=self._m=self._l=0               # Number of translations along the first, second and third lattice vectors
        self._maxn=self._maxm=self._maxl=0      # Maximum number of translations for resizing purposes
        self._origin=array([0.0,0.0,0.0])       # Position with minimum values of X, Y, and Z coordinates of the structure.        
        
        type(self)._append(self)
        
    @dispatch(int,loc=(tuple,list),update=bool)
    def add_atom(self,atomid,loc=(0,0,0),update=True):
        '''
        Adds an atom to the structure.

        Parameters
        ----------
        atomid : integer
            ID of the atom to be added to the structure.
        loc : Python tuple, optional
            Location of the atom in the supercell. The default is (0,0,0).
        update : logical, optional
            Whether or not update some structure attributes if atomic coordinates 
            or cell shape or cell size change.
        '''
        if atomid>=0 and atomid<Atom._curratomid:
            for atom in Atom._instances:
                if atomid==atom._id and atom._belongs_to is None:
                    self._atoms.append(atom)
                    atom._belongs_to=self
                    self._loc[atom._id]=tuple(loc)
                    
                    if loc==(0,0,0):
                        self._ucell.append(atom)
                        
                    if not atom._symbol in self._species:
                        self._species.append(atom._symbol)
                        
                    if update:
                        self._update()
                            
                    break
                else:
                    raise AtomCollectionError("Atom already belongs to another structure!")
            else:
                raise AtomCollectionError("Atom ID not found!")
        else:
            raise AtomCollectionError("'atomid' must be an integer greater than or equal to zero!")
                
    @dispatch(Atom,loc=(tuple,list),update=bool)
    def add_atom(self,atom,loc=(0,0,0),update=True):
        '''
        Adds an atom to the structure.

        Parameters
        ----------
        atom : Atom object
            Atom to be added to the structure.
        loc : Python tuple, optional
            Location of the atom in the supercell. The default is (0,0,0).
        update : logical, optional
            Whether or not update some structure attributes if atomic coordinates 
            or cell shape or cell size change.
        '''
        if atom._belongs_to is None:
            self._atoms.append(atom)
            atom._belongs_to=self
            self._loc[atom._id]=tuple(loc)
            
            if loc==(0,0,0):
                self._ucell.append(atom)
                
            if not atom._symbol in self._species:
                self._species.append(atom._symbol)
            
            if update:
                self._update()
        else:
            raise AtomCollectionError("Atom already belongs to another structure!")

    @dispatch(int,update=bool)            
    def remove_atom(self,atomid,update=True):
        '''
        Removes an atom from the structure.

        Parameters
        ----------
        atomid : integer
            ID of the atom to be removed from the structure.
        update : logical, optional
            Whether or not update some structure attributes if atomic coordinates 
            or cell shape or cell size change.
        '''     
        if atomid>=0 and atomid<Atom._curratomid:
            for atom in self._atoms:
                if atom._id==atomid:
                    self._atoms.remove(atom)
                    
                    if self._loc[atom._id]==(0,0,0):
                        self._ucell.remove(atom)
                    
                    self._loc.pop(atom._id)
                    
                    if update:
                        self._update()
                    
                    atom._belongs_to=None
                                        
                    for other_atom in self._atoms:
                        if atom._symbol==other_atom._symbol:
                            break
                    else:
                        self._species.remove(atom._symbol)
                    
                    break
            else:
                raise AtomCollectionError("Atom ID not found!")
        else:
            raise AtomCollectionError("'atomid' must be an integer greater than or equal to zero!")
                
    @dispatch(Atom,update=bool)            
    def remove_atom(self,atom,update=True):
        '''
        Removes an atom from the structure.

        Parameters
        ----------
        atom : Atom object
            Atom to be removed from the structure.
        '''     
        if atom in self._atoms:
            self._atoms.remove(atom)
            
            if self._loc[atom._id]==(0,0,0):
                self._ucell.remove(atom)           
            
            self._loc.pop(atom._id)
            
            if update:
                self._update()
            
            atom._belongs_to=None
                                
            for other_atom in self._atoms:
                if atom._symbol==other_atom._symbol:
                    break
            else:
                self._species.remove(atom._symbol)                
        else:
            raise AtomCollectionError("Atom not found!")
            
    def displace(self,disp):
        '''
        displace(disp) -> rigidly displaces the atomic structure.

        Parameters
        ----------
        disp : Numpy array
            Displacement vector. It can also be provided as a Python list or tuple.
        '''
        if isinstance(disp,(list,tuple)):
            disp=array(disp)
        
        if isinstance(disp,ndarray) and disp.shape[0]==3:
            disp=disp.astype(float)
            
            for atom in self._atoms:
                atom._x+=disp

            self._update()
        else:
            raise AtomCollectionError("'disp' must be an array with three elements!")

    @dispatch(int)
    def location(self,atomid):
        '''
        Returns three indices correspoding to the location of the atom in the 
        supercell.

        Parameters
        ----------
        atomid : integer
            ID of an atom in the supercell.

        Returns
        -------
        Python tuple
            Indices corresponding to the location of the atom in the supercell. 
        '''
        if atomid in self._loc.keys():
            return self._loc[atomid]
        else:
            raise AtomCollectionError("Atom ID not found!")

    @dispatch(Atom)
    def location(self,atom):
        '''
        Returns three indices correspoding to the location of the atom in the 
        supercell.

        Parameters
        ----------
        atom : Atom object
            An atom in the supercell.

        Returns
        -------
        Python tuple
            Indices corresponding to the location of the atom in the supercell. 
        '''
        if atom._id in self._loc.keys():
            return self._loc[atom._id]
        else:
            raise AtomCollectionError("Atom not found!")
            
    def write_xyz(self,file_name="coords.xyz",ucell=False):
        '''
        Saves the atomic coordinates of the structure into an XYZ file.

        Parameters
        ----------
        file_name : string, optional
            Name of the XYZ file. The default is "coords.xyz".
        ucell : logical, optional
            Write only the coordinates of atoms in the unit cell. The default is 
            False.
        '''
        if not (isinstance(file_name,str) and len(file_name)>0):
            raise AtomCollectionError("'file_name' must be a non-empty string!")
        
        natoms=0
        atompos=""
        n=self._n+1 if not ucell else 1
        m=self._m+1 if not ucell else 1
        l=self._l+1 if not ucell else 1
        
        for atom in self._atoms:
            i,j,k=self._loc[atom._id]
            
            if atom._active and i<n and j<m and k<l:
                if not ucell or atom in self._ucell:
                    natoms+=1
                    atompos+="%s %.6f %.6f %.6f\n" % (atom._symbol,atom._x[0],
                                                      atom._x[1],atom._x[2])
                
        if natoms==0:
            raise AtomCollectionError("There is no active atom! Nothing to do!")   

        with open(file_name,"w") as f:
            f.write("%d\n" % natoms)
            f.write('Lattice="%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f" Properties=species:S:1:pos:R:3\n'
                    % (self._a0*n*self._latvec[0][0],
                       self._a0*n*self._latvec[0][1],
                       self._a0*n*self._latvec[0][2],
                       self._a0*m*self._latvec[1][0],
                       self._a0*m*self._latvec[1][1],
                       self._a0*m*self._latvec[1][2],
                       self._a0*l*self._latvec[2][0],
                       self._a0*l*self._latvec[2][1],
                       self._a0*l*self._latvec[2][2]))
            f.write(atompos)

    def write_pw_input(self,file_name="pw.in",ucell=False,pseudopotentials={},
                       pwargs={}):
        '''
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
        '''
        if not (isinstance(file_name,str) and len(file_name)>0):
           raise AtomCollectionError("'file_name' must be a non-empty string!")
           
        ang2bohr=1.8897259886
        natoms=0
        inpstr=""
        atompos=""
        n=self._n+1 if not ucell else 1
        m=self._m+1 if not ucell else 1
        l=self._l+1 if not ucell else 1
        
        for atom in self._atoms:
            i,j,k=self._loc[atom._id]
            
            if atom._active and i<n and j<m and k<l:
                if not ucell or atom in self._ucell:
                    natoms+=1
                    atompos+="%s %.6f %.6f %.6f %d %d %d\n" % (atom._symbol,
                                                               atom._x[0],
                                                               atom._x[1],
                                                               atom._x[2],
                                                               int(not(atom._fixed[0])),
                                                               int(not(atom._fixed[1])),
                                                               int(not(atom._fixed[2])))
                
        if natoms==0:
            raise AtomCollectionError("There is no active atom! Nothing to do!")
        
        inpstr="&CONTROL\n"
        
        if "calculation" in pwargs:            
            if pwargs["calculation"] in ("relax","vc-relax"):
                calculation=pwargs["calculation"]                
                inpstr+="    calculation='%s',\n" % pwargs["calculation"]
            else:
                raise AtomCollectionError("Only 'relax' and 'vc-relax' calculations are allowed!")
        else:
            calculation="relax"
            inpstr+="    calculation='relax',\n"
            
        inpstr+="    restart_mode='from_scratch',\n"
        inpstr+="    pseudo_dir='.',\n"
        inpstr+="    prefix='%s',\n" % (os.path.splitext
                                        (os.path.basename(file_name))[0])
        inpstr+="    nstep=500,\n/\n"
        inpstr+="&SYSTEM\n"
        
        if "ibrav" in pwargs and pwargs["ibrav"]!=0:
            ibrav=int(pwargs["ibrav"])
            inpstr+="    ibrav=%d,\n" % int(pwargs["ibrav"])
            latpar=self._a0*sqrt(self._latvec[0].dot(self._latvec[0]))*n*ang2bohr
            inpstr+="    celldm(1)=%.10f,\n" % latpar
            
            if pwargs["ibrav"] in (4,6,7):
                if "celldm(3)" in pwargs:
                    inpstr+="    celldm(3)=%.10f,\n" % pwargs["celldm(3)"]
                else:
                    raise AtomCollectionError("'celldm(3)' must be provided!")
            elif pwargs["ibrav"] in (-5,5):
                if "celldm(4)" in pwargs:
                    inpstr+="    celldm(4)=%.10f,\n" % pwargs["celldm(4)"]
                else:
                    raise AtomCollectionError("'celldm(4)' must be provided!")
            elif pwargs["ibrav"] in (-9,8,9,10,11,91):
                if all(key in pwargs["ibrav"] for key in 
                       ("celldm(2)","celldm(3)")):
                    inpstr+="    celldm(2)=%.10f,\n" % pwargs["celldm(2)"]
                    inpstr+="    celldm(3)=%.10f,\n" % pwargs["celldm(3)"]
                else:
                    raise AtomCollectionError("'celldm(2)' and 'celldm(3)' must be provided!")
            elif pwargs["ibrav"] in (12,13):
                if all(key in pwargs["ibrav"] for key in 
                       ("celldm(2)","celldm(3)","celldm(4)")):
                    inpstr+="    celldm(2)=%.10f,\n" % pwargs["celldm(2)"]
                    inpstr+="    celldm(3)=%.10f,\n" % pwargs["celldm(3)"]
                    inpstr+="    celldm(4)=%.10f,\n" % pwargs["celldm(4)"]
                else:
                    raise AtomCollectionError("'celldm(2)', 'celldm(3)' and 'celldm(4)' must be provided!")
            elif pwargs["ibrav"] in (-13,-12):
                if all(key in pwargs["ibrav"] for key in 
                       ("celldm(2)","celldm(3)","celldm(5)")):
                    inpstr+="    celldm(2)=%.10f,\n" % pwargs["celldm(2)"]
                    inpstr+="    celldm(3)=%.10f,\n" % pwargs["celldm(3)"]
                    inpstr+="    celldm(5)=%.10f,\n" % pwargs["celldm(5)"]
                else:
                    raise AtomCollectionError("'celldm(2)', 'celldm(3)' and 'celldm(5)' must be provided!")
            elif pwargs["ibrav"]==14:
                if all(key in pwargs["ibrav"] for key in 
                       ("celldm(2)","celldm(3)","celldm(4)","celldm(5)","celldm(6)")):
                    inpstr+="    celldm(2)=%.10f,\n" % pwargs["celldm(2)"]
                    inpstr+="    celldm(3)=%.10f,\n" % pwargs["celldm(3)"]
                    inpstr+="    celldm(4)=%.10f,\n" % pwargs["celldm(4)"]
                    inpstr+="    celldm(5)=%.10f,\n" % pwargs["celldm(5)"]
                    inpstr+="    celldm(6)=%.10f,\n" % pwargs["celldm(6)"]
                else:
                    raise AtomCollectionError("'celldm(2)', 'celldm(3)', 'celldm(4)', 'celldm(5)' and 'celldm(6)' must be provided!")
        else:
            ibrav=0
            inpstr+="    ibrav=0,\n"
            inpstr+="    celldm(1)=%.10f,\n" % (self._a0*ang2bohr)
            
        inpstr+="    nat=%d,\n" % natoms
        inpstr+="    ntyp=%d,\n" % (len(self._species))
        
        if "ecutwfc" in pwargs:
            inpstr+="    ecutwfc=%.2f,\n" % pwargs["ecutwfc"]
        else:
            inpstr+="    ecutwfc=32.00,\n"
        
        if "ecutrho" in pwargs:
            inpstr+="    ecutrho=%.2f,\n" % pwargs["ecutrho"]
        else:
            inpstr+="    ecutrho=128.00,\n"
            
        if "occupations" in pwargs:
            inpstr+="    occupations='%s',\n" % pwargs["occupations"]
            
            if pwargs["occupations"]=="smearing":
                if "smearing" in pwargs:
                    inpstr+="    smearing='%s',\n" % pwargs["smearing"]
                else:
                    inpstr+="    smearing='gaussian',\n"
                    
                if "degauss" in pwargs:
                    inpstr+="    degauss=%.4f,\n" % pwargs["degauss"]
                else:
                    inpstr+="    degauss=0.02"
        else:
            inpstr+="    occupations='fixed',\n"
            
        if "nspin" in pwargs and pwargs["nspin"]==2:
            inpstr+="    nspin=2,\n"
                
            if "starting_magnetization" in pwargs and \
                len(pwargs["starting_magnetization"])>0:
                    for i in range(len(pwargs["starting_magnetization"])):
                        if i<len(self._species):                           
                            inpstr+="    starting_magnetization(%d)=%.2f,\n" \
                                % (i+1,pwargs["starting_magnetization"][i])
                        else:
                            break
            else:
                inpstr+="    starting_magnetization(1)=1.00,\n"
        else:
            inpstr+="    nspin=1,\n"
            
        if "input_dft" in pwargs and pwargs["input_dft"]=="vdw-df":
            inpstr+="    input_dft='vdw-df',\n"
        
        inpstr+="/\n"
        inpstr+="&ELECTRONS\n"
        
        if "mixing_mode" in pwargs:
            inpstr+="    mixing_mode='%s',\n" % pwargs["mixing_mode"]
        else:
            inpstr+="    mixing_mode='plain',\n"
            
        if "mixing_beta" in pwargs and pwargs["mixing_beta"]>0.0:
            inpstr+="    mixing_beta=%.2f,\n" % pwargs["mixing_beta"]
        else:
            inpstr+="    mixing_beta=0.7,\n"
        
        inpstr+="    electron_maxstep=200,\n/\n"
        inpstr+="&IONS\n/\n"
        
        if calculation=="vc-relax":
            inpstr+="&CELL\n/\n"
        
        if ibrav==0:
            inpstr+="CELL_PARAMETERS alat\n"
            inpstr+="%.4f   %.4f   %4f\n" % (n*self._latvec[0][0],
                                             n*self._latvec[0][1],
                                             n*self._latvec[0][2])
            inpstr+="%.4f   %.4f   %4f\n" % (m*self._latvec[1][0],
                                             m*self._latvec[1][1],
                                             m*self._latvec[1][2])
            inpstr+="%.4f   %.4f   %4f\n" % (l*self._latvec[2][0],
                                             l*self._latvec[2][1],
                                             l*self._latvec[2][2])

        inpstr+="K_POINTS automatic\n"

        if "kvec" in pwargs and len(pwargs["kvec"])==6:
            inpstr+="%d %d %d %d %d %d\n" % (pwargs["kvec"][0],
                                             pwargs["kvec"][1],
                                             pwargs["kvec"][2],
                                             pwargs["kvec"][3],
                                             pwargs["kvec"][4],
                                             pwargs["kvec"][5])
        else:
            inpstr+="1 1 1 0 0 0\n"
        
        inpstr+="ATOMIC_SPECIES\n"
        
        for symbol in self._species:
            if symbol in pseudopotentials:
                pseudopotential=pseudopotentials[symbol]
            else:
                pseudopotential=symbol+".UPF"
                
            inpstr+="%s %.4f %s\n" % (symbol,Species[symbol]["atomic mass"],
                                      pseudopotential)
                
        inpstr+="ATOMIC_POSITIONS angstrom\n"
        inpstr+=atompos
            
        with open(file_name,"w") as f:
            f.write(inpstr)

    def resize(self,n,m,l):
        '''
        Resizes the atom collection structure.

        Parameters
        ----------
        n : integer
            Number of repetitions of the structure's unit cell along the first 
            lattice vector.
        m : integer
            Number of repetitions of the structure's unit cell along the second 
            lattice vector.
        l : integer
            Number of repetitions of the structure's unit cell along the third 
            lattice vector.
        '''
        if not (isinstance(n,int) and isinstance(m,int) and isinstance (l,int)
                and n>=0 and m>=0 and l>=0):
            raise AtomCollectionError("'n', 'm' and 'l' values must be integers greater than or equal to zero!")            
        elif n==self._n and m==self._m and l==self._l:
            raise AtomCollectionError("The current structure size was not changed!")
        elif n>self._maxn or m>self._maxm and l>self._maxl:
            print("WARNING: Atoms newly created will not be removed if the supercell is subsequently reduced!")
                         
            for i in range(n+1):
                for j in range(m+1):
                    for k in range(l+1):
                        if i<=self._maxn and j<=self._maxm and k<=self._maxl:
                            continue
                        
                        disp=dot(self._a0*array([i,j,k]),self._latvec)
                        
                        for atom in self._ucell:
                            newatom=atom.copy()
                            newatom._belongs_to=self                            
                            newatom._x+=disp
                            self._loc[newatom._id]=(i,j,k)
                            
                            self._atoms.append(newatom)
            
            self._maxn=max([self._maxn,n])
            self._maxm=max([self._maxm,m])
            self._maxl=max([self._maxl,l])
                
        self._n=n
        self._m=m
        self._l=l
        
        self._update()

    def copy(self):
        '''
        Returns a copy of the atom collection object.

        Returns
        -------
        Atom collection object
            Copy of the atom collection object.
        '''        
        newobj=copy(self)
        newobj._origin=self._origin.copy()
        newobj._species=self._species.copy()
        newobj._latvec=self._latvec.copy()
        newobj._atoms=[]
        newobj._ucell=[]
        newobj._loc={}
        newobj._belongs_to=None
               
        for atom in self._atoms:
            newatom=atom.copy()
            newatom._belongs_to=newobj
            
            newobj._atoms.append(newatom)
            
            newobj._loc[newatom._id]=copy(self._loc[atom._id])
            
            if atom in self._ucell:
                newobj._ucell.append(newatom)
        
        type(self)._append(newobj)
        
        return newobj
            
    def _read_xyz(self,file_name):
        '''
        Reads the atomic coordinates of the unit cell from a file in XYZ format.

        Parameters
        ----------
        filename : string
            Name of the XYZ file.
        '''
        if not (isinstance(file_name,str) and len(file_name)>0):
            raise AtomCollectionError("'file_name' must be a valid file name!")        
        
        with open(file_name,"r") as f:
            lines=f.readlines()
            
        count=0
        extxyz=-1

        for line in lines:
            if count==0:
                count=1
                
                continue
            elif count==1:
                extxyz=line.find("Lattice")
                
                if extxyz>=0:
                    line=line.replace('Lattice','')
                    line=line.replace('=','')
                    line=line.replace('"','')
                    line=line.replace("'","")
                    l=line.split()
                    self._latvec=array([[float(l[0]),float(l[1]),float(l[2])],
                                        [float(l[3]),float(l[4]),float(l[5])],
                                        [float(l[6]),float(l[7]),float(l[8])]])
                else:
                    print("WARNING: Lattice vectors not found in the XYZ file!")
                    
                count=2
                    
                continue
                
            l=line.split()
        
            if len(l)>=4:
                if l[0] in Species:
                    self.add_atom(Atom(l[0],x=array([float(l[1]),float(l[2]),
                                                     float(l[3])])),loc=(0,0,0),
                                  update=False)
                else:                    
                    print("WARNING: Species '%s' and atoms of this species must be added manually!" % l[0])                              
                    
        print("INFO: %d Atom objects have just been created from file '%s'!" 
              % (len(self._atoms),file_name))
        
        self._update()
    
    def __getitem__(self,idx):
        '''
        Returns the Atom object corresponding to the specified index in the 
        structure's atom list.

        Parameters
        ----------
        idx : integer
            Atom index in the structure's atom list.

        Returns
        -------
        Atom
            Atom object corresponding to the specified index in the structure's 
            atom list.
        '''
        maxidx=len(self._atoms)
        
        if isinstance(idx,int) and idx>=-maxidx and idx<maxidx:
            return self._atoms[idx]
        else:
            raise AtomCollectionError("Invalid index!")

    @dispatch(int,int)
    def __setitem__(self,idx,atomid):
        '''
        Assigns an Atom object to the corresponding index in the structure's atom 
        list.

        Parameters
        ----------
        idx : integer
            Index in the structure's atom list.
        atomid : integer
            ID of the atom to be assigned to the corresponding index in the 
            structure's atom list.
        '''        
        maxidx=len(self._atoms)
        
        if idx>=-maxidx and idx<maxidx:
            if atomid>=0 and atomid<Atom._curratomid:
                for atom in Atom._instances:
                    if atomid==atom._id and atom._belongs_to is None:
                        self._atoms[idx]._belongs_to=None
                        
                        for other_atom in self._atoms:
                            if self._atoms[idx]._symbol==other_atom._symbol:
                                break
                        else:
                            self._species.remove(self._atoms[idx]._symbol)
                        
                        self._atoms[idx]=atom
                        atom._belongs_to=self
                        
                        if not atom._symbol in self._species:
                            self._species.append(atom._symbol)
                                
                        self._update()  
                        
                        break
                    else:
                        raise AtomCollectionError("Atom already belongs to another structure!")
                else:
                    raise AtomCollectionError("Atom ID not found!")
            else:
                raise AtomCollectionError("'atomid' must be an integer greater than or equal to zero!")            
        else:
            raise AtomCollectionError("Invalid index!")
    
    @dispatch(int,Atom)
    def __setitem__(self,idx,atom):
        '''
        Assigns an Atom object to the corresponding index in the structure's 
        atom list.

        Parameters
        ----------
        idx : integer
            Index in the structure's atom list.
        atom : Atom
            Atom object to be assigned to the corresponding index in the structure's 
            atom list.
        '''
        maxidx=len(self._atoms)
        
        if idx>=-maxidx and idx<maxidx:
            if atom._belongs_to is None:
                self._atoms[idx]._belongs_to=None
                
                for other_atom in self._atoms:
                    if self._atoms[idx]._symbol==other_atom._symbol:
                        break
                else:
                    self._species.remove(self._atoms[idx]._symbol)
                
                self._atoms[idx]=atom
                atom._belongs_to=self
                
                if not atom._symbol in self._species:
                    self._species.append(atom._symbol)
                        
                self._update()
            else:
                raise AtomCollectionError("Atom already belongs to another structure!")
        else:
            raise AtomCollectionError("Invalid index!")
    
    def __len__(self):
        '''
        Returns the number of Atom objects (active or not) in the structure.

        Returns
        -------
        integer
            Number of Atom objects (active or not) in the structure.
        '''
        return len(self._atoms)
    
    def __iter__(self):
        '''
        Returns an iterator that allows the user to iterate over the Atom objects 
        in the structure's list of atoms.

        Returns
        -------
        iterator
            An iterator that allows the user to iterate over the Atom objects in 
            the structure's list of atoms.
        '''
        return iter(self._atoms)
    
    @classmethod
    def _append(cls,obj):
        '''
        Appends an AtomCollection object to the list of AtomCollection objects 
        created so far.

        Parameters
        ----------
        obj : AtomCollection object
            AtomCollection object to be added to the list of AtomCollection objects 
            created so far.
        '''        
        if cls._instances is None:
            cls._currid=0
            cls._instances=[]
        
        obj._id=cls._currid
        cls._instances.append(obj)
        cls._currid+=1
        
    @abstractmethod    
    def _update(self):
        '''
        To be implemented in a derived atomic structure, this method is expected 
        to update specific attributes of that structure when atomic coordinates 
        or cell shape or cell size change.
        '''
        pass

    '''
    Properties
    ----------
   label : string
       Label used to identify the atomic structure.
    belongs_to : Hybrid-derived object, readonly
        Hybrid atomic structure to which the atom collection object belongs.
    unit_cell : Python tuple, readonly
        Tuple of Atom objects in the atomic structure's unit cell.
    active_atoms : Python tuple, readonly
        Tuple of atoms that are active and contained in the current supercell, 
        the boundaries of which are determined by the n, m and l attributes.
    species : Python tuple, readonly
        Tuple of atomic species in the atomic structure.
    a0 : float
        Lattice constant.
    lattice_vectors : Numpy array
        Lattice vectors.
    origin : Numpy array
        Position with minimum values of X, Y, and Z coordinates of the atomic 
        structure.
    n, m, l : integer, readonly
        Number of repetitions of the unit cell along the first, second and third 
        lattice vectors, respectively.
    ID : integer, readonly
        Unique atomic structure identifier.
    '''
    @property
    def label(self):
        return self._label
    
    @label.setter
    def label(self,val):
        if isinstance(val,str) and len(val)>0:
            self._label=val
        else:
            raise AtomCollectionError("The label must be a non-empty string!")
    
    @property
    def belongs_to(self):
        return self._belongs_to

    @property
    def unit_cell(self):
        return tuple(self._ucell)
    
    @property
    def active_atoms(self):
        return tuple(atom for atom in self._atoms if atom._active and 
                     self._loc[atom._id][0]<=self._n and 
                     self._loc[atom._id][1]<=self._m and 
                     self._loc[atom._id][2]<=self._l)
    
    @property
    def species(self):
        return tuple(self._species)
    
    @property
    def a0(self):
        return self._a0
    
    @a0.setter
    def a0(self,val):      
        if isinstance(val,(float,int)) and val>0.0:
            print("WARNING: Changing the lattice parameter rescales the atomic coordinates accordingly!")
            
            old,self._a0=self._a0,val
            fac=val/old
            
            for atom in self._atoms:
                atom._x*=fac
            
            self._update()
        else:
            raise AtomCollectionError("The lattice constant must be a positive number greater than zero!")

    @property
    def lattice_vectors(self):
        return self._latvec
    
    @lattice_vectors.setter
    def lattice_vectors(self,val):
        if isinstance(val,(list,tuple)):
            val=array(val)
        
        if isinstance(val,ndarray) and val.shape[0]==3 and val.size==9:
            self._latvec=val.astype(float)
        else:
            raise AtomCollectionError("Lattice vectors must be a Numpy array with three vectors!")

    @property
    def origin(self):
        return self._origin
    
    @origin.setter
    def origin(self,val):
        if isinstance(val,(list,tuple)):
            val=array(val)
            
        if isinstance(val,ndarray) and val.shape[0]==3:
            print("WARNING: Moving the origin displaces all atoms accordingly!")
            
            disp=val-self._origin
            self._origin=val.astype(float)
            
            self.displace(disp)
        else:
            raise AtomCollectionError("The new origin must be a Numpy array with three components!")

    @property
    def n(self):
        return self._n
    
    @property
    def m(self):
        return self._m
    
    @property
    def l(self):
        return self._l
        
    @property
    def ID(self):
        return self._id