from __future__ import print_function
from MolAdsPy.__species__ import Species
from MolAdsPy.__atom__ import Atom
from MolAdsPy.__exception__ import BasicException
from math import cos,sin,radians
from copy import deepcopy
from multipledispatch import dispatch
import os.path

class Molecule:
    __currmolid__=0   # ID available for a new molecule
    
    @dispatch(str)
    def __init__(self,moltype):
        '''
        __init__(moltype) -> class constructor.

        Parameters
        ----------
        moltype : string
            A label allowing you to identify the type of the molecule, for 
            instance, "benzene" or "C6H12".        
        '''
        if isinstance(moltype,str) and len(moltype)>0:
            self.__moltype__=moltype
        else:
            raise MoleculeError("'moltype' must be a string!")
                    
        self.__atoms__=[]                               # List of Atom objects that constitute the molecule
        self.__species__=[]                             # List of atomic species in the molecule
        self.__maxx__=self.__maxy__=self.__maxz__=None  # Maximum Cartesian coordinates
        self.__minx__=self.__miny__=self.__minz__=None  # Minimum Cartesian coordinates
        self.__anchors__={"com":[]}                     # Anchor points for translations and rotations
        self.__id__=Molecule.__currmolid__
        Molecule.__currmolid__+=1                

    @dispatch(str,str,str)    
    def __init__(self,moltype,filename,filetype):
        '''
        __init__(moltype,filename,filetype) -> class constructor.

        Parameters
        ----------
        moltype : string
            A label allowing you to identify the type of the molecule, for 
            instance, "benzene" or "C6H12".
        filename : string
            Name of an XYZ file containing the molecule structure.
        filetype : string
            Type of the file containing the molecule structure. For now, only 
            'XYZ' is accepted.

        Returns
        -------
        None.
        '''
        if isinstance(moltype,str) and len(moltype)>0:
            self.__moltype__=moltype
        else:
            raise MoleculeError("'moltype' must be a string!")
                    
        self.__atoms__=[]                               # List of Atom objects that constitute the molecule
        self.__species__=[]                             # List of atomic species in the molecule
        self.__maxx__=self.__maxy__=self.__maxz__=None  # Maximum Cartesian coordinates
        self.__minx__=self.__miny__=self.__minz__=None  # Minimum Cartesian coordinates
        self.__anchors__={"com":[]}                     # Anchor points for translations and rotations

        if isinstance(filetype,str) and filetype.lower()=='xyz':       
            if isinstance(filename,str) and len(filename)>0:
                self.__read_xyz__(filename)
            else:
                raise MoleculeError("'filename' must be a valid file name!")
        else:
            raise MoleculeError("'filetype' must be 'XYZ'!")
                                    
        self.__id__=Molecule.__currmolid__
        Molecule.__currmolid__+=1
        
    @dispatch(str,list)    
    def __init__(self,moltype,atomlist):
        '''
        __init__(moltype,atomlist) -> class constructor.

        Parameters
        ----------
        moltype : string
            A label allowing you to identify the type of the molecule, for 
            instance, "benzene" or "C6H12".
        atomlist : Python list
            List of Atom objects to be added to the molecule.

        Returns
        -------
        None.
        '''
        if isinstance(moltype,str) and len(moltype)>0:
            self.__moltype__=moltype
        else:
            raise MoleculeError("'moltype' must be a string!")
                    
        self.__atoms__=[]                               # List of Atom objects that constitute the molecule
        self.__species__=[]                             # List of atomic species in the molecule
        self.__maxx__=self.__maxy__=self.__maxz__=None  # Maximum Cartesian coordinates
        self.__minx__=self.__miny__=self.__minz__=None  # Minimum Cartesian coordinates
        self.__anchors__={"com":[]}                     # Anchor points for translations and rotations
                        
        if isinstance(atomlist,list):
            for atom in atomlist:
                self.add_atom(atom)
        else:
            raise MoleculeError("'atomlist' must be a non-empyt list!")
                                    
        self.__id__=Molecule.__currmolid__
        Molecule.__currmolid__+=1
        
    def add_atom(self,atom):
        '''
        add_atom(atom) -> adds an atom to the molecule.

        Parameters
        ----------
        atom : Atom object
            New atom added to the molecule.

        Returns
        -------
        None.
        '''
        if isinstance(atom,Atom):
            self.__atoms__.append(atom)
            
            if not atom.__symbol__ in self.__species__:
                self.__species__.append(atom.__symbol__)
            
            if len(self.__atoms__)==1:
                self.__maxx__=self.__minx__=atom.__x__[0]
                self.__maxy__=self.__miny__=atom.__x__[1]
                self.__maxz__=self.__minz__=atom.__x__[2]
            else:
                if atom.__x__[0]>self.__maxx__:
                    self.__maxx__=atom.__x__[0]
                elif atom.__x__[0]<self.__minx__:
                    self.__minx__=atom.__x__[0]
                
                if atom.__x__[1]>self.__maxy__:
                    self.__maxy__=atom.__x__[1]
                elif atom.__x__[1]<self.__miny__:
                    self.__miny__=atom.__x__[1]
                
                if atom.__x__[2]>self.__maxz__:
                    self.__maxz__=atom.__x__[2]
                elif atom.__x__[2]<self.__minz__:
                    self.__minz__=atom.__x__[2]
                    
            self.__calc_com__()
        else:
            raise MoleculeError("'atom' must be an Atom object!")

    @dispatch(int)            
    def remove_atom(self,atomid):
        '''
        remove_atom(atomid) -> removes an atom from the molecule.

        Parameters
        ----------
        atomid : integer
            ID of the atom to be removed.

        Returns
        -------
        None.
        '''     
        if isinstance(atomid,int) and atomid>=0 and atomid<Atom.__curratomid__:
            for atom in self.__atoms__:
                if atom.__id__==atomid:
                    self.__atoms__.remove(atom)
                    self.__calc_com__()
                    self.__get_extremes__()
                    
                    for other_atom in self.__atoms__:
                        if atom.__symbol__==other_atom.__symbol__:
                            break
                    else:
                        self.__species__.remove(atom.__symbol__)
                                        
                    break
            else:
                raise MoleculeError("Atom ID not found!")                
        else:
            raise MoleculeError("'atomid' must be an integer greater than or equal to zero!")
            
    @dispatch(Atom)            
    def remove_atom(self,atom):
        '''
        remove_atom(atom) -> removes an atom from the molecule.

        Parameters
        ----------
        atom : Atom object
            Atom to be removed.

        Returns
        -------
        None.
        '''     
        if isinstance(atom,Atom):
            if atom in self.__atoms__:
                self.__atoms__.remove(atom)
                self.__calc_com__()
                self.__get_extremes__()
                
                for other_atom in self.__atoms__:
                    if atom.__symbol__==other_atom.__symbol__:
                        break
                else:
                    self.__species__.remove(atom.__symbol__)                                    
            else:
                raise MoleculeError("Atom not found!")                
        else:
            raise MoleculeError("'atom' must be an Atom object!")

    def add_anchor(self,anchor,pos):
        '''
        add_anchor(anchor,pos) -> adds a new anchor point that can be used as 
            the reference point of the molecule for translations and rotations.

        Parameters
        ----------
        anchor : string
            Name of the anchor point.
        pos : Python list
            Cartesian coordinates of the anchor point.

        Returns
        -------
        None.
        '''
        if isinstance(anchor,str) and len(anchor)>0:
            if isinstance(pos,list) and len(pos)==3 and \
                isinstance(pos[0],(int,float)) and \
                isinstance(pos[1],(int,float)) and \
                isinstance(pos[2],(int,float)):
                self.__anchors__[anchor]=pos.copy()
            else:
                raise MoleculeError("'pos' must be a list with three components!")
        else:
            raise MoleculeError("'anchor' must be a non-empty string!")
            
    def remove_anchor(self,anchor):
        '''
        remove_anchor(anchor) -> removes an anchor point.

        Parameters
        ----------
        anchor : string
            Name of the anchor point.

        Returns
        -------
        None.
        '''
        if isinstance(anchor,str) and anchor in self.__anchors__.keys():
            self.__anchors__.pop(anchor)
        else:
            raise MoleculeError("'anchor' is not a valid anchor point!")
        
    def write_xyz(self,filename="molecule.xyz"):
        '''
        write_xyz(filename) -> saves the atomic coordinates of the molecule
            into an XYZ file.

        Parameters
        ----------
        filename : string, optional
            Name of the XYZ file. The default is "molecule.xyz".

        Returns
        -------
        None.
        '''
        with open(filename,"w") as f:
            f.write("%d\n\n" % (len(self.__atoms__)))
            
            for atom in self.__atoms__:
                f.write("%s %.6f %.6f %.6f\n" % (atom.__symbol__,atom.__x__[0],
                                                 atom.__x__[1],atom.__x__[2]))
                
    def write_pw_input(self,filename="molecule.in",vacuum=10.0,**kwargs):
        '''
        write_pw_input(filename,vacuum,**kwargs) -> creates a basic input file
            for geometry relaxation of the molecule using the pw.x code in the 
            Quantum Espresso package.

        Parameters
        ----------
        filename : string, optional
            Name of the input file. The default is "molecule.in".
        vacuum : float, optional
            Vacuum size separating the molecule from its images in all 
            directions. The default is 10.0 angstroms.
        **kwargs : Python dictionary.
            Dictionary containing key-value pairs that allow some customization 
            of the input file. At the moment, those are the keys accepted:
                ecutwfc: float, optional
                    Plane-wave energy cutoff. The default is 24.0 Ry.
                ecutrho: float, optional
                    Charge density cutoff. The default is 96.0 Ry.
                nspin: int, optional
                    Spin polarization. It can be either 1 (non-polarized) or 
                    2 (polarized, magnetization along z axis). The default is 1.
                starting_magnetization : Python list, optional
                    List of starting magnetic moments for the different atom 
                    types in the molecule. The default is 1.0 Bohr magneton 
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

        Returns
        -------
        None.
        '''
        inpstr="&CONTROL\n"
        inpstr+="    calculation='relax',\n"
        inpstr+="    restart_mode='from_scratch',\n"
        inpstr+="    pseudo_dir='.',\n"
        inpstr+="    prefix='"+os.path.splitext(os.path.basename
                                                (filename))[0]+"',\n"
        inpstr+="    nstep=500,\n/\n"
        inpstr+="&SYSTEM\n"
        inpstr+="    ibrav=0,\n"
        celldm=self.__maxx__-self.__minx__+vacuum
        inpstr+="    celldm(1)=%.4f,\n" % (celldm*1.8897259886)
        inpstr+="    nat=%d,\n" % (len(self.__atoms__))
        inpstr+="    ntyp=%d,\n" % (len(self.__species__))
        
        if "ecutwfc" in kwargs:
            inpstr+="    ecutwfc=%.4f,\n" % kwargs["ecutwfc"]
        else:
            inpstr+="    ecutwfc=24.00,\n"
        
        if "ecutrho" in kwargs:
            inpstr+="    ecutrho=%.4f,\n" % kwargs["ecutrho"]
        else:
            inpstr+="    ecutrho=96.00,\n"
            
        if "nspin" in kwargs and kwargs["nspin"]==2:
            inpstr+="    nspin=2,\n"
            inpstr+="    occupations='smearing',\n"
            inpstr+="    smearing='gaussian',\n"
            inpstr+="    degauss=0.02,\n"
                
            if "starting_magnetization" in kwargs and \
                len(kwargs["starting_magnetization"])>0:
                    for i in range(len(kwargs["starting_magnetization"])):
                        if i<len(self.__species__):                           
                            inpstr+="    starting_magnetization(%d)=%.2f,\n" \
                                % (i+1,kwargs["starting_magnetization"][i])
                        else:
                            break
            else:
                inpstr+="    starting_magnetization(1)=1.00,\n"
        else:
            inpstr+="    nspin=1,\n"
            
        if "input_dft" in kwargs and kwargs["input_dft"]=="vdw-df":
            inpstr+="    input_dft='vdw-df',\n"
        
        inpstr+="/\n"
        inpstr+="&ELECTRONS\n"
        
        if "mixing_mode" in kwargs and (kwargs["mixing_mode"]=="TF" or\
            kwargs["mixing_mode"]=="local-TF"):
                inpstr+="    mixing_mode='%s',\n" % kwargs["mixing_mode"]
        else:
            inpstr+="    mixing_mode='plain',\n"
        
        if "mixing_beta" in kwargs and kwargs["mixing_beta"]>0.0:
            inpstr+="    mixing_beta=%.4f,\n" % kwargs["mixing_beta"]
        else:
            inpstr+="    mixing_beta=0.7,\n"
        
        inpstr+="    electron_maxstep=200,\n/\n"
        inpstr+="&IONS\n/\n"
        inpstr+="CELL_PARAMETERS alat\n"
        inpstr+="1.0000   0.0000   0.0000\n"
        inpstr+="0.0000   %.4f   0.0000\n" % ((self.__maxy__-
                                               self.__miny__+vacuum)
                                              /celldm)
        inpstr+="0.0000   0.0000   %.4f\n" % ((self.__maxz__-
                                               self.__minz__+vacuum)
                                              /celldm)
        inpstr+="K_POINTS gamma\n"
        inpstr+="ATOMIC_SPECIES\n"
        
        for symbol in self.__species__:
            for species in Species.__species__:
                if species.__symbol__==symbol:
                    if species.__pseudopotential__ is None or\
                        len(species.__pseudopotential__)==0:
                        pseudopotential=symbol+".UPF"
                    else:
                        pseudopotential=species.__pseudopotential__
                    
                    inpstr+="%s %.4f %s\n" % (symbol,species.__atomicmass__,
                                              pseudopotential)
                    
                    break
                
        inpstr+="ATOMIC_POSITIONS angstrom\n"
        
        for atom in self.__atoms__:                    
            inpstr+="%s %.6f %.6f %.6f %d %d %d\n" % (atom.__symbol__,
                                                      atom.__x__[0],
                                                      atom.__x__[1],
                                                      atom.__x__[2],
                                                      int(not(atom.__fixed__[0])),
                                                      int(not(atom.__fixed__[1])),
                                                      int(not(atom.__fixed__[2])))
            
        with open(filename,"w") as f:
            f.write(inpstr)
        
    def displace(self,disp):
        '''
        displace(disp) -> rigidly displaces the molecule.

        Parameters
        ----------
        disp : Python list
            Displacement vector.

        Returns
        -------
        None.
        '''
        if isinstance(disp,list) and len(disp)==3 and \
            isinstance(disp[0],(float,int)) and \
            isinstance(disp[1],(float,int)) and \
            isinstance(disp[2],(float,int)):
            for atom in self.__atoms__:
                atom.__x__[0]+=disp[0]
                atom.__x__[1]+=disp[1]
                atom.__x__[2]+=disp[2]
                
            for key in self.__anchors__.keys():
                self.__anchors__[key][0]+=disp[0]
                self.__anchors__[key][1]+=disp[1]
                self.__anchors__[key][2]+=disp[2]

            self.__maxx__+=disp[0]
            self.__minx__+=disp[0]
            self.__maxy__+=disp[1]
            self.__miny__+=disp[1]
            self.__maxz__+=disp[2]
            self.__minz__+=disp[2]
        else:
            raise MoleculeError("'disp' must be a list with three components!")
            
    def move_to(self,x=[0.0,0.0,0.0],anchor="com"):
        '''
        move_to(x,anchor) -> rigidly moves the molecule such that the anchor 
        point is located at 'x'.

        Parameters
        ----------
        x : Python list, optional
            New position of the molecule's anchor point. The default is 
            [0.0,0.0,0.0].
        anchor : string, optional
            Anchor point used as reference for translations. The default is 
            'com', which means the molecule's center of mass.

        Returns
        -------
        None.
        '''
        disp=[0.0,0.0,0.0]
        
        if isinstance(anchor,str) and anchor in self.__anchors__.keys():
            anchorpos=[self.__anchors__[anchor][0],
                       self.__anchors__[anchor][1],
                       self.__anchors__[anchor][2]]
        else:
            raise MoleculeError("'anchor' is not a valid anchor point!")
        
        if isinstance(x,list) and len(x)==3 and \
            isinstance(x[0],(float,int)) and \
            isinstance(x[1],(float,int)) and \
            isinstance(x[2],(float,int)):
            disp[0]=x[0]-anchorpos[0]
            disp[1]=x[1]-anchorpos[1]
            disp[2]=x[2]-anchorpos[2]
        
            self.displace(disp)
        else:
            raise MoleculeError("'disp' must be a list with three components!")

    def rotate(self,theta,phi,psi,anchor="com"):
        '''
        rotate(theta,phi,psi) -> rotates the molecule around an anchor point.

        Parameters
        ----------
        theta : float, optional
            First Euler angle, in degrees.
        phi : float, optional
            Second Euler angle, in degrees.
        psi : float, optional
            Third Euler angle, in degrees.
        anchor : string, optional
            Anchor point around which the molecule will be rotated. The default 
            is 'com', which means the molecule's center of mass.

        Returns
        -------
        None.
        '''
        if isinstance(theta,(float,int)) and \
            isinstance(phi,(float,int)) and \
            isinstance(psi,(float,int)):       
            theta=radians(theta)
            phi=radians(phi)
            psi=radians(psi)
        else:
            raise MoleculeError("Euler angles must be provided as numbers!")
            
        if isinstance(anchor,str) and anchor in self.__anchors__.keys():
            anchorpos=[self.__anchors__[anchor][0],
                       self.__anchors__[anchor][1],
                       self.__anchors__[anchor][2]]
        else:
            raise MoleculeError("'anchor' is not a valid anchor point!")
            
        a11=cos(theta)*cos(psi)
        a12=-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)
        a13=sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)
        a21=cos(theta)*sin(psi)
        a22=cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)
        a23=-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)
        a31=-sin(theta)
        a32=sin(phi)*cos(theta)
        a33=cos(phi)*cos(theta)
        
        for atom in self.__atoms__:
            x1=atom.__x__[0]-anchorpos[0]
            y1=atom.__x__[1]-anchorpos[1]
            z1=atom.__x__[2]-anchorpos[2]
            x2=x1*a11+y1*a12+z1*a13
            y2=x1*a21+y1*a22+z1*a23
            z2=x1*a31+y1*a32+z1*a33
            atom.__x__[0]=x2+anchorpos[0]
            atom.__x__[1]=y2+anchorpos[1]
            atom.__x__[2]=z2+anchorpos[2]
            
        for key in self.__anchors__.keys():
            ax1=self.__anchors__[key][0]-anchorpos[0]
            ay1=self.__anchors__[key][1]-anchorpos[1]
            az1=self.__anchors__[key][2]-anchorpos[2]
            ax2=ax1*a11+ay1*a12+az1*a13
            ay2=ax1*a21+ay1*a22+az1*a23
            az2=ax1*a31+ay1*a32+az1*a33
            self.__anchors__[key][0]=ax2+anchorpos[0]
            self.__anchors__[key][1]=ay2+anchorpos[1]
            self.__anchors__[key][2]=az2+anchorpos[2]
        
        self.__calc_com__()
        self.__get_extremes__()
                    
    def center(self):
        '''
        center() -> displaces the molecule such that its center of mass will be 
        located at [0.0,0.0,0.0].

        Returns
        -------
        None.
        '''
        for atom in self.__atoms__:
            atom.__x__[0]-=self.__anchors__["com"][0]
            atom.__x__[1]-=self.__anchors__["com"][1]
            atom.__x__[2]-=self.__anchors__["com"][2]
        
        self.__maxx__-=self.__anchors__["com"][0]
        self.__minx__-=self.__anchors__["com"][0]
        self.__maxy__-=self.__anchors__["com"][1]
        self.__minx__-=self.__anchors__["com"][1]
        self.__maxz__-=self.__anchors__["com"][2]
        self.__minz__-=self.__anchors__["com"][2]
        self.__anchors__["com"]=[0.0,0.0,0.0]
        
    def copy(self):
        '''
        copy() -> returns a copy of the current Molecule object.

        Returns
        -------
        Molecule object
            Copy of the molecule.
        '''
        newmol=deepcopy(self)
        newmol.__id__=Molecule.__currmolid__
        Molecule.__currmolid__+=1
        
        for atom in newmol.atoms:
            atom.__id__=Atom.__curratomid__
            Atom.__curratomid__+=1
        
        return newmol
        
    @dispatch(int)
    def displace_atom(self,atomid,disp):
        '''
        displace_atom(atomid,disp) -> displaces an atom in the molecule from 
            its current position.

        Parameters
        ----------
        atomid : int
            ID of the atom that will be moved to a new position.
        disp : Python list
            Displacement vector.

        Returns
        -------
        None.
        '''
        if isinstance(atomid,int) and atomid>=0 and atomid<Atom.__curratomid__:
            for atom in self.__atoms__:
                if atom.__id__==atomid:
                    atom.displace(disp)
                    self.__calc_com__()
                    self.__get_extremes__()
                    
                    break
            else:
                raise MoleculeError("Atom ID not found!")               
        else:
            raise MoleculeError("'atomid' must be an integer greater than or equal to zero!")
            
    @dispatch(Atom)
    def displace_atom(self,atom,disp):
        '''
        displace_atom(atom,disp) -> displaces an atom in the molecule from 
            its current position.

        Parameters
        ----------
        atom : Atom object
            Atom that will be moved to a new position.
        disp : Python list
            Displacement vector.

        Returns
        -------
        None.
        '''
        if isinstance(atom,Atom):
            if atom in self.__atoms__:
                atom.displace(disp)
                self.__calc_com__()
                self.__get_extremes__()                    
            else:
                raise MoleculeError("Atom not found!")               
        else:
            raise MoleculeError("'atom' must be an Atom object!")
                    
    def list_atoms(self):
        '''
        list_atoms() -> prints a list of the atoms in the molecule.

        Returns
        -------
        None.
        '''
        print("AtomID, Symbol, Element, X, Y, Z\n------")
        
        for atom in self.__atoms__:
            print("%d, %s, %s, %f, %f, %f" % (atom.__id__,atom.__symbol__,
                                              atom.__element__,atom.__x__[0],
                                              atom.__x__[1],atom.__x__[2]))

    def list_anchors(self):
        '''
        list_anchors() -> lists the anchor points in the molecule.

        Returns
        -------
        None.
        '''
        print("Anchor, X, Y, Z\n------")
        
        for key in self.__anchors__.keys():
            print("%s, %f, %f, %f" % (key,self.__anchors__[key][0],
                                      self.__anchors__[key][1],
                                      self.__anchors__[key][2]))  

    def __read_xyz__(self,filename):
        '''
        __read_xyz__(filename) -> reads the structure of the molecule from a file
            in XYZ format.

        Parameters
        ----------
        filename : string
            Name of the XYZ file.

        Returns
        -------
        None.
        '''        
        with open(filename,"r") as f:
            lines=f.readlines()

        count=0

        for line in lines:
            if count<2:
                count+=1
                
                continue
            
            l=line.split()
            
            if len(l)>=4:
                species=None
            
                for obj in Species.__species__:
                    if l[0]==obj.__symbol__:
                        species=obj
                        
                        break
                else:
                    try:
                        species=Species(l[0])
                    except:                    
                        print("WARNING: Species '%s' and atoms of this species must be added manually!" % l[0])
    
                if species is not None:                
                    x=[float(l[1]),float(l[2]),float(l[3])]
                    
                    self.__atoms__.append(Atom(species,x))
                    
                    if not species.__symbol__ in self.__species__:
                        self.__species__.append(species.__symbol__)
                        
        print("INFO: %d Atom objects have just been created from file '%s'!" 
              % (len(self.__atoms__),filename))

        self.__calc_com__()
        self.__get_extremes__()

    def __calc_com__(self):
        '''
        __calc_com__() -> sets the molecule's center of mass.

        Returns
        -------
        None.
        '''
        if len(self.__atoms__)>0:
            X=Y=Z=0.0
            totmass=0.0
        
            for atom in self.__atoms__:
                totmass=totmass+atom.__atomicmass__
                X+=atom.__atomicmass__*atom.__x__[0]
                Y+=atom.__atomicmass__*atom.__x__[1]
                Z+=atom.__atomicmass__*atom.__x__[2]

            self.__anchors__["com"]=[X/totmass,Y/totmass,Z/totmass]
            
    def __get_extremes__(self):
        '''
        __get_extremes__() -> gets the extreme coordinates of the molecule in 
        all directions.

        Returns
        -------
        None.
        '''
        natoms=len(self.__atoms__)
        
        if natoms==1:
            self.__maxx__=self.__minx__=self.__atoms__[0].__x__[0]
            self.__maxy__=self.__miny__=self.__atoms__[0].__x__[1]
            self.__maxz__=self.__minz__=self.__atoms__[0].__x__[2]
        elif natoms>1:
            self.__maxx__=max([atom.__x__[0] for atom in self.__atoms__])
            self.__minx__=min([atom.__x__[0] for atom in self.__atoms__])
            self.__maxy__=max([atom.__x__[1] for atom in self.__atoms__])
            self.__miny__=min([atom.__x__[1] for atom in self.__atoms__])
            self.__maxz__=max([atom.__x__[2] for atom in self.__atoms__])
            self.__minz__=min([atom.__x__[2] for atom in self.__atoms__])
            
    '''
    Properties
    ----------
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
    '''
    @property
    def moleculetype(self):
        return self.__moltype__
    
    @property
    def atoms(self):
        return self.__atoms__
    
    @property
    def anchors(self):
        return self.__anchors__
    
    @property
    def atomtypes(self):
        return self.__species__
    
    @property
    def centerofmass(self):
        return self.__anchors__["com"]
    
    @centerofmass.setter
    def centerofmass(self,val):
        if isinstance(val,list) and len(val)==3 and \
            isinstance(val[0],(int,float)) and \
            isinstance(val[1],(int,float)) and \
            isinstance(val[2],(int,float)):
            self.move_to(val,"com")
        else:
            raise MoleculeError("Center of mass position must be a list with three components!")            
    
    @property
    def maxx(self):
        return self.__maxx__
    
    @property
    def maxy(self):
        return self.__maxy__
    
    @property
    def maxz(self):
        return self.__maxz__
    
    @property
    def minx(self):
        return self.__minx__
    
    @property
    def miny(self):
        return self.__miny__
    
    @property
    def minz(self):
        return self.__minz__
    
    @property
    def ID(self):
        return self.__id__
    
class MoleculeError(BasicException):
    pass