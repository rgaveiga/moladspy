from __future__ import print_function
from MolAdsPy.__species__ import Species
from MolAdsPy.__atom__ import Atom
from MolAdsPy.__exception__ import BasicException
from math import sqrt
from copy import deepcopy
from multipledispatch import dispatch
import os.path

class Slab:
    __currslabid__=0    # ID available for a new slab
    
    @dispatch(str)
    def __init__(self,slabname):
        '''
        __init__(slabname) -> class constructor.

        Parameters
        ----------
        slabname : string
            A label allowing you to identify the slab, for instance, "graphene".
        '''
        if isinstance(slabname,str) and len(slabname)>0:
            self.__slabname__=slabname
        else:
            raise SlabError("'slabname' must be a string!")
            
        self.__atoms__=[]       # List of Atom objects that constitute the slab
        self.__species__=[]     # List of atomic species in the slab
        self.__ads_sites__={}   # Dictionary of adsorption sites
        self.__top__=None       # Maximum Z coordinate of the slab
        self.__bottom__=None    # Miniumum Z coordinate of the slab
        self.__origin__=[]      # Origin of the slab
        
        self.origin=[0.0,0.0,0.0]              
        self.__id__=Slab.__currslabid__
        Slab.__currslabid__+=1
    
    @dispatch(str,str,str)
    def __init__(self,slabname,filename,filetype):
        '''
        __init__(slabname,filename,filetype) -> class constructor.

        Parameters
        ----------
        slabname : string
            A label allowing you to identify the slab, for instance, "graphene".
        filename : string
            Name of the file containing the slab structure.
        filetype : string
            Type of the file containing the slab structure. For now, only 'XYZ' 
            is accepted.

        Returns
        -------
        None.
        '''
        if isinstance(slabname,str) and len(slabname)>0:
            self.__slabname__=slabname
        else:
            raise SlabError("'slabname' must be a string!")
            
        self.__atoms__=[]       # List of Atom objects that constitute the slab
        self.__species__=[]     # List of atomic species in the slab
        self.__ads_sites__={}   # Dictionary of adsorption sites
        self.__top__=None       # Maximum Z coordinate of the slab
        self.__bottom__=None    # Miniumum Z coordinate of the slab
        self.__origin__=[]      # Origin of the slab

        if isinstance(filetype,str) and filetype.lower()=='xyz':        
            if isinstance(filename,str) and len(filename)>0:
                self.__read_xyz__(filename)
            else:
                raise SlabError("'filename' must be a valid file name!")
        else:
            raise SlabError("'filetype' must be 'XYZ'!")
            
        self.origin=[0.0,0.0,0.0]              
        self.__id__=Slab.__currslabid__
        Slab.__currslabid__+=1

    @dispatch(str,list,list)
    def __init__(self,slabname,atomlist,lattice_vectors):
        '''
        __init__(slabname,atomlist,lattice_vectors) -> class constructor.

        Parameters
        ----------
        slabname : string
            A label allowing you to identify the slab, for instance, "graphene".
        atomlist : Python list
            List of Atom objects to be added to the slab.
        lattice_vectors : Python list
            List of slab lattice vectors in the XY plane.

        Returns
        -------
        None.
        '''
        if isinstance(slabname,str) and len(slabname)>0:
            self.__slabname__=slabname
        else:
            raise SlabError("'slabname' must be a string!")
                        
        self.__atoms__=[]       # List of Atom objects that constitute the slab
        self.__species__=[]     # List of atomic species in the slab
        self.__ads_sites__={}   # Dictionary of adsorption sites
        self.__top__=None       # Maximum Z coordinate of the slab
        self.__bottom__=None    # Miniumum Z coordinate of the slab
        self.__origin__=[]      # Origin of the slab
        
        if isinstance(lattice_vectors,list) and len(lattice_vectors)==2:
            if isinstance(lattice_vectors[0],list) and \
                isinstance(lattice_vectors[1],list) and \
                len(lattice_vectors[0])==len(lattice_vectors[1])==2:
                self.__lattice_vectors__=lattice_vectors.copy()
            else:
                raise SlabError("Each component of 'lattice_vectors' must be a list with two elements!")
        else:
            raise SlabError("'lattice_vectors' must be a list with two elements!")
            
        if isinstance(atomlist,list) and len(atomlist)>0:
            for atom in atomlist:
                self.add_atom(atom)
        else:
            raise SlabError("'atomlist' must be a non-empyt list!")
            
        self.origin=[0.0,0.0,0.0]                
        self.__id__=Slab.__currslabid__
        Slab.__currslabid__+=1
                
    def add_atom(self,atom):
        '''
        add_atom(atom) -> adds an atom to the slab.

        Parameters
        ----------
        atom : Atom object
            New atom added to the slab.

        Returns
        -------
        None.
        '''
        if isinstance(atom,Atom):
            self.__atoms__.append(atom)
            
            if not atom.__symbol__ in self.__species__:
                self.__species__.append(atom.__symbol__)
            
            if len(self.__atoms__)==1:
                self.__top__=self.__bottom__=atom.__x__[2]
                self.__origin__=[atom.__x__[0],atom.__x__[1],atom.__x__[2]]
            else:
                if not (len(self.__origin__)==3 and \
                        isinstance(self.__origin__[0],(float,int)) and \
                        isinstance(self.__origin__[1],(float,int)) and \
                        isinstance(self.__origin__[2],(float,int))):
                    self.__get_origin__()
                    
                if self.__top__ is None or self.__bottom__ is None:
                    self.__get_z__()
                    
                if atom.__x__[0]<self.__origin__[0]:
                    self.__origin__[0]=atom.__x__[0]
                    
                if atom.__x__[1]<self.__origin__[1]:
                    self.__origin__[1]=atom.__x__[1]
                
                if atom.__x__[2]>self.__top__:
                    self.__top__=atom.__x__[2]
                elif atom.__x__[2]<self.__bottom__:
                    self.__bottom__=self.__origin__[2]=atom.__x__[2]
        else:
            raise SlabError("'atom' must be an Atom object!")
    
    @dispatch(int)    
    def remove_atom(self,atomid):
        '''
        remove_atom(atomid) -> removes an atom from the slab.

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
                    self.__get_z__()
                    self.__get_origin__()
                                        
                    for other_atom in self.__atoms__:
                        if atom.__symbol__==other_atom.__symbol__:
                            break
                    else:
                        self.__species__.remove(atom.__symbol__)
                    
                    break
            else:
                raise SlabError("Atom ID not found!")
        else:
            raise SlabError("'atomid' must be an integer greater than or equal to zero!")

    @dispatch(Atom)    
    def remove_atom(self,atom):
        '''
        remove_atom(atom) -> removes an atom from the slab.

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
                self.__get_z__()
                self.__get_origin__()
                                    
                for other_atom in self.__atoms__:
                    if atom.__symbol__==other_atom.__symbol__:
                        break
                else:
                    self.__species__.remove(atom.__symbol__)                
            else:
                raise SlabError("Atom not found!")
        else:
            raise SlabError("'atomid' must be an Atom object!")

    def add_adsorption_site(self,ads_site,pos):
        '''
        add_adsorption_site(ads_site,pos) -> adds an adsorption site on the 
            slab.

        Parameters
        ----------
        ads_site : string
            Label of the adsorption site.
        pos : Python list
            Coordinates of the adsorption site on the XY plane.

        Returns
        -------
        None.
        '''
        ads_site=str(ads_site)
                
        if isinstance(pos,list) and len(pos)==2:
            if isinstance(pos[0],(int,float)) and \
                isinstance(pos[1],(int,float)):
                    
                if not ads_site in self.__ads_sites__.keys():
                    self.__ads_sites__[ads_site]=[]
                    
                self.__ads_sites__[ads_site].append(pos)
            else:                
                raise SlabError("'pos' must contain XY numeric coordinates!")                
        else:
            raise SlabError("'pos' must be a list with two elements!")
                            
    def remove_adsorption_site(self,ads_site,index):
        '''
        remove_adsorption_site(ads_site,index) -> removes an adsorption site.

        Parameters
        ----------
        ads_site : string
            Label of the adsorption site.
        index : integer
            Index of the adsorption site.

        Returns
        -------
        None.
        '''
        ads_site=str(ads_site)
        
        if not ads_site in self.__ads_sites__.keys():
            raise SlabError("'ads_site' is not a valid label for an adsorption site!")
                    
        if isinstance(index,int) and index>=0 and \
            index<len(self.__ads_sites__[ads_site]):
            self.__ads_sites__[ads_site].pop(index)
            
            if len(self.__ads_sites__[ads_site])==0:
                self.__ads_sites__.pop(ads_site)
        else:
            raise SlabError("'index' must be an integer greater than or equal to zero!")
                    
    def replicate(self,n,m,replicate_ads_sites=True):
        '''
        replicate(n,m,replicate_ads_sites) -> replicates the slab in the 
            XY plane.

        Parameters
        ----------
        n : integer
            Number of replications along the first lattice vector.
        m : integer
            Number of replications along the second lattice vector.
        replicate_ads_sites : logical, optional
            Whether also replicate the adsorption sites or not. The default is 
            True.

        Returns
        -------
        None.
        '''
        if not (isinstance(n,int) and isinstance(m,int) and n>0 and m>0):
            raise SlabError("Number of repetitions must be integers greater than zero!")            
        elif n==1 and m==1:
            raise SlabError("Supercell is equal to the unit cell; nothing to do!")
                    
        ucell=deepcopy(self.__atoms__)
        
        if replicate_ads_sites and bool(self.__ads_sites__):
            ucell_ads_sites=deepcopy(self.__ads_sites__)
        else:
            ucell_ads_sites={}
        
        for i in range(n):
            for j in range(m):
                if i==0 and j==0:
                    continue
                
                for atom in ucell:
                    newatom=atom.copy()
                    newatom.__x__[0]+=i*self.__lattice_vectors__[0][0]+\
                        j*self.__lattice_vectors__[1][0]
                    newatom.__x__[1]+=i*self.__lattice_vectors__[0][1]+\
                        j*self.__lattice_vectors__[1][1]

                    self.add_atom(newatom)
                    
                if bool(ucell_ads_sites):                    
                    for key,sitelist in ucell_ads_sites.items():
                        for site in sitelist:
                            newsite=site.copy()
                            newsite[0]+=i*self.__lattice_vectors__[0][0]+\
                                j*self.__lattice_vectors__[1][0]
                            newsite[1]+=i*self.__lattice_vectors__[0][1]+\
                                j*self.__lattice_vectors__[1][1]
                        
                            self.add_adsorption_site(key,newsite)
                    
        self.__lattice_vectors__[0][0]*=n
        self.__lattice_vectors__[0][1]*=n
        self.__lattice_vectors__[1][0]*=m
        self.__lattice_vectors__[1][1]*=m
        
    def copy(self):
        '''
        copy() -> returns a copy of the current Slab object.

        Returns
        -------
        Slab object
            Copy of the slab.
        '''
        newslab=deepcopy(self)
        newslab.__id__=Slab.__currslabid__
        Slab.__currslabid__+=1
        
        for atom in newslab.atoms:
            atom.__id__=Atom.__curratomid__
            Atom.__curratomid__+=1
        
        return newslab
                                    
    def write_xyz(self,filename="slab.xyz"):
        '''
        write_xyz(filename) -> saves the atomic coordinates of the slab into
            an XYZ file.

        Parameters
        ----------
        filename : string, optional
            Name of the XYZ file. The default is "slab.xyz".

        Returns
        -------
        None.
        '''
        dimz=self.__top__-self.__bottom__+20.0
        
        with open(filename,"w") as f:
            f.write("%d\n" % (len(self.__atoms__)))
            f.write('Lattice="%.6f %.6f 0.0 %.6f %.6f 0.0 0.0 0.0 %.6f" Properties=species:S:1:pos:R:3\n' 
                    % (self.__lattice_vectors__[0][0],self.__lattice_vectors__[0][1],
                       self.__lattice_vectors__[1][0],self.__lattice_vectors__[1][1],
                       dimz))
            
            for atom in self.__atoms__:
                f.write("%s %.6f %.6f %.6f\n" % (atom.__symbol__,atom.__x__[0],
                                                 atom.__x__[1],atom.__x__[2]))
        
    def write_pw_input(self,filename="slab.in",vacuum=20.0,**kwargs):
        '''
        write_pw_input(filename,vacuum,**kwargs) -> creates a basic input file
            for geometry relaxation of the slab using the pw.x code in the 
            Quantum Espresso package.

        Parameters
        ----------
        filename : string, optional
            Name of the input file. The default is "slab.in".
        vacuum : float, optional
            Vacuum size separating the top of the slab from the bottom of its
            image. The default is 20.0 angstroms.
        **kwargs : Python dictionary.
            Dictionary containing key-value pairs that allow some customization 
            of the input file. At the moment, those are the keys accepted:
                ecutwfc: float, optional
                    Plane-wave energy cutoff. The default is 24.0 Ry.
                ecutrho: float, optional
                    Charge density cutoff. The default is 96.0 Ry.
                nspin: integer, optional
                    Spin polarization. It can be either 1 (non-polarized) or 
                    2 (polarized, magnetization along z axis). The default is 1.
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
        celldm=sqrt(pow((self.__lattice_vectors__[0][0]),2)+
                    pow((self.__lattice_vectors__[0][1]),2))
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
        inpstr+="%.4f   %.4f   0.0000\n" % (self.__lattice_vectors__[0][0]
                                            /celldm,
                                            self.__lattice_vectors__[0][1]
                                            /celldm)
        inpstr+="%.4f   %.4f   0.0000\n" % (self.__lattice_vectors__[1][0]
                                            /celldm,
                                            self.__lattice_vectors__[1][1]
                                            /celldm)
        inpstr+="0.0000   0.0000   %.4f\n" % ((self.__top__-self.__bottom__
                                               +vacuum)/celldm)

        inpstr+="K_POINTS automatic\n"

        if "kvec" in kwargs and len(kwargs["kvec"])==6:
            inpstr+="%d %d %d %d %d %d\n" % (kwargs["kvec"][0],
                                             kwargs["kvec"][1],
                                             kwargs["kvec"][2],
                                             kwargs["kvec"][3],
                                             kwargs["kvec"][4],
                                             kwargs["kvec"][5])
        else:
            inpstr+="1 1 1 0 0 0\n"
        
        inpstr+="ATOMIC_SPECIES\n"
        
        for symbol in self.__species__:
            for atomtype in Species.__species__:
                if atomtype.__symbol__==symbol:
                    if atomtype.__pseudopotential__ is None or\
                        len(atomtype.__pseudopotential__)==0:
                        pseudopotential=symbol+".UPF"
                    else:
                        pseudopotential=atomtype.__pseudopotential__
                    
                    inpstr+="%s %.4f %s\n" % (symbol,atomtype.__atomicmass__,
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

    @dispatch(int)       
    def displace_atom(self,atomid,disp):
        '''
        displace_atom(atomid,x) -> displaces an atom in the slab from its 
            current position.

        Parameters
        ----------
        atomid : integer
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
                    self.__get_z__()
                    self.__get_origin__()
                    
                    break
            else:
                raise SlabError("Atom ID not found!")                
        else:
            raise SlabError("'atomid' must be an integer greater than or equal to zero!")
            
    @dispatch(Atom)       
    def displace_atom(self,atom,disp):
        '''
        displace_atom(atom,x) -> displaces an atom in the slab from its current 
            position.

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
                self.__get_z__()
                self.__get_origin__()                
            else:
                raise SlabError("Atom not found!")                
        else:
            raise SlabError("'atom' must be an Atom object!")
                    
    def list_atoms(self):
        '''
        list_atoms() -> prints a list of the atoms in the slab.

        Returns
        -------
        None.
        '''
        print("AtomID, Symbol, Element, X, Y, Z\n------")
        
        for atom in self.__atoms__:
            print("%d, %s, %s, %f, %f, %f" % (atom.__id__,atom.__symbol__,
                                              atom.__element__,atom.__x__[0],
                                              atom.__x__[1],atom.__x__[2]))
            
    def list_adsorption_sites(self):
        '''
        list_adsorption_sites() -> prints a list of the adsorption sites in 
            the slab.

        Returns
        -------
        None.
        '''
        print("Label, Index, X, Y\n------")
        
        for key in self.__ads_sites__.keys():
            for i in range(len(self.__ads_sites__[key])):
                print("%s, %d, %f, %f" % (key,i,self.__ads_sites__[key][i][0],
                                          self.__ads_sites__[key][i][1]))

    def __read_xyz__(self,filename):
        '''
        __read_xyz__(filename) -> reads the structure of the slab from a file
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
                else:
                    self.__lattice_vectors__=[[0.0,0.0],[0.0,0.0]]
                    
                    print("WARNING: Lattice vectors not found in the XYZ file!")
                
            l=line.split()
            
            if count==1 and extxyz>=0:
                self.__lattice_vectors__=[[float(l[0]),float(l[1])],
                                          [float(l[3]),float(l[4])]]
                count=2
                
                continue
        
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

        self.__get_origin__()
        self.__get_z__()
            
    def __get_origin__(self):
        '''
        __get_origin__() -> assigns the origin of the slab.
        
        Returns
        -------
        None.
        '''        
        if len(self.__atoms__)==0:
            self.__origin__=[0.0,0.0,0.0]
        elif len(self.__atoms__)==1:
            self.__origin__=[self.__atoms__[0].__x__[0],
                             self.__atoms__[0].__x__[1],
                             self.__atoms__[0].__x__[2]]
        else:
            self.__origin__=[min([atom.__x__[0] for atom in self.__atoms__]),
                             min([atom.__x__[1] for atom in self.__atoms__]),
                             min([atom.__x__[2] for atom in self.__atoms__])]

    def __get_z__(self):
        '''
        __get_z__() -> assigns the minimum and maximum Z coordinates of the 
        slab, which corresponds to minimum and maximum positions of the two 
        surfaces.

        Returns
        -------
        None.
        '''       
        if len(self.__atoms__)==0:
            self.__top__=self.__bottom__=None
        elif len(self.__atoms__)==1:
            self.__top__=self.__bottom__=self.__atoms__[0].__x__[2]
        else:
            self.__top__=max([atom.__x__[2] for atom in self.__atoms__])
            self.__bottom__=min([atom.__x__[2] for atom in self.__atoms__])
            
    '''
    Properties
    ----------
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
    '''
    @property
    def label(self):
        return self.__slabname__
    
    @property
    def atoms(self):
        return self.__atoms__
    
    @property
    def atomtypes(self):
        return self.__species__
    
    @property
    def adsorptionsites(self):
        return self.__ads_sites__
    
    @property
    def top(self):
        return self.__top__
    
    @property
    def bottom(self):
        return self.__bottom__
    
    @property
    def origin(self):
        return self.__origin__
     
    @origin.setter    
    def origin(self,val):
        disp=[0.0,0.0,0.0]
        
        if not (isinstance(val,list) and len(val)==3 and \
            isinstance(val[0],(float,int)) and \
            isinstance(val[1],(float,int)) and \
            isinstance(val[2],(float,int))):
            raise SlabError("A list with three numeric components must be provided!")
                    
        if not (len(self.__origin__)==3 and \
                isinstance(self.__origin__[0],(float,int)) and \
                isinstance(self.__origin__[1],(float,int)) and \
                isinstance(self.__origin__[2],(float,int))):
            if len(self.__atoms__)>0:
                self.__get_origin__()
            else:
                self.__origin__=[0.0,0.0,0.0]
                            
        for i in range(3):
            disp[i]=val[i]-self.__origin__[i]
            
        self.__origin__=val

        if len(self.__atoms__)>0:
            for atom in self.__atoms__:
                    atom.displace(disp)
                    
            for sitelist in self.__ads_sites__.values():
                for site in sitelist:
                    site[0]+=disp[0]
                    site[1]+=disp[1]
            
            self.__get_z__()
        else:
            print("WARNING: No atom in the slab!")

    @property
    def latticevectors(self):
        return self.__lattice_vectors__
    
    @latticevectors.setter
    def latticevectors(self,val):
        if isinstance(val,list) and len(val)==2:
            if isinstance(val[0],list) and isinstance(val[1],list) and \
                len(val[0])==len(val[1])==2:
                self.__lattice_vectors__=val.copy()
            else:
                raise SlabError("Each component of the lattice vectors must be a list with two elements!")                
        else:
            raise SlabError("Lattice vectors must be provided as a list with two elements!")
            
    @property
    def ID(self):
        return self.__id__
    
class SlabError(BasicException):
    pass