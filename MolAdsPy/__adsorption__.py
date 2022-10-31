from __future__ import print_function
from MolAdsPy.__species__ import Species
from MolAdsPy.__molecule__ import Molecule
from MolAdsPy.__slab__ import Slab
from MolAdsPy.__exception__ import BasicException
from math import inf,sqrt
from multipledispatch import dispatch
import os.path

class Adsorption:
    def __init__(self,substrate):
        '''
        __init__(substrate) -> Class constructor.
        
        Parameters
        ----------
        substrate : Slab object
            Substrate where one or more molecules can be adsorbed.

        Returns
        -------
        None.
        '''
        self.__species__=[]
        
        if isinstance(substrate,Slab):
            self.__substrate__=substrate
            
            for symbol in substrate.__species__:
                if not symbol in self.__species__:
                    self.__species__.append(symbol)
        else:
            raise AdsorptionError("Substrate must be a Slab object!")
        
        self.__molecules__=[]
        self.__minsep__=1.0
    
    @dispatch(Molecule,str,int,str,(int,float))    
    def add_molecule(self,molecule,ads_site,ads_site_index,anchor,vert_sep):
        '''
        add_molecule(molecule,ads_site,ads_site_index,anchor,vert_sep) -> 
            adsorbs a molecule onto the substrate.

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

        Returns
        -------
        None.
        '''
        xnew=[0.0,0.0,0.0]
        
        if self.__substrate__ is None:
            raise AdsorptionError("No substrate has been defined!")
        elif not isinstance(molecule,Molecule):
            raise AdsorptionError("'molecule' must be a Molecule object!")
        elif not isinstance(anchor,str) or len(anchor)==0 or \
            anchor not in molecule.__anchors__.keys():
            raise AdsorptionError("'anchor' must be a valid anchor point!")
        elif not isinstance(ads_site,str):
            raise AdsorptionError("'ads_site' must be a string!")        
        elif not isinstance(ads_site_index,int):
            raise AdsorptionError("'ads_site_index' must be an integer!")              
        elif (not isinstance(vert_sep,(int,float))) or vert_sep<=0.0:
            raise AdsorptionError("'vert_sep' must be a number greater than zero!")
        
        if len(ads_site)>0 and \
            (ads_site in self.__substrate__.__ads_sites__.keys()) and \
            ads_site_index>=0 and \
            ads_site_index<len(self.__substrate__.__ads_sites__[ads_site]):
            xnew[0]=self.__substrate__.__ads_sites__[ads_site][ads_site_index][0]
            xnew[1]=self.__substrate__.__ads_sites__[ads_site][ads_site_index][1]
        else:
            raise AdsorptionError("A valid adsorption site must be provided!")
                   
        self.__molecules__.append(molecule)
        
        for symbol in molecule.__species__:
            if not symbol in self.__species__:
                self.__species__.append(symbol)
        
        zmin=self.__substrate__.top+vert_sep
        dz=zmin-molecule.minz
        xnew[2]=molecule.__anchors__[anchor][2]+dz
                
        molecule.move_to(xnew,anchor)
        
    @dispatch(Molecule,list,str,(int,float))
    def add_molecule(self,molecule,pos,anchor,vert_sep):
        '''
        add_molecule(molecule,pos,anchor,vert_sep) -> adsorbs a molecule onto 
            the substrate.

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

        Returns
        -------
        None.
        '''
        xnew=[0.0,0.0,0.0]
        
        if self.__substrate__ is None:
            raise AdsorptionError("No substrate has been defined!")
        elif not isinstance(molecule,Molecule):
            raise AdsorptionError("'molecule' must be a Molecule object!")
        elif not isinstance(anchor,str) or len(anchor)==0 or \
            anchor not in molecule.__anchors__.keys():
            raise AdsorptionError("'anchor' must be a valid anchor point!")
        elif (not isinstance(vert_sep,(int,float))) or vert_sep<=0.0:
            raise AdsorptionError("'vert_sep' must be a number greater than zero!")
        
        if isinstance(pos,list) and len(pos)==2 and \
            isinstance(pos[0],(int,float)) and isinstance(pos[1],(int,float)):
            xnew[0]=pos[0]
            xnew[1]=pos[1]
        else:
            raise AdsorptionError("'pos' must be provided as a list with two components!")
                   
        self.__molecules__.append(molecule)
        
        for symbol in molecule.__species__:
            if not symbol in self.__species__:
                self.__species__.append(symbol)
        
        zmin=self.__substrate__.top+vert_sep
        dz=zmin-molecule.minz
        xnew[2]=molecule.__anchors__[anchor][2]+dz
                
        molecule.move_to(xnew,anchor)
    
    @dispatch(Molecule)    
    def remove_molecule(self,molecule):
        '''
        remove_molecule(molecule) -> removes an adsorbed molecule.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule to be removed.

        Returns
        -------
        None.
        '''
        if isinstance(molecule,Molecule):
            if molecule in self.__molecules__:
                self.__molecules__.remove(molecule)
                
                for symbol in molecule.__species__:
                    remove_species=True
                    
                    if symbol in self.__substrate__.__species__:
                        remove_species=False
                    else:
                        for other_molecule in self.__molecules__:
                            if symbol in other_molecule.__species__:
                                remove_species=False
                                
                                break
                    
                    if remove_species:
                        self.__species__.remove(symbol)
            else:
                raise AdsorptionError("Molecule not found!")
        else:
            raise AdsorptionError("'molecule' must be a Molecule object!")
            
    @dispatch(int)
    def remove_molecule(self,molid):
        '''
        remove_molecule(molid) -> removes an adsorbed molecule.    

        Parameters
        ----------
        molid : int
            ID of the adsorbed molecule.

        Returns
        -------
        None.
        '''
        if isinstance(molid,int) and molid>=0 and molid<Molecule.__currmolid__:
            for molecule in self.__molecules__:
                if molecule.__id__==molid:
                    self.__molecules__.remove(molecule)
                    
                    for symbol in molecule.__species__:
                        remove_species=True
                        
                        if symbol in self.__substrate__.__species__:
                            remove_species=False
                        else:
                            for other_molecule in self.__molecules__:
                                if symbol in other_molecule.__species__:
                                    remove_species=False
                                    
                                    break
                        
                        if remove_species:
                            self.__species__.remove(symbol)
                                                
                    break
            else:
                raise AdsorptionError("Molecule ID not found!")
        else:
            raise AdsorptionError("'molid' must be an integer greater than or equal to zero!")

    @dispatch()
    def enforce_minsep(self):
        '''
        enforce_minsep() -> if necessary, displaces vertically a molecule such 
            that the minimum separation of the molecule with respect to the 
            top of the substrate is the prescribed minimum distance.

        Returns
        -------
        None.
        '''
        for molecule in self.__molecules__:
            if molecule.minz-self.__substrate__.top<self.__minsep__:
                dispz=self.__minsep__-(molecule.minz-self.__substrate__.top)
                
                molecule.displace([0.0,0.0,dispz])

    @dispatch(int)
    def enforce_minsep(self,molid):
        '''
        enforce_minsep(molid) -> if necessary, displaces vertically a molecule 
            such that the minimum separation of the molecule with respect to 
            the top of the substrate is the prescribed minimum distance.

        Parameters
        ----------
        molid : int
            ID of the molecule.

        Returns
        -------
        None.
        '''       
        if isinstance(molid,int) and molid>=0 and molid<Molecule.__currmolid__:
            for molecule in self.__molecules__:
                if molecule.__id__==molid:
                    if molecule.minz-self.__substrate__.top<self.__minsep__:
                        dispz=self.__minsep__-(molecule.minz-
                                                self.__substrate__.top)
                
                        molecule.displace([0.0,0.0,dispz])
                        
                        break
            else:
                raise AdsorptionError("Molecule ID not found!")
        else:
            raise AdsorptionError("'molid' must be an integer greater than or equal to zero!")
    
    @dispatch(Molecule)                
    def enforce_minsep(self,molecule):
        '''
        enforce_minsep(molecule) -> if necessary, displaces vertically a 
            molecule such that the minimum separation of the molecule with 
            respect to the top of the substrate is the prescribed minimum 
            distance.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.

        Returns
        -------
        None.
        '''
        if isinstance(molecule,Molecule):
            if molecule in self.__molecules__:
                if molecule.minz-self.__substrate__.top<self.__minsep__:
                    dispz=self.__minsep__-(molecule.minz-
                                            self.__substrate__.top)
                
                    molecule.displace([0.0,0.0,dispz])
            else:
                raise AdsorptionError("Molecule not found!")
        else:
            raise AdsorptionError("'molecule' must be a Molecule object!")
                
    def write_xyz(self,filename="adsorbed.xyz",vacuum=10.0):
        '''
        write_xyz(filename,vacuum) -> saves the atomic coordinates of the 
            adsorbed system into an XYZ file.

        Parameters
        ----------
        filename : string, optional
            Name of the XYZ file. The default is "adsorbed.xyz".
        vacuum : float, optional
            Height of the vacuum layer to be added to the top of the adsorbed 
            system. The default is 10.0 (Angstroms).

        Returns
        -------
        None.
        '''
        maxz=-inf
        natoms=len(self.__substrate__.__atoms__)
        
        for molecule in self.__molecules__:
            natoms+=len(molecule.__atoms__)
            
            if molecule.maxz>maxz:
                maxz=molecule.maxz
                
        with open(filename,"w") as f:
            f.write("%d\n" % natoms)
            f.write('Lattice="%.6f %.6f 0.0 %.6f %.6f 0.0 0.0 0.0 %.6f" Properties=species:S:1:pos:R:3\n' 
                    % (self.__substrate__.__lattice_vectors__[0][0],
                       self.__substrate__.__lattice_vectors__[0][1],
                       self.__substrate__.__lattice_vectors__[1][0],
                       self.__substrate__.__lattice_vectors__[1][1],
                       maxz+vacuum))
            
            for atom in self.__substrate__.__atoms__:
                f.write("%s %.6f %.6f %.6f\n" % (atom.__symbol__,atom.__x__[0],
                                                 atom.__x__[1],atom.__x__[2]))
                
            for molecule in self.__molecules__:
                for atom in molecule.__atoms__:
                    f.write("%s %.6f %.6f %.6f\n" % (atom.__symbol__,
                                                     atom.__x__[0],
                                                     atom.__x__[1],
                                                     atom.__x__[2]))
                    
    def write_pw_input(self,filename="adsorbed.in",vacuum=20.0,**kwargs):
        '''
        write_pw_input(filename,vacuum,**kwargs) -> creates a basic input file
            for geometry relaxation of the adsorbed system using the pw.x code 
            in the Quantum Espresso package.

        Parameters
        ----------
        filename : string, optional
            Name of the input file. The default is "adsorbed.in".
        vacuum : float, optional
            Vacuum size separating the top of the highest molecule from the
            bottom of the image of the slab. The default is 20.0 angstroms.
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
                    types in the system. The default is 1.0 Bohr magneton 
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
        celldm=sqrt(pow((self.__substrate__.__lattice_vectors__[0][0]),2)+
                    pow((self.__substrate__.__lattice_vectors__[0][1]),2))
        inpstr+="    celldm(1)=%.4f,\n" % (celldm*1.8897259886)
        maxz=-inf
        natoms=len(self.__substrate__.__atoms__)
        
        for molecule in self.__molecules__:
            natoms+=len(molecule.__atoms__)
            
            if molecule.maxz>maxz:
                maxz=molecule.maxz
        
        inpstr+="    nat=%d,\n" % natoms
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
        inpstr+="%.4f   %.4f   0.0000\n" % (self.__substrate__.__lattice_vectors__[0][0]
                                            /celldm,
                                            self.__substrate__.__lattice_vectors__[0][1]
                                            /celldm)
        inpstr+="%.4f   %.4f   0.0000\n" % (self.__substrate__.__lattice_vectors__[1][0]
                                            /celldm,
                                            self.__substrate__.__lattice_vectors__[1][1]
                                            /celldm)
        inpstr+="0.0000   0.0000   %.4f\n" % ((maxz+vacuum)/celldm)

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
            for species in Species.__species__:
                if species.__symbol__==symbol:
                    if species.__pseudopotential__=="":
                        pseudopotential=symbol+".UPF"
                    else:
                        pseudopotential=species.__pseudopotential__
                    
                    inpstr+="%s %.4f %s\n" % (symbol,species.__atomicmass__,
                                              pseudopotential)
                    
                    break
                
        inpstr+="ATOMIC_POSITIONS angstrom\n"
        
        for atom in self.__substrate__.__atoms__:                    
            inpstr+="%s %.6f %.6f %.6f %d %d %d\n" % (atom.__symbol__,
                                                      atom.__x__[0],
                                                      atom.__x__[1],
                                                      atom.__x__[2],
                                                      int(not(atom.__fixed__[0])),
                                                      int(not(atom.__fixed__[1])),
                                                      int(not(atom.__fixed__[2])))
            
        for molecule in self.__molecules__:
            for atom in molecule.__atoms__:
                inpstr+="%s %.6f %.6f %.6f %d %d %d\n" % (atom.__symbol__,
                                                          atom.__x__[0],
                                                          atom.__x__[1],
                                                          atom.__x__[2],
                                                          int(not(atom.__fixed__[0])),
                                                          int(not(atom.__fixed__[1])),
                                                          int(not(atom.__fixed__[2])))
            
        with open(filename,"w") as f:
            f.write(inpstr)
 
    @dispatch(Molecule,list)
    def displace_molecule(self,molecule,disp):
        '''
        displace_molecule(molecule,disp) -> rigidly displaces the molecule.

        Parameters
        ----------
        molecule : Molecule object
            Molecule that will be displaced to a new position.
        disp : Python list
            Displacement vector.

        Returns
        -------
        None.
        '''
        if isinstance(molecule,Molecule):
            if molecule in self.__molecules__:
                molecule.displace(disp)
                self.enforce_minsep(molecule)
            else:
                raise AdsorptionError("Molecule not found!")
        else:
            raise AdsorptionError("'molecule' must be a Molecule object!")                
    
    @dispatch(int,list)        
    def displace_molecule(self,molid,disp):
        '''
        displace_molecule(molid,disp) -> rigidly displaces the molecule.

        Parameters
        ----------
        molid : int
            ID of the molecule that will be displaced to a new position.
        disp : Python list
            Displacement vector.

        Returns
        -------
        None.
        '''
        if isinstance(molid,int) and molid>=0 and molid<Molecule.__currmolid__:
            for molecule in self.__molecules__:
                if molecule.__id__==molid:
                    molecule.displace(disp)
                    self.enforce_minsep(molecule)
                    
                    break
            else:
                raise AdsorptionError("Molecule ID not found!")
        else:
            raise AdsorptionError("'molid' must be an integer greater than or equal to zero!")

    @dispatch(int,list,str)            
    def move_molecule_to(self,molid,x,anchor):
        '''
        move_molecule_to(molid,x,anchor) -> rigidly moves the molecule such 
            that the anchor point is located at 'x'.

        Parameters
        ----------
        molid : int
            ID of the adsorbed molecule.
        x : Python list
            New position of the molecule's anchor point.
        anchor : string
            Anchor point used as reference for translations.

        Returns
        -------
        None.
        '''
        if isinstance(molid,int) and molid>=0 and molid<Molecule.__currmolid__:
            for molecule in self.__molecules__:
                if molecule.__id__==molid:
                    if not isinstance(anchor,str) or len(anchor)==0 or \
                        anchor not in molecule.__anchors__.keys():
                        raise AdsorptionError("'anchor' must be a valid anchor point!")
                    
                    molecule.move_to(x,anchor)
                    self.enforce_minsep(molecule)
                    
                    break
            else:
                raise AdsorptionError("Molecule ID not found!")
        else:
            raise AdsorptionError("'molid' must be an integer greater than or equal to zero!")
            
    @dispatch(Molecule,list,str)            
    def move_molecule_to(self,molecule,x,anchor):
        '''
        move_molecule_to(molecule,x,anchor) -> rigidly moves the molecule such 
            that the anchor point is located at 'x'.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        x : Python list
            New position of the molecule's anchor point.
        anchor : string
            Anchor point used as reference for translations.

        Returns
        -------
        None.
        '''
        if isinstance(molecule,Molecule):
            if molecule in self.__molecules__:
                if not isinstance(anchor,str) or len(anchor)==0 or \
                    anchor not in molecule.__anchors__.keys():
                    raise AdsorptionError("'anchor' must be a valid anchor point!")
                
                molecule.move_to(x,anchor)
                self.enforce_minsep(molecule)                    
            else:
                raise AdsorptionError("Molecule not found!")
        else:
            raise AdsorptionError("'molecule' must be a Molecule object!")
            
    @dispatch(int,str,int,str,(int,float))
    def move_molecule_to(self,molid,ads_site,ads_site_index,anchor,vert_sep):
        '''
        move_molecule_to(molid,x,anchor) -> rigidly moves the molecule such 
            that the anchor point is located at the specified adsorption site.

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

        Returns
        -------
        None.
        '''
        xnew=[0.0,0.0,0.0]
        
        if not isinstance(ads_site,str):
            raise AdsorptionError("'ads_site' must be a string!")        
        elif not isinstance(ads_site_index,int):
            raise AdsorptionError("'ads_site_index' must be an integer!")              
        elif (not isinstance(vert_sep,(int,float))) or vert_sep<=0.0:
            raise AdsorptionError("'vert_sep' must be a number greater than zero!")
        
        if len(ads_site)>0 and \
            (ads_site in self.__substrate__.__ads_sites__.keys()) and \
            ads_site_index>=0 and \
            ads_site_index<len(self.__substrate__.__ads_sites__[ads_site]):
            xnew[0]=self.__substrate__.__ads_sites__[ads_site][ads_site_index][0]
            xnew[1]=self.__substrate__.__ads_sites__[ads_site][ads_site_index][1]
        else:
            raise AdsorptionError("A valid adsorption site must be provided!")
                           
        if isinstance(molid,int) and molid>=0 and molid<Molecule.__currmolid__:
            for molecule in self.__molecules__:
                if molecule.__id__==molid:
                    if not isinstance(anchor,str) or len(anchor)==0 or \
                        anchor not in molecule.__anchors__.keys():
                        raise AdsorptionError("'anchor' must be a valid anchor point!")
                    
                    zmin=self.__substrate__.top+vert_sep
                    dz=zmin-molecule.minz
                    xnew[2]=molecule.__anchors__[anchor][2]+dz
                
                    molecule.move_to(xnew,anchor)
                    self.enforce_minsep(molecule)
                    
                    break
            else:
                raise AdsorptionError("Molecule ID not found!")
        else:
            raise AdsorptionError("'molid' must be an integer greater than or equal to zero!")
            
    @dispatch(Molecule,str,int,str,(int,float))
    def move_molecule_to(self,molecule,ads_site,ads_site_index,anchor,vert_sep):
        '''
        move_molecule_to(molecule,x,anchor) -> rigidly moves the molecule such 
            that the anchor point is located at the specified adsorption site.

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

        Returns
        -------
        None.
        '''
        xnew=[0.0,0.0,0.0]
        
        if not isinstance(ads_site,str):
            raise AdsorptionError("'ads_site' must be a string!")        
        elif not isinstance(ads_site_index,int):
            raise AdsorptionError("'ads_site_index' must be an integer!")              
        elif (not isinstance(vert_sep,(int,float))) or vert_sep<=0.0:
            raise AdsorptionError("'vert_sep' must be a number greater than zero!")
        
        if len(ads_site)>0 and \
            (ads_site in self.__substrate__.__ads_sites__.keys()) and \
            ads_site_index>=0 and \
            ads_site_index<len(self.__substrate__.__ads_sites__[ads_site]):
            xnew[0]=self.__substrate__.__ads_sites__[ads_site][ads_site_index][0]
            xnew[1]=self.__substrate__.__ads_sites__[ads_site][ads_site_index][1]
        else:
            raise AdsorptionError("A valid adsorption site must be provided!")
                           
        if isinstance(molecule,Molecule):
            if molecule in self.__molecules__:
                if not isinstance(anchor,str) or len(anchor)==0 or \
                    anchor not in molecule.__anchors__.keys():
                    raise AdsorptionError("'anchor' must be a valid anchor point!")
                    
                zmin=self.__substrate__.top+vert_sep
                dz=zmin-molecule.minz
                xnew[2]=molecule.__anchors__[anchor][2]+dz
            
                molecule.move_to(xnew,anchor)
                self.enforce_minsep(molecule)                    
            else:
                raise AdsorptionError("Molecule not found!")
        else:
            raise AdsorptionError("'molecule' must be a Molecule object!")
            
    @dispatch(Molecule,(int,float),(int,float),(int,float),str)        
    def rotate_molecule(self,molecule,theta,phi,psi,anchor):
        '''
        rotate_molecule(molecule,theta,phi,psi,anchor) -> rotates the molecule 
            around an anchor point.

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

        Returns
        -------
        None.
        '''        
        if isinstance(molecule,Molecule):
            if molecule in self.__molecules__:
                molecule.rotate(theta,phi,psi,anchor)
                self.enforce_minsep(molecule)                    
            else:
                raise AdsorptionError("Molecule not found!")
        else:
            raise AdsorptionError("'molecule' must be a Molecule object!")
    
    @dispatch(int,(int,float),(int,float),(int,float),str)        
    def rotate_molecule(self,molid,theta,phi,psi,anchor):
        '''
        rotate_molecule(molid,theta,phi,psi,anchor) -> rotates the molecule 
            around an anchor point.

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

        Returns
        -------
        None.
        '''        
        if isinstance(molid,int) and molid>=0 and molid<Molecule.__currmolid__:
            for molecule in self.__molecules__:
                if molecule.__id__==molid:
                    molecule.rotate(theta,phi,psi,anchor)
                    self.enforce_minsep(molid)
                    
                    break
            else:
                raise AdsorptionError("Molecule ID not found!")
        else:
            raise AdsorptionError("'molid' must be an integer greater than or equal to zero!")

    @dispatch(Molecule,(int,float))
    def set_separation(self,molecule,sep):
        '''
        set_separation(molecule,sep) -> rigidly displaces the molecule along 
            the direction perpendicular to the substrate, such that the 
            sepatation between the top of the substrate and the nearest atom 
            in the molecule be equal to 'sep'.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.

        Returns
        -------
        None.
        '''
        if isinstance(molecule,Molecule):
            if molecule in self.__molecules__:
                currsep=molecule.minz-self.__substrate__.top
                dz=sep-currsep
                molecule.displace([0.0,0.0,dz])
                self.enforce_minsep(molecule)                    
            else:
                raise AdsorptionError("Molecule not found!")
        else:
            raise AdsorptionError("'molecule' must be a Molecule object!")
            
    @dispatch(int,(int,float))
    def set_separation(self,molid,sep):
        '''
        set_separation(molid,sep) -> rigidly displaces the molecule along 
            the direction perpendicular to the substrate, such that the 
            sepatation between the top of the substrate and the nearest atom 
            in the molecule be equal to 'sep'.

        Parameters
        ----------
        molid : int
            ID of the adsorbed molecule.
        sep : float
            Nearest distance separating the molecule from the substrate.

        Returns
        -------
        None.
        '''
        if isinstance(molid,int) and molid>=0 and molid<Molecule.__currmolid__:
            for molecule in self.__molecules__:
                if molecule.__id__==molid:
                    currsep=molecule.minz-self.__substrate__.top
                    dz=sep-currsep
                    molecule.displace([0.0,0.0,dz])
                    self.enforce_minsep(molecule)                    
            else:
                raise AdsorptionError("Molecule ID not found!")
        else:
            raise AdsorptionError("'molid' must be an integer greater than or equal to zero!")
            
    def list_molecules(self):
        '''
        list_molecules() -> prints a list of the adsorbed molecules.

        Returns
        -------
        None.
        '''
        print("MoleculeID, MoleculeType, XCM, YCM, ZCM, Separation\n------")
        
        for molecule in self.__molecules__:
            sep=molecule.minz-self.__substrate__.top
            print("%d, %s, %f, %f, %f, %f" % (molecule.ID,molecule.moleculetype,
                                              molecule.centerofmass[0],
                                              molecule.centerofmass[1],
                                              molecule.centerofmass[2],
                                              sep))

    '''
    Properties
    ----------
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
    '''            
    @property
    def substrate(self):
        return self.__substrate__
    
    @property
    def adsorbedmolecules(self):
        return self.__molecules__
                
    @property
    def atomtypes(self):
        return self.__species__
        
    @property
    def minimumseparation(self):
        return self.__minsep__
    
    @minimumseparation.setter
    def minimumseparation(self,val):
        if isinstance(val,(float,int)) and val>0.0:
            self.__minsep__=val
        else:
            raise AdsorptionError("Minimum separation must be a number greater than zero!")
                           
class AdsorptionError(BasicException):
    pass