from MolAdsPy.__species__ import Species
from MolAdsPy.__exception__ import BasicException
from copy import deepcopy        

class Atom:
    __curratomid__=0     # ID available for a new atom
    
    def __init__(self,species,x=[0.0,0.0,0.0],fixed=[False,False,False]):
        '''
        __init__(species,x,fixed) -> class constructor.

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

        Returns
        -------
        None.
        '''
        if isinstance(species,Species):
            self.__symbol__=species.__symbol__
            self.__element__=species.__element__
            self.__atomicnumber__=species.__atomicnumber__
            self.__atomicmass__=species.__atomicmass__
        else:
            raise AtomError("'species' must be a Species object!")
            
        if isinstance(x,list) and len(x)==3 and \
            isinstance(x[0],(int,float)) and \
            isinstance(x[1],(int,float)) and \
            isinstance(x[2],(int,float)):
            self.__x__=x.copy()
        else:
            raise AtomError("'x' must be a list with three components!")
                        
        if isinstance(fixed,list) and len(fixed)==3 and \
            isinstance(fixed[0],bool) and \
            isinstance(fixed[1],bool) and \
            isinstance(fixed[2],bool):
            self.__fixed__=fixed.copy()
        else:
            raise AtomError("'fixed' must be a list with three components!")
                        
        self.__id__=Atom.__curratomid__
        Atom.__curratomid__+=1        
            
    def displace(self,disp):
        '''
        displace(disp) -> displaces the atom from its current position.

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
            for i in range(3):
                self.__x__[i]+=disp[i]
        else:
            raise AtomError("'disp' must be a list with three components!")
            
    def move_to(self,x):
        '''
        move_to(x) -> moves the atom to a new position.

        Parameters
        ----------
        x : Python list
            New atom position.

        Returns
        -------
        None.
        '''       
        if isinstance(x,list) and len(x)==3 and \
            isinstance(x[0],(float,int)) and \
            isinstance(x[1],(float,int)) and \
            isinstance(x[2],(float,int)):
            for i in range(3):
                self.__x__[i]=x[i]
        else:
            raise AtomError("'x' must be a list with three components!")
            
    def copy(self):
        '''
        copy() -> returns a copy of the current Atom object.

        Returns
        -------
        Atom object
            Copy of the atom.
        '''
        newatom=deepcopy(self)
        newatom.__id__=Atom.__curratomid__
        Atom.__curratomid__+=1
        
        return newatom

    '''
    Properties
    ----------
    symbol : string, readonly
        Symbol of the chemical element.
    element : string, readonly
        Name of the atomic element.
    atomicnumber : int, readonly
        Atomic number of the element.
    atomicmass : float, readonly
        Atomic mass of the element, usually in atomic units.
    coords : Python list
        Cartesian coordinates of the atom.
    isfixed : Python list
        Defines if an atom component should remain unaltered in the course
        of a geometry optimization.
    ID : int, readonly
        Unique atom identifier.
    '''    
    @property
    def symbol(self):
        return self.__symbol__
    
    @property
    def element(self):
        return self.__element__
    
    @property
    def atomicnumber(self):
        return self.__atomicnumber__
    
    @property
    def atomicmass(self):
        return self.__atomicmass__
    
    @property
    def coords(self):
        return self.__x__
    
    @coords.setter
    def coords(self,val):
        if isinstance(val,list) and len(val)==3 and \
            isinstance(val[0],(int,float)) and \
            isinstance(val[1],(int,float)) and \
            isinstance(val[2],(int,float)):
            self.__x__=val.copy()
        else:
            raise AtomError("'x' must be a list with three components!") 
            
    @property
    def isfixed(self):
        return self.__fixed__
    
    @isfixed.setter
    def isfixed(self,val):
        if isinstance(val,list) and len(val)==3 and \
            isinstance(val[0],bool) and \
            isinstance(val[1],bool) and \
            isinstance(val[2],bool):
            self.__fixed__=val
        else:
            raise AtomError("'fixed' must be a list with three components!")
            
    @property
    def ID(self):
        return self.__id__
    
class AtomError(BasicException):
    pass