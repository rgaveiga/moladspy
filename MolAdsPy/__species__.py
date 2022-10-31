from __future__ import print_function
from MolAdsPy.__exception__ import BasicException

SpecDict={"H":["hydrogen",1,1.0078],"He":["helium",2,4.0026],
          "Li":["lithium",3,6.941],"Be":["berylium",4,9.0122],
          "B":["boron",5,10.811],"C":["carbon",6,12.011],
          "N":["nitrogen",7,14.007],"O":["oxygen",8,15.999],
          "F":["fluorine",9,18.998],"Ne":["neon",10,20.18],
          "Na":["sodium",11,22.99],"Mg":["magnesium",12,24.305],
          "Al":["aluminium",13,26.982],"Si":["silicon",14,28.086],
          "P":["phosphorus",15,30.974],"S":["sulfur",16,32.065],
          "Cl":["chlorine",17,35.453],"Ar":["argon",18,39.948],
          "K":["potassium",19,39.098],"Ca":["calcium",20,40.078],
          "Sc":["scandium",21,44.956],"Ti":["titanium",22,47.867],
          "V":["vanadium",23,50.942],"Cr":["chromium",24,51.996],
          "Mn":["manganese",25,54.938],"Fe":["iron",26,55.845],
          "Co":["cobalt",27,58.933],"Ni":["nickel",28,58.693],
          "Cu":["copper",29,63.546],"Zn":["zinc",30,65.380],
          "Ga":["gallium",31,69.723],"Ge":["germanium",32,72.64],
          "As":["arsenic",33,74.922],"Se":["selenium",34,78.96],
          "Br":["bromine",35,79.904],"Kr":["krypton",36,83.709],
          "Rb":["rubidium",37,85.468],"Sr":["strontium",38,87.620],
          "Y":["yttrium",39,88.906],"Zr":["zirconium",40,91.224],
          "Nb":["niobium",41,92.906],"Mo":["molybdenum",42,95.95],
          "Tc":["technetium",43,98],"Ru":["ruthenium",44,101.07],
          "Rh":["rhodium",45,102.91],"Pd":["palladium",46,106.42],
          "Ag":["silver",47,107.87],"Cd":["cadmium",48,112.41],
          "In":["indium",49,114.82],"Sn":["tin",50,118.71],
          "Sb":["antimony",51,121.76],"Te":["tellurium",52,127.6],
          "I":["iodine",53,126.9],"Xe":["xenon",54,131.29],
          "Cs":["cesium",55,132.91],"Ba":["barium",56,137.33],
          "La":["lanthanum",57,138.91],"Ce":["cerium",58,140.12],
          "Pr":["praseodymium",59,140.91],"Nd":["neodymium",60,144.24],
          "Pm":["promethium",61,145],"Sm":["samarium",62,150.36],
          "Eu":["europium",63,151.96],"Gd":["gadolinium",64,157.25],
          "Tb":["terbium",65,158.93],"Dy":["dysprosium",66,162.5],
          "Ho":["holmium",67,164.93],"Er":["erbium",68,167.26],
          "Tm":["thulium",69,168.93],"Yb":["ytterbium",70,173.04],
          "Lu":["lutetium",71,174.97],"Hf":["hafnium",72,178.49],
          "Ta":["tantalum",73,180.95],"W":["tungsten",74,183.84],
          "Re":["rhenium",75,186.21],"Os":["osmium",76,190.23],
          "Ir":["iridium",77,192.22],"Pt":["platinum",195.08],
          "Au":["gold",79,196.97],"Hg":["mercury",80,200.59],
          "Tl":["thallium",81,204.38],"Pb":["lead",82,207.20],          
          "Bi":["bismuth",83,208.98]}

class Species:
    __species__=[]     # List of atomic species
    
    def __init__(self,symbol,element=None,atomicnumber=None,atomicmass=None,
                 pseudopotential=None):
        '''
        __init__(symbol,element,atomicnumber,atomicmass,pseudopotential) ->
            class constructor.

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

        Returns
        -------
        None.
        '''
        self.__pseudopotential__=None
        
        if isinstance(symbol,str) and len(symbol)<=4 and len(symbol)>0:
            self.__symbol__=symbol
            
            if symbol in SpecDict.keys():
                self.__element__=SpecDict[symbol][0]
                self.__atomicnumber__=SpecDict[symbol][1]
                self.__atomicmass__=SpecDict[symbol][2]
            else:
                if isinstance(element,str) and len(element)>0:
                    self.__element__=element
                else:
                    raise SpeciesError("'element' must be provided as a string!")
                    
                if isinstance(atomicnumber,int) and atomicnumber>0:
                    self.__atomicnumber__=atomicnumber
                else:
                    raise SpeciesError("'atomicnumber' must be provided as a positive integer!")                

                if isinstance(atomicmass,(int,float)) and atomicmass>0.0:
                    self.__atomicmass__=atomicmass
                else:
                    raise SpeciesError("'atomicmass' must be provided as a positive number!")                
        else:
            raise SpeciesError("'symbol' must be a provided as a string with up to 4 characters!")

        if pseudopotential is not None and isinstance(pseudopotential,str) and \
            len(pseudopotential)>0:
            self.__pseudopotential__=pseudopotential
        else:
            print("WARNING: 'pseudopotential' for element '%s', if any, must be provided as a string!"
                  % self.__element__)

        Species.__species__.append(self)

    '''
    Properties
    ----------
    symbol : string, readonly
        Symbol of the chemical element.
    element : string, readonly
        Name of the atomic element.
    atomicnumber : integer, readonly
        Atomic number of the element.
    atomicmass : float, readonly
        Atomic mass of the element, usually in atomic units.
    pseudopotential : string
        Name of the file containing the pseudopotential for the element,
        to be employed in DFT calculations.
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
    def pseudopotential(self):
        return self.__pseudopotential__
    
    @pseudopotential.setter
    def pseudopotential(self,val):
        if isinstance(val,str):
            self.__pseudopotential__=val
        else:
            raise SpeciesError("'pseudopotential' must be a string!")

class SpeciesError(BasicException):
    pass