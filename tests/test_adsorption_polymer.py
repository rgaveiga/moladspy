import sys

# TODO: Make this computer dettached
path = "E:\\PBC\\GitHub\\moladspy\\MolAdsPy"

sys.path.append(path)
from core import Atom
from core import Polymer
from core import Slab
from core import Molecule
from core import Adsorption
from numpy import array
import numpy as np
from inspect import signature
import unittest
import os
import sys
import logging

# logging.basicConfig(level=logging.DEBUG, filename='logs.log', format = -)


'''

assertEqual(x, y, msg=None)	x == y
assertNotEqual(x,y,msg=None)	x != y
assertTrue(x, msg=None)	bool(x) is True
assertFalse(x, msg=None)	bool(x) is False
assertIs(x, y , msg=None)	x is y
assertIsNot(x, y, msg=None)	x is not y
assertIsNone(x, msg=None)	x is None
assertIsNotNone(x , msg=None)	x is not None
assertIn(x, y, msg=None)	x in y
assertNotIn(x, y, msg=None)	x not in y
assertIsInstance(x, y, msg=None)	isinstance(x, y)
assertNotIsInstance(x, y, msg=None)	not isinstance(x, y)

'''

class Singleton(object):
    _instance = None
    def __new__(class_, *args, **kwargs):
        if not isinstance(class_._instance, class_):
            class_._instance = object.__new__(class_, *args, **kwargs)
        return class_._instance

class GlobalsForTest(Singleton):
    def __init__(self):
        self.__loc__ = os.path.realpath(os.path.join(os.getcwd(),os.path.dirname(__file__)))
        self.results_folder = os.path.join(self.__loc__,"TestAdorptionPolymerResults")
        if not os.path.isdir(self.results_folder):
            os.makedirs(self.results_folder)
        self.filename_polymer = "pla1"
        self.filename_slab = "graphene"

        self.filepath_polymer = os.path.join(path,self.filename_polymer + ".xyz")
        self.filepath_slab = os.path.join(path,self.filename_slab + ".xyz")

        self.name_polymer = "pla"
        self.name_slab = "graphene"

        self.vacuum = 5.0

        self.graphene = Slab(self.name_slab,self.filename_slab + ".xyz")
        self.pol = Polymer(self.name_polymer,self.filepath_polymer,vaccuum=self.vacuum)


class TestPolymerSlabAdsorption(unittest.TestCase):
        # TODO: Correct align Y Resize
    def test_align_polymer_in_adsorption(self):
        '''
        This test will create a rectangular Slab with resize different for each size.

        - First it will be verified if the polymer is correctly aligned over and under the Slab.
        - Then, it will be calculated the minimum number of repetitions for the slab and polymer to fit the box.

        '''
        print("\n\n### Test Align Polymer Outputs ###\n\n")
        glob = GlobalsForTest()
        
        
        chosen_polymer = glob.pol       # Chosen Polymer for the test (Asserts will always be with already generated files)
        chosen_slab = glob.graphene     # Chosen Slab for the test
        initial_min_sep = 1.0           # min_sep that is initialized with Adsorption obj
        test_vertsep = 2.0              # vertsep that will be used during set_separation calls
        initial_vaccuum = 10.0          # Vaccumm that will be initialized in Adorption objs
        initial_slab_resize_x = 4       # X Resize for Slabs
        initial_slab_resize_y = 2       # Y Resize for Slabs
        
        # First Test : Check if Adsorption raises the correct error when trying to add a polymer incorrectly
        #   aligned with the Slab (z)
        first_test_axis_correct_expected = ['x','y']
        first_test_axis_error_expected = 'z'
        verbose = False
        
        axis = ['x','y','z']

        # Initialize Polymers with initial align in each axis
        pols = {ax:chosen_polymer.copy() for ax in axis}    # Create one polymer for each axis in a dictionary
        [pols[ax].align(ax) for ax in axis]                 # Align each polymer to its respective axis
        
        slabs = {ax:chosen_slab.copy() for ax in axis}      # Create individual copies of slab for each align case

        # Initialize a dictionary with one Adsorption obj for each axis
        adsorbeds =  {ax:Adsorption(slabs[ax],minsep=initial_min_sep,vaccuum=initial_vaccuum,verbose=verbose) for ax in axis}
        # Resize slab in adsorbed objs
        [adsorbeds[ax].resize(initial_slab_resize_x,initial_slab_resize_y) for ax in axis]
        
        
        # Test 1 - test is if polymers are accepted with wrong alignment

        for ax,ads in adsorbeds.items():
            try:
                ads.add_component(pols[ax],vertsep=test_vertsep,side="top")
                ads.write_xyz((os.path.join(glob.results_folder,"align_test.init_align_" + str(ax) + ".xyz")))
                print("Adsorption obj successfully accepted polymer initially aligned in " + ax)
            except Exception as e:
                print(e)
                print("Adsorption obj didn't accepted polymer initially aligned in " + ax)
                #self.assertTrue( any(ax in x for x in first_test_axis_error_expected),"Raised error correctly.")
                

        # TODO: Second test is if it is possible to realign wrong with polymer already in Adsorption object

        print()
        print()
        print("Entered Align in Ads")
        print()
        print()
        for _ax,pol in pols.items():
            for ax in axis:
                try:
                    pol.align(ax)
                    adsorbeds[_ax].write_xyz(os.path.join(glob.results_folder,"align_test." + str(_ax) + ".to_" + str(ax) + ".xyz"))
                    print("Polymer initial aligned with " + _ax + " was successfully aligned with " + ax + " axis.")
                except Exception as e:
                    adsorbeds[_ax].write_xyz(os.path.join(glob.results_folder,"align_test." + str(_ax) + ".to_" + str(ax) + ".xyz"))
                    print("Polymer initial aligned with " + _ax + " when tried to align with " + ax + " axis, it raised the error below:")
                    print(e)
                    #self.assertEqual(first_test_axis_error_expected,ax,"Raised error correctly.")

    def test_resize(self):
        glob = GlobalsForTest()
        tmp_pol = glob.pol.copy()
        tmp_pol_sized = glob.pol.copy()

        tmp_pol.write_xyz(os.path.join(glob.results_folder,"original.xyz"))
        tmp_slab = glob.graphene.copy()

        n = 5
        tmp_pol_sized.resize(n)
        for i in range(n):
            tmp_pol.resize(i)
            tmp_slab.resize(i,0)
            tmp_pol_sized.resize(n-i)
            tmp_pol.write_xyz(os.path.join(glob.results_folder,"resize_" + str(i) + ".xyz"))
            tmp_slab.write_xyz(os.path.join(glob.results_folder,"slab_resize_" + str(i) + ".xyz"))
            tmp_pol_sized.write_xyz(os.path.join(glob.results_folder,"resize_reverse_" + str(i) + ".xyz"))

    # TODO: Verify 3 Aligns as input
    # TODO: Test for Top and Bottom, several set_separation, without violation, and with violation.
    # TODO: Put verbose mode to say if violates or not
    def test_set_separation_polymer_in_adsorption(self):

        '''
        This test verifies if all loading parameters are correct for a polymer with a slab.
        It will be verified if the polymer pass over defined adsorption sites.
        '''
        print("\n\n### Test Load Polymer Outputs ###\n\n")
        glob = GlobalsForTest()
        
        pol1 = glob.pol.copy()
        graphene1 = glob.graphene.copy()

        separations = [0.5,1.0,2.0,3.0,4.0,5.0]

        adsorbed = Adsorption(graphene1,minsep=1.0,vaccuum=10.0)
        # Arbitrary initial Resize for Slab
        adsorbed.resize(4,2)
        adsorbed.add_component(pol1,vertsep=2.0,side="top")
        # Check in asssert if resize is always correct for certain mismatch and deformation



        adsorbed.set_separation(pol1,1.0)
        print(len(adsorbed._atoms),len(adsorbed._species))
        adsorbed.write_xyz("polymer_adsorbed.1.xyz")
    

    
    # TODO: Define flags of auto_correction of set_separation, align, etc...
    # TODO: Move "y" error
    def test_move_polymer_in_adsorption(self):
        '''
        This test will verify if the polymer can move freely, but respecting proximity constraints to the slab.
        '''
        print("\n\n### Test Move Polymer Outputs ###\n\n")
        glob = GlobalsForTest()
        pol_x = glob.pol.copy()
        pol_y = glob.pol.copy()
        pol_z = glob.pol.copy()
        pol_x.align('x')
        pol_y.align('y')
        pol_z.align('z')

        dic_pols = {
            'x': pol_x,
            'y': pol_y,
            'z': pol_z
        }

        move_x = array([1.0,0.0,0.0])
        move_y = array([0.0,1.0,0.0])
        move_z1 = array([0.0,0.0,1.0])
        move_z2 = array([0.0,0.0,-1.0])

        graphene_x = glob.graphene.copy()
        graphene_y = glob.graphene.copy()
        graphene_z = glob.graphene.copy()

        adsorbed_x = Adsorption(graphene_x,minsep=1.0,vaccuum=10.0)
        adsorbed_y = Adsorption(graphene_y,minsep=1.0,vaccuum=10.0)
        adsorbed_z = Adsorption(graphene_z,minsep=1.0,vaccuum=10.0)
        adsorbed_x.resize(4,2)
        adsorbed_y.resize(4,2)
        adsorbed_z.resize(4,2)

        dic_adsorbed = {
            'x': adsorbed_x,
            'y': adsorbed_y,
            'z': adsorbed_z
        }

        adsorbed_x.add_component(pol_x,vertsep=2.0,side="top")
        adsorbed_y.add_component(pol_y,vertsep=2.0,side="top")
        try:
            adsorbed_z.add_component(pol_z,vertsep=2.0,side="top")
        except:
            pass

        # Move x
        for key,_pol in dic_pols.items():
            _pol.move_to(move_x)
            dic_adsorbed[key].write_xyz('move_' + key + '_to_x.xyz')
        
        # Move y
        for key,_pol in dic_pols.items():
            _pol.move_to(move_y)
            dic_adsorbed[key].write_xyz('move_' + key + '_to_y.xyz')

        # Move z away
        for key,_pol in dic_pols.items():
            _pol.move_to(move_z1)
            dic_adsorbed[key].write_xyz('move_' + key + '_to_z_away.xyz')

        # Move z forward
        for key,_pol in dic_pols.items():
            _pol.move_to(move_z2)
            dic_adsorbed[key].write_xyz('move_' + key + '_to_x.xyz')

    def test_different_sep_polymer_in_adsorption(self):
        pass

    def test_rotate_polymer_in_adsorption(self):
        '''
        This test will rotate the polymer near the slab and verify if constraints are respected and the polymer readjusted to a correct height.
        '''
        print("\n\n### Test Rotate Polymer Outputs ###\n\n")
        glob = GlobalsForTest()
        filename = "pla_rotate_in_graphene"
        pol_x = glob.pol.copy()
        pol_y = glob.pol.copy()
        pol_z = glob.pol.copy()
        pol_x.align('x')
        pol_y.align('y')
        pol_z.align('z')

        graphene_x = glob.graphene.copy()
        graphene_y = glob.graphene.copy()
        graphene_z = glob.graphene.copy()

        adsorbed_x = Adsorption(graphene_x,minsep=1.0,vaccuum=10.0)
        adsorbed_y = Adsorption(graphene_y,minsep=1.0,vaccuum=10.0)
        adsorbed_z = Adsorption(graphene_z,minsep=1.0,vaccuum=10.0)
        adsorbed_x.resize(4,2)
        adsorbed_y.resize(4,2)
        adsorbed_z.resize(4,2)

        adsorbed_x.add_component(pol_x,vertsep=2.0,side="top")
        adsorbed_y.add_component(pol_y,vertsep=2.0,side="top")
        try:
            adsorbed_z.add_component(pol_z,vertsep=2.0,side="top")
        except:
            pass

        degrees = np.arange(60,361,60)
        print('Angles are: ' + str(degrees))
        # Rotate in x
        for degree in degrees:
            pol_x.rotate(degree)
            adsorbed_x.write_xyz(filename+".align_x.degree_" + str(degree) + ".xyz")
            pol_y.rotate(degree)
            adsorbed_y.write_xyz(filename+".align_y.degree_" + str(degree) + ".xyz")
            pol_z.rotate(degree)
            adsorbed_z.write_xyz(filename+".align_z.degree_" + str(degree) + ".xyz")

# TODO: Quickly implement to test
class TestPolymerMoleculeSlabAdsorption(unittest.TestCase):
    pass

# TODO: Trabalho de mestrado: Polipropileno em duas direções, duas rotaçoes com metil pra baixo, 5 separações.
class CasePolypropileneGraphene(unittest.TestCase):
    @unittest.skip("Just when everything is working")
    def GenerateSystems(self):
        filename_polymer = "pp"
        filename_slab = "graphene"

        filepath_polymer = os.path.join(path,filename_polymer + ".xyz")
        filepath_slab = os.path.join(path,filename_slab + ".xyz")

        name_polymer = "polypropilene"
        name_slab = "graphene"

        vacuum = 10.0

        slab = Slab(name_slab,filename_slab + ".xyz")
        pol = Polymer(name_polymer,filepath_polymer,vaccuum=vacuum)
        adsorbed = Adsorption(slab,minsep=1.0,vaccuum=10.0)
        adsorbed.resize(5,5)

        adsorbed.add_component(pol,vertsep=1.0,side="top")
        # Check  in asssert if resize is always correct for certain mismatch and deformation

        # Parameters for combinations
        angles = [0.0,180.0]
        separations = [1.0,2.0,3.0,4.0,5.0]
        aligns = ['x','y']

        for align in aligns:
            pol.align(align)
            for sep in separations:
                adsorbed.set_separation(pol,sep)
                for angle in angles:
                    pol.rotate(angle)
                    adsorbed.write_xyz("ads_pp_graphene.align_" + align + ".sep_" + str(sep) + ".angle_" + str(angle) + ".xyz")


if __name__=='__main__':
     #unittest.main()
    tmp = TestPolymerSlabAdsorption()
    tmp.test_align_polymer_in_adsorption()
    #tmp.test_resize()




