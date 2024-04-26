from MolAdsPy import Polymer
import numpy as np
import unittest

# TODO: Make a virtual simple polymer and test align, rotate and resize
# TODO: Test same tests on PLA
# TODO: Test same tests on Polypropilene
# TODO: Make Assert Files on results and Unit Tests

class TestPolymerClassMethods(unittest.TestCase):
    def test_readwrite(self):
        pass

    def test_rotate(self):
        pass

    def test_align(self):
        '''
        print("1)")

        pol=Polymer("pla",filepath,vaccuum=20.0)


        pol.write_xyz(filename + ".1.xyz")

        print("\n9)")

        pol3=pol.copy()
        pol2=pol.copy()

        pol2._rotate(-30,-30,-30)
        pol3._rotate(30,30,30)
        pol2.write_xyz(filename + ".2.xyz")
        pol3.write_xyz(filename + ".3.xyz")

        print("Max ",pol.maxx,pol.maxy,pol.maxz)
        print("Min ",pol.minx,pol.miny,pol.minz)
        print("COM ",pol.center_of_mass)

        pol_xx = pol.copy()
        pol_xx.align('x')
        pol_xx.write_xyz(filename + ".align_xx.xyz")

        pol_xy = pol.copy()
        pol_xy.align('y')
        pol_xy.write_xyz(filename + ".align_xy.xyz")

        pol_xz = pol.copy()
        pol_xz.align('z')
        pol_xz.write_xyz(filename + ".align_xz.xyz")
        '''
    def test_align_axis(self):
        pol=Polymer("pla","pla1.xyz",vaccuum=20.0)
        pol1 = pol.copy()
        
        pol1.align('x')
        self.assertEqual(pol1._orientation,0,msg="Wrong orientation when aligning to x.")
        pol1.write_xyz("pol.align_1.1.xyz")

        pol1.align('y')
        self.assertEqual(pol1._orientation,1,msg="Wrong orientation when aligning to y.")
        pol1.write_xyz("pol.align_1.2.xyz")

        pol1.align('z')
        self.assertEqual(pol1._orientation,2,msg="Wrong orientation when aligning to z.")
        pol1.write_xyz("pol.align_1.3.xyz")
        
        pol1.align('x')
        self.assertEqual(pol1._orientation,0,msg="Wrong orientation when aligning to x.")
        pol1.write_xyz("pol.align_1.4.xyz")

    def test_rotate_axis(self):

        # This test will rotate 60 degrees on each align 5 times
        pol=Polymer("pla","pla1.xyz",vaccuum=20.0)
        pol1 = pol.copy()
        
        degrees = np.arange(60,361,60)
        print(degrees)

        pol1.align('x')
        for degree in degrees:
            pol1.rotate(degree)
            pol1.write_xyz("pol.align_1.1.degree_" + str(degree) + ".xyz")
        
        pol1.align('y')
        for degree in degrees:
            pol1.rotate(degree)
            pol1.write_xyz("pol.align_1.2.degree_" + str(degree) + ".xyz")
        
        pol1.align('z')
        for degree in degrees:
            pol1.rotate(degree)
            pol1.write_xyz("pol.align_1.3.degree_" + str(degree) + ".xyz")





if __name__ == '__main__':
    unittest.main()