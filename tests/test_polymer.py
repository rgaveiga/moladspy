import sys

# TODO: Make this computer dettached
path = "E:\\PBC\\GitHub\\moladspy\\MolAdsPy"

sys.path.append(path)
from core import Atom
from core import Polymer
from numpy import array
import unittest
import os

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
        filename = "pla1"
        filepath = os.path.join(path,filename + ".xyz")
        name_polymer = "pla"
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
    def test_rotate_axis(self):
        filename = "pla1"
        filepath = os.path.join(path,filename + ".xyz")
        pol=Polymer("pla",filepath,vaccuum=20.0)
        pol1 = pol.copy()
        #pol1.align('x')
        #pol1.write_xyz(filename+".align_1.1.xyz")
        pol1.align('y')
        print()
        print("Final orientation: " + str(pol1._orientation))
        pol1.write_xyz(filename+".align_1.2.xyz")
        print()
        print("Final orientation: " + str(pol1._orientation))
        pol1.align('z')
        print()
        print("Final orientation: " + str(pol1._orientation))

        pol1.write_xyz(filename+".align_1.3.xyz")
        print()
        print("Final orientation: " + str(pol1._orientation))
        pol1.align('x')
        pol1.write_xyz(filename+".align_1.4.xyz")
        print()
        print("Final orientation: " + str(pol1._orientation))


if __name__ == '__main__':
    unittest.main()