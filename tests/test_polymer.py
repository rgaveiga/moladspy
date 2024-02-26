import sys

# TODO: Make this computer dettached
path = "E:\\PBC\\GitHub\\moladspy\\MolAdsPy"

sys.path.append(path)
from core import Atom
from core import Polymer
from numpy import array
import os

# TODO: Make a virtual simple polymer and test align, rotate and resize
# TODO: Test same tests on PLA
# TODO: Test same tests on Polypropilene
# TODO: Make Assert Files on results and Unit Tests

filename = "pla1"
filepath = os.path.join(path,filename + ".xyz")
name_polymer = "pla"
print("1)")

pol=Polymer("pla",filepath,vaccuum=20.0)

pol.write_xyz(filename + ".1.xyz")

print("\n9)")

pol2=pol.copy()

pol2.rotate(-30,-30,-30)
pol.rotate(30,30,30)
pol.write_xyz(filename + ".2.xyz")
pol2.write_xyz(filename + ".3.xyz")

print("Max ",pol.maxx,pol.maxy,pol.maxz)
print("Min ",pol.minx,pol.miny,pol.minz)
print("COM ",pol.center_of_mass)

pol.align('y')
pol.write_xyz(filename + ".align_y.xyz")

pol.align('x')
pol.write_xyz(filename + ".align_x.xyz")

pol.align('z')
pol.write_xyz(filename + ".align_z.xyz")



