import MolAdsPy
from MolAdsPy import Atom
from MolAdsPy import Molecule
from numpy import array

print("MolAdsPy version: %s" % MolAdsPy.__version__)

print("1)")

h2_mol = Molecule("h2")

print(h2_mol)

for atom in h2_mol:
    print("Atom ", atom.ID, atom.coords)

print("\n2)")

h1 = Atom("H", x=[0.5, 0.0, 0.0])
h2 = h1.copy()
h2.coords[0] = 1.5

h2_mol.add_atom(h1)
h2_mol.add_atom(h2)
h2_mol.write_xyz("h2.1.xyz")

for atom in h2_mol:
    print("Atom ", atom.ID, atom.coords)

print("Max ", h2_mol.maxx, h2_mol.maxy, h2_mol.maxz)
print("Min ", h2_mol.minx, h2_mol.miny, h2_mol.minz)
print("COM ", h2_mol.center_of_mass)

print("\n3)")

h2_mol.rotate(90, 0, 0)
h2_mol.write_xyz("h2.2.xyz")

print("Max ", h2_mol.maxx, h2_mol.maxy, h2_mol.maxz)
print("Min ", h2_mol.minx, h2_mol.miny, h2_mol.minz)
print("COM ", h2_mol.center_of_mass)

print("\n4)")

h2_mol.add_anchor("top1", h2_mol[0].coords)
h2_mol.add_anchor("top2", h2_mol[1].coords)
h2_mol.displace([1, 1, 1])
h2_mol.write_xyz("h2.3.xyz")

for key in h2_mol.anchors.keys():
    print(key, h2_mol.anchors[key])

print("\n5)")

h2_mol.rotate(0, 90, 0, "top1")
h2_mol.write_xyz("h2.4.xyz")

print("Max ", h2_mol.maxx, h2_mol.maxy, h2_mol.maxz)
print("Min ", h2_mol.minx, h2_mol.miny, h2_mol.minz)
print("COM ", h2_mol.center_of_mass)

for key in h2_mol.anchors.keys():
    print(key, h2_mol.anchors[key])

print("\n6)")

h2_mol.rotate(0, 0, 90)
h2_mol.write_xyz("h2.5.xyz")

print("Max ", h2_mol.maxx, h2_mol.maxy, h2_mol.maxz)
print("Min ", h2_mol.minx, h2_mol.miny, h2_mol.minz)
print("COM ", h2_mol.center_of_mass)

for key in h2_mol.anchors.keys():
    print(key, h2_mol.anchors[key])

print("\n7)")

h2_mol.move_to([2, 2, 2], "top2")
h2_mol.write_xyz("h2.6.xyz")

print("Max ", h2_mol.maxx, h2_mol.maxy, h2_mol.maxz)
print("Min ", h2_mol.minx, h2_mol.miny, h2_mol.minz)
print("COM ", h2_mol.center_of_mass)

for key in h2_mol.anchors.keys():
    print(key, h2_mol.anchors[key])

print("\n8)")

tcnq = Molecule("tcnq", "tcnq.xyz", vaccuum=20.0)

tcnq.write_xyz("tcnq.1.xyz")

print("\n9)")

tcnq2 = tcnq.copy()

tcnq2.rotate(-30, -30, -30)
tcnq.rotate(30, 30, 30)
tcnq.write_xyz("tcnq.2.xyz")
tcnq2.write_xyz("tcnq.3.xyz")

print("Max ", tcnq.maxx, tcnq.maxy, tcnq.maxz)
print("Min ", tcnq.minx, tcnq.miny, tcnq.minz)
print("COM ", tcnq.center_of_mass)

print("\n10)")

print(tcnq2)

tcnq2.add_anchor("com+disp", tcnq2.anchors["com"] + array([1, 1, 2]))

for anchor in tcnq2.anchors:
    print(anchor, tcnq2.anchors[anchor])

tcnq2.rotate(10, -70, 90, "com+disp")

tcnq2.write_xyz("tcnq.4.xyz")

for anchor in tcnq2.anchors:
    print(anchor, tcnq2.anchors[anchor])

tcnq2.write_pw_input("tcnq.in")
