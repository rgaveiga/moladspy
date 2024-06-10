from MolAdsPy import Atom
from MolAdsPy import Molecule
from MolAdsPy import Slab
from MolAdsPy import Adsorption
from numpy import array

print("1)")

graphene = Slab("graphene", "graphene.xyz")
tcnq = Molecule("tcnq", "tcnq.xyz")

print(graphene)
print(tcnq)
print(graphene.label, len(graphene))
print(tcnq.label, len(tcnq), tcnq.center_of_mass)

print("\n2)")

adsorbed = Adsorption(graphene)
top1_pos = array([graphene[0].coords[0], graphene[0].coords[1]])
bridge1_pos = array(
    [
        0.0,
        (graphene[1].coords[1] + (graphene[2].coords[1] - graphene[1].coords[1]) / 2.0),
    ]
)

graphene.add_adsorption_site("top1", top1_pos)
graphene.add_adsorption_site("bridge1", bridge1_pos)

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))

print("\n3)")

graphene.resize(4, 4)

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)

print("\n4)")

adsorbed.add_component(tcnq, "top1", n=2, m=2, anchor="com", vertsep=2.0, side="top")

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))
print(tcnq.label, tcnq.center_of_mass)

adsorbed.write_xyz("adsorbed.1.xyz")

print("\n5)")

tcnq.rotate(0, 0, 90)

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))
print(tcnq.label, tcnq.center_of_mass)

adsorbed.write_xyz("adsorbed.2.xyz")

print("\n6)")

tcnq2 = tcnq.copy()

adsorbed.add_component(tcnq2, [6.7, 4.5], anchor="com", vertsep=2.5, side="bottom")

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))
print(tcnq2.label, tcnq2.center_of_mass)

adsorbed.write_xyz("adsorbed.3.xyz")

print("\n7)")

tcnq2.displace([0, 0, 5])

print(graphene.bottom - tcnq2.maxz)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))

adsorbed.write_xyz("adsorbed.4.xyz")

print("\n8)")

tcnq2.rotate(90, 90, 0)

print(graphene.bottom - tcnq2.maxz)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))

adsorbed.write_xyz("adsorbed.5.xyz")

print("\n9)")

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)

adsorbed.remove_component(tcnq2)

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))

print("\n10)")

h1 = Atom("H", x=[0.5, 0.0, 0.0])
h2 = h1.copy()
h2.coords[0] = 1.5
h2_mol = Molecule("h2")
h2_mol.add_atom(h1)
h2_mol.add_atom(h2)

adsorbed.add_component(h2_mol, [6.7, 4.5], anchor="com", vertsep=1.0, side="bottom")

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))

adsorbed.write_xyz("adsorbed.6.xyz")
adsorbed.write_pw_input(
    pseudopotentials={"C": "C.paw.UPF", "H": "H.paw.UPF"},
    pwargs={"ibrav": 8, "kvec": [4, 4, 4, 0, 0, 0], "occupations": "smearing"},
)

print("\n11)")

h2_mol.displace([-2, 2, 5])

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))

adsorbed.write_xyz("adsorbed.7.xyz")

print("\n12)")

adsorbed.set_separation(h2_mol, 4)

print(graphene.bottom - h2_mol.maxz)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))

adsorbed.write_xyz("adsorbed.8.xyz")

print("\n13)")

adsorbed.remove_component(tcnq)

print(
    "# of atoms, # of species: %d %d" % (len(adsorbed._atoms), len(adsorbed._species))
)
print("Top: %f; Bottom: %f" % (adsorbed.top, adsorbed.bottom))

adsorbed.write_xyz("adsorbed.9.xyz")
