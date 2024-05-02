from MolAdsPy import Atom
from MolAdsPy import Slab
import numpy as np

print("1)")
graphene = Slab("graphene", "graphene.xyz", file_type="xyz", vaccuum=16)

print(graphene)

for atom in graphene:
    print(atom.ID, atom.symbol, atom.element, atom.coords)

graphene.write_xyz("graphene.1.xyz")

print("\n2)")

graphene2 = graphene.copy()

print(graphene2)

for atom in graphene2:
    print(atom.ID, atom.symbol, atom.element, atom.coords)

print("\n3)")
graphene.resize(2, 2)

graphene.write_xyz("graphene.2.xyz")

print("\n4)")
graphene.resize(1, 1)

print("Total number of atoms: %d" % len(graphene))
print("Number of active atoms only: %d" % len(graphene.active_atoms))

graphene.write_xyz("graphene.3.xyz")
graphene.write_pw_input("graphene_ucell.in", ucell=True)
graphene.write_pw_input("graphene.1.in")

print("\n5)")

graphene2.resize(3, 3)

graphene2.write_xyz("graphene.4.xyz")

print("\n6)")

print(
    "Total number of active atoms before modification: %d" % len(graphene2.active_atoms)
)

graphene2[8].active = False
graphene2[24].active = False
graphene2.remove_atom(66)

print(
    "Total number of active atoms after modification: %d" % len(graphene2.active_atoms)
)

graphene2.write_xyz("graphene.5.xyz")
graphene2.write_pw_input(
    "graphene2.1.in", pseudopotentials={"C": "C.paw.UPF", "H": "H.paw.UPF"}
)

print("\n7)")

graphene2.resize(2, 3)

print(
    "Total number of active atoms before modification: %d" % len(graphene2.active_atoms)
)

graphene2[8].active = True

print(
    "Total number of active atoms before modification: %d" % len(graphene2.active_atoms)
)

graphene2.write_xyz("graphene.6.xyz")

print("\n8)")

S = Atom("S")
S.coords = graphene2[24].coords
S.coords[2] = 1.5

graphene2.add_atom(S, loc=graphene2.location(graphene2[24]))
graphene2[0].displace([0, -0.1, 0])
graphene2[3].displace([0, 0, -1])
graphene2.write_xyz("graphene.7.xyz")

print(graphene2.origin, graphene2.top, graphene2.bottom)

print("\n9)")

for atom in graphene2.unit_cell:
    print(atom.ID, atom.coords)

graphene2.a0 = 1.2

for atom in graphene2.unit_cell:
    print(atom.ID, atom.coords)

graphene2.write_xyz("graphene.8.xyz", ucell=True)

print("\n10)")

graphene2.a0 = 1.0

for atom in graphene2.unit_cell:
    print(atom.ID, atom.coords)

graphene2.write_xyz("graphene.9.xyz", ucell=True)

print("\n11)")

top1_pos = np.array([graphene[0].coords[0], graphene[0].coords[1]])
bridge1_pos = np.array(
    [
        0.0,
        (graphene[1].coords[1] + (graphene[2].coords[1] - graphene[1].coords[1]) / 2.0),
    ]
)

graphene.add_adsorption_site("top1", top1_pos)
graphene.add_adsorption_site("bridge1", bridge1_pos)

h1 = Atom(
    "H",
    x=np.array(
        [
            graphene.adsorption_sites["top1"][0],
            graphene.adsorption_sites["top1"][1],
            1.0,
        ]
    ),
)
h2 = Atom(
    "H",
    x=np.array(
        [
            graphene.adsorption_sites["bridge1"][0],
            graphene.adsorption_sites["bridge1"][1],
            1.0,
        ]
    ),
)

graphene.add_atom(h1)
graphene.add_atom(h2)

graphene.write_xyz("graphene.10.xyz")

print(graphene.adsorption_sites)
print(graphene2.adsorption_sites)

print("\n12)")

graphene2.write_pw_input("graphene2.2.in", pwargs={"ecutwfc": 48, "ecutrho": 480})
graphene.write_pw_input("graphene.2.in", pwargs={"kvec": [4, 4, 1, 0, 0, 0]})

print("\n13)")

graphene2.resize(5, 0)
graphene2.write_xyz("graphene.11.xyz")

print("\n14)")

print(graphene2.lattice_vectors)

graphene2.resize(3, 3)
graphene2.write_pw_input("graphene2.3.in")

graphene2.lattice_vectors *= 1.01

print(graphene2.lattice_vectors)

graphene2.write_pw_input("graphene2.4.in")

print("\n15)")

print(graphene.origin)
graphene.write_xyz("graphene.12.xyz")

graphene.origin = [-5, -10, 10]

print(graphene.origin)
graphene.write_xyz("graphene.13.xyz")

print("\n16)")

b1 = Atom("B")
b2 = b1.copy()
n1 = Atom("N")
n2 = n1.copy()
b1.coords = graphene[0].coords
n1.coords = graphene[1].coords
b2.coords = graphene[2].coords
n2.coords = graphene[3].coords
atomlist = [b1, n1, b2, n2]
bn = Slab("boron nitride", atomlist, graphene.a0, graphene.lattice_vectors)
bn.origin = [0.0, 0.0, 0.0]
bn.a0 = 1.0162

bn.write_xyz("bn.1.xyz")

print("\n17)")

bn.resize(3, 3)

bn.write_xyz("bn.2.xyz")

for atom in bn:
    if atom.symbol == "N":
        atom.active = False

bn.write_xyz("bn.3.xyz")
