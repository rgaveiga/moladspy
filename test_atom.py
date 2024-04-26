import MolAdsPy
from MolAdsPy import Atom
import numpy as np

print("MolAdsPy version: %s" % MolAdsPy.__version__)

print("1)")
h_atom=Atom("H",x=[1,0,0])

print(h_atom)
print(h_atom.coords,h_atom.isfixed,h_atom.active)

print(Atom._instances)

print("\n2)")

h_atom2=h_atom.copy()
print(h_atom2)

print(h_atom.coords,h_atom.isfixed,h_atom.active)
print(h_atom2.coords,h_atom2.isfixed,h_atom.active)

h_atom2.coords[0]=0
h_atom2.coords[2]=-1
h_atom2.isfixed[2]=True
h_atom2.active=False

print(h_atom.coords,h_atom.isfixed,h_atom.active)
print(h_atom2.coords,h_atom2.isfixed,h_atom2.active)

print(Atom._instances)

print("\n3)")

he_atom=Atom("He")

print(he_atom)

he_atom.coords=(0,1,2)
he_atom.isfixed=np.array([True,True,False])
he_atom.active=False

print(he_atom.coords,he_atom.isfixed,he_atom.active)

print("\n4)")
print(Atom._curratomid)
print(Atom._instances)

print("\n5)")
print(h_atom.coords)

h_atom.displace([1,2,3])

print(h_atom.coords)

h_atom.coords[0]=-1

print(h_atom.coords)

print("\n6)")
print(h_atom2.isfixed)

h_atom2.isfixed[1]=True

print(h_atom2.isfixed)
print(h_atom2.oxidation_states[0],h_atom2.valency[0],h_atom2.valence_electrons)

print("\n7)")
os_atom=Atom("Os")

print(os_atom)
print(os_atom.element,os_atom.valence_electrons,os_atom.valency,
      os_atom.atomic_mass,os_atom.atomic_radius,os_atom.oxidation_states)

os_atom.coords[1]=-196

print(os_atom.coords)
