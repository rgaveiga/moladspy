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
import unittest
import os

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
class TestAdsorptionClass(unittest.TestCase):
    def test_load(self):
        graphene=Slab("graphene","graphene.xyz")
        tcnq=Molecule("tcnq","tcnq.xyz")

        self.assertIsNotNone(graphene,msg='Error loading graphene.xyz')
        self.assertIsNotNone(tcnq, msg='Error loading tcnq.xyz')

        self.assertEqual(graphene.label,"graphene",msg="Slab label is not graphene.")

        print(graphene)
        print(tcnq)
        print(graphene.label,len(graphene))
        print(tcnq.label,len(tcnq),tcnq.center_of_mass)

        print("\n2)")

    def test_roberto(self):
        graphene=Slab("graphene","graphene.xyz")
        tcnq=Molecule("tcnq","tcnq.xyz")

        adsorbed=Adsorption(graphene)
        top1_pos=array([graphene[0].coords[0],graphene[0].coords[1]])
        bridge1_pos=array([0.0,(graphene[1].coords[1]+(graphene[2].coords[1]-
                                                        graphene[1].coords[1])/2.0)])

        graphene.add_adsorption_site("top1",top1_pos)
        graphene.add_adsorption_site("bridge1",bridge1_pos)
        adsorbed.resize(4,4)

        print(len(adsorbed._atoms),len(adsorbed._species))

        adsorbed.add_component(tcnq,"top1",n=2,m=2,anchor="com",vertsep=2.0,side="top")

        print(len(adsorbed._atoms),len(adsorbed._species))
        print(tcnq.label,tcnq.center_of_mass)

        adsorbed.write_xyz("adsorbed.1.xyz")

        print("\n3)")

        tcnq.rotate(0,0,90)

        adsorbed.write_xyz("adsorbed.2.xyz")

        print("\n4)")

        print("Top: %f; Bottom: %f" % (adsorbed.top,adsorbed.bottom))

        tcnq2=tcnq.copy()

        adsorbed.add_component(tcnq2,[6.7,4.5],anchor="com",vertsep=2.5,side="bottom")

        print(len(adsorbed._atoms),len(adsorbed._species))
        print(tcnq2.label,tcnq2.center_of_mass)

        print("Top: %f; Bottom: %f" % (adsorbed.top,adsorbed.bottom))

        adsorbed.write_xyz("adsorbed.3.xyz")

        print("\n5)")

        tcnq2.displace([0,0,5])

        print(graphene.bottom-tcnq2.maxz)

        adsorbed.write_xyz("adsorbed.4.xyz")

        print("\n6)")

        tcnq2.rotate(90,90,0)

        print(graphene.bottom-tcnq2.maxz)

        adsorbed.write_xyz("adsorbed.5.xyz")

        print("\n7)")

        print(len(adsorbed._atoms))

        adsorbed.remove_component(tcnq2)

        print(len(adsorbed._atoms))

        h1=Atom("H",x=[0.5,0.0,0.0])
        h2=h1.copy()
        h2.coords[0]=1.5
        h2_mol=Molecule("h2")
        h2_mol.add_atom(h1)
        h2_mol.add_atom(h2)

        adsorbed.add_component(h2_mol,[6.7,4.5],anchor="com",vertsep=1.0,side="bottom")

        adsorbed.write_xyz("adsorbed.6.xyz")
        adsorbed.write_pw_input(pseudopotentials={"C":"C.paw.UPF","H":"H.paw.UPF"},
                                pwargs={"ibrav":8,"kvec":[4,4,4,0,0,0],
                                        "occupations":"smearing"})

        h2_mol.displace([-2,2,5])

        print(len(adsorbed._atoms))

        adsorbed.write_xyz("adsorbed.7.xyz")

        print("\n8)")

        adsorbed.set_separation(h2_mol,4)

        print(graphene.bottom-h2_mol.maxz)

        adsorbed.write_xyz("adsorbed.8.xyz")


if __name__=='__main__':
     unittest.main()