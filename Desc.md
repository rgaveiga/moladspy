# MolAdsPy

This is a package written entirely in Python to place and manipulate molecules on a substrate and then generate coordinates that can be used as inputs to first-principles calculation codes like [Quantum Espresso](https://www.quantum-espresso.org/) or [VASP](https://www.vasp.at/).

In more specific terms, the package allows the user to create objects such as molecules and slabs, which are made of atomic objects, and combine them into a system in which one or more molecular objects are adsorbed onto a slab object that acts as a substrate. An adsorbed molecule can be moved to a specific location on the substrate (e.g., a predefined adsorption site), can be shifted a certain distance from its current position, and can be rotated in three dimensions. Last but not least, an adsorbed molecule can be desorbed (removed) from the substrate.

Both the source code and documentation, whenever they are updated, will be available on the [project's GitHub site](https://github.com/rgaveiga/moladspy). Particularly, the package documentation is still at a very early stage. It's something that takes a lot of time and time ends up being a precious commodity that I don't have much of.

> In any case, in the project's *tests* directory on GitHub, use cases that show the library's capabilities are made available in the form of Jupyter notebooks.

Since the adsorbed structures are built algorithmically from Python scripts, hundreds or thousands of them can be created in a fraction of a second. This gives the modeler a lot of flexibility, especially if he/she is working on the high-throughput screening of this kind of system. That being said, the code is provided as is. The author makes no guarantee that its results are accurate and is not responsible for any losses caused by the use of the code.

If you use *MolAdsPy* in your own projects and if you don't mind, when publishing the results, please cite the following papers of mine:

* [*Structural, electronic, and magnetic properties of pristine and oxygen-adsorbed graphene nanoribbons*](https://www.sciencedirect.com/science/article/abs/pii/S0169433210003855)
* [*Interfacial properties of polyethylene/Ti3C2Tx mxene nanocomposites investigated by first-principles calculations*](https://www.sciencedirect.com/science/article/abs/pii/S0169433222028720)
* [*Self-assembly of NiTPP on Cu(111): a transition from disordered 1D wires to 2D chiral domains*](https://pubs.rsc.org/en/content/articlelanding/2015/cp/c5cp01288k)

If you have any comments, suggestions or corrections, please feel free to [drop me a message](mailto:roberto.veiga@ufabc.edu.br).