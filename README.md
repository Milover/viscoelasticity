## General

This is a somewhat modified port of the `viscoelasticFluidFoam` solver for
OpenFOAM v2106&ndash;v2112.

Both a steady state and a transient solver are available, along with
some rheological models (Newtonian, Carreau-Yasuda, PTT, multi-mode). The solver
assumes a laminar, isothermal flow of an incompressible fluid.

More stuff to be added at a later point.

## TODO

- [ ] add links and references
- [x] implement a transient solver
	- [ ] set up transient case template
- [ ] finish 3D mesh
	- [ ] 3D mesh independence test
- [ ] extract data from the [paper][ChahuanSasmal]
- [ ] set up automated testing for different versions

[ChahuanSasmal]: https://doi.org/10.1016/j.ijengsci.2021.103565
