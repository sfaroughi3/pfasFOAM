pfasFOAM
========

**pfasFOAM** is an open-source OpenFOAM solver for the simulation of PFAS fate and transport in variably saturated (saturated/unsaturated) porous media. It is built on [porousMultiphaseFoam (PMF)](https://github.com/phorgue/porousMultiphaseFoam).

The link to publication:; 
- Baker, N.; Faroughi, S. *pfasFOAM: An OpenFOAM Solver for PFAS Fate and Transport in Variably Saturated Porous Media.* [doi:10.2139/ssrn.7358518](https://dx.doi.org/10.2139/ssrn.7358518)

About
-----

pfasFOAM is developed by Nicholas Baker and Salah Faroughi at the Energy & Intelligence Lab (EiLab), University of Utah. It is based on OpenFOAM 9 ([openfoam.org](https://openfoam.org) branch) and solves the governing equations for PFAS transport derived by Guo et al. (2020).

The toolbox provides the following solvers:

| Solver | Description |
|---|---|
| `pfasFOAM` | Main solver. Solves the governing equations with the empirical air-water interfacial area expression. |
| `pfasFOAMthermoaaw` | Uses the thermodynamic air-water interfacial area model found in newer formulations such as PFAS-LEACH. |

More details on the implementation will be available in the pfasFOAM paper (see [Citation](#citation)).

pfasFOAM is built on, and makes extensive use of, porousMultiphaseFoam, an excellent toolbox for subsurface modeling. The PMF contributors and the original PMF README are listed at the end of this document.

Requirements
------------

- OpenFOAM 9 (openfoam.org)

Installation
------------

Make sure the OpenFOAM 9 environment is loaded, then from the top-level directory run:

```bash
./Allwmake
```

This compiles the libraries, solvers and utilities into `$FOAM_USER_LIBBIN` and `$FOAM_USER_APPBIN`.

To clean the build:

```bash
./Allwclean            # remove build files
./Allwclean --purge    # also remove compiled binaries and libraries
```

Tutorials
---------

All examples from the pfasFOAM paper are provided in the `tutorials` directory:

- `1DExperiments`: Water flow validation, PFAS transport equation validation, and Comparison with PFAS-LEACH
- `2DExperiments`: 2D Water Flow, 2D Coupled Flow, and 2D Coupled Heterogeneous Flow
- `3DExperiments`: Mesh refinement study and Heterogeneous 3D case

To run a case, execute the `run` script in its directory; use `clean` to reset it. All tutorial cases run in serial by default.

```bash
cd tutorials/1DExperiments/PFASNonlinear/C1DPFAS
./run
```

Repository structure
--------------------

```
libraries/    porous media models (relative permeability, etc.), boundary conditions (including darcyGradPressure), numerical methods and GIS tools (from PMF)
solvers/      pfasFOAM solvers
utilities/    pre-processing utilities (from PMF)
tutorials/    example cases
tools/        helper scripts
extras/       Python Validation Code
```

Future work
-----------

Planned developments include:

- validation of the new coupling method
- support for multiple soil types
- run-time selectable air-water relationships
- a stable Newton method for the PFAS transport equation

Citation
--------

If you use pfasFOAM in a scientific publication, please cite:

- Baker, N.; Faroughi, S. *pfasFOAM: An OpenFOAM Solver for PFAS Fate and Transport in Variably Saturated Porous Media.* [doi:10.2139/ssrn.7358518](https://dx.doi.org/10.2139/ssrn.7358518)

License
-------

pfasFOAM is distributed under the GNU General Public License v3. See `COPYING_OPENFOAM` for details.

Acknowledgements
----------------

pfasFOAM is built on porousMultiphaseFoam. The original PMF README is reproduced below.

---

# Original README: porousMultiphaseFoam (PMF) for OpenFOAM

About
-----

**PMF** is an open-source toolbox dedicated to simulation of flow and transport processes in porous media.

It is based on the OpenFOAM environment and therefore benefits from its multiple feature such as pre-and post-processing tools or parallel efficiency for example.

Read the [PMF documentation](https://porousmultiphasefoam.readthedocs.io/en/dev/) for more information.

This toolbox needs only a standard OpenFOAM installation and is compatible with the two main forks: [www.openfoam.com](https://www.openfoam.com) and [www.openfoam.org](https://www.openfoam.org).

Documentation
-------------

For information about installation, usage or numerical details, you should read the [PMF documentation](https://porousmultiphasefoam.readthedocs.io/en/dev).

Improvements and bugs fixes are registered in the **ReleaseNotes.txt** file.

How to cite PMF
---------------

If your are using the PMF toolbox (directly or as a base-code) for a scientific publication, please cite at least one of the main references.

The initial paper:

- [An open-source toolbox for multiphase flow in porous media](https://www.sciencedirect.com/science/article/pii/S0010465514003403)
  P Horgue, C Soulaine, J Franc, R Guibert, G Debenest
  *Computer Physics Communications 187, 217-226*, **2015**

and/or the hydrological module:

- [porousMultiphaseFoam v2107: An open-source tool for modeling saturated/unsaturated water flows and solute transfers at watershed scale](https://www.sciencedirect.com/science/article/pii/S0010465521003908)
  P Horgue, F Renard, GS Gerlero, R Guibert, G Debenest
  *Computer Physics Communications 273, 108278*, **2022**

Developers
----------

Pierre Horgue (phorgue@imft.fr)  
Jacques Franc  
Gabriel Gerlero  
Gillien Latour  

Bug track/fix and beta tester
-----------------------------

avlmachado (github)  
Mohammed Elwardi Fadeli  
Romain Guibert  
Matteo Icardi  
Lachlan Mason  
François Renard  
Cyprien Soulaine
