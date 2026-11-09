/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Implementation of a dispersion model with constant diffusion coefficients
      Dx, Dy and Dz read from the constantDiffusionCoeffs dictionary.

  License:
      SPDX-License-Identifier: GPL-3.0-or-later

      This file is part of pfasFOAM.

      OpenFOAM is free software: you can redistribute it and/or modify it
      under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      (at your option) any later version.

      OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
      ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
      FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
      for more details.

      You should have received a copy of the GNU General Public License
      along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#include "constantDiffusion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
defineTypeNameAndDebug(constantDiffusion, 0);

addToRunTimeSelectionTable
(
    dispersionModel,
    constantDiffusion,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModels::constantDiffusion::constantDiffusion
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
    :
    dispersionModel(name, transportProperties, mesh),
    constantDiffusionCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    Dx_
    (
        IOobject
        (
            "Dx",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(constantDiffusionCoeffs_.lookup("Dx"))
    ),
    Dy_
    (
        IOobject
        (
            "Dy",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(constantDiffusionCoeffs_.lookup("Dy"))
    ),
    Dz_
    (
        IOobject
        (
            "Dz",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(constantDiffusionCoeffs_.lookup("Dz"))
    )
{
    Info << "Parameters for constant diffusion model" << nl << "{" << endl;
    Info << "    Dx ";
    if (Dx_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Dx_).value() << endl;}
    Info <<  "   Dy ";
    if (Dy_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Dy_).value() << endl;}
    Info << "    Dz ";
    if (Dz_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Dz_).value() << endl;}
    Info << "}" << endl;
}

// ************************************************************************* //
