/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Implementation of the Celia relative permeability model,
      kr = A / (A + |100 h|^gamma).

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

#include "krCelia.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krCelia, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krCelia,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krCelia::krCelia
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname,
    const word porousRegion
)
    :
    relativePermeabilityModel(mesh, transportProperties.subDict(typeName + "Coeffs"), Sname, porousRegion),
    gammapow_
    (
        IOobject
        (
            "gamma"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.lookupOrDefault<scalar>("gamma"+porousRegion,0))
    ),
    Aparam_
    (
        IOobject
        (
            "A"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.lookupOrDefault<scalar>("A"+porousRegion,1))
    ),
    alpha_
    (
        IOobject
        (
            "alpha"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.lookupOrDefault<scalar>("alpha"+porousRegion,1))
    ),
    beta_
    (
        IOobject
        (
            "beta"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.lookupOrDefault<scalar>("beta"+porousRegion,1))
    )
{
    Info << "Celia parameters for relative permeability model" << nl << "{" << endl;
    Info <<  "    gamma" << porousRegion  << " ";
    if (gammapow_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(gammapow_).value() << endl;}
    Info << "    A" << porousRegion << " ";
    if (Aparam_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Aparam_).value() << endl;}
    Info << "} \n" << endl;
}

// ************************************************************************* //
