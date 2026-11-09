/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Implementation of the Celia capillarity model (saturation and capillary
      capacity from head pressure).

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

#include "pcCelia.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{
defineTypeNameAndDebug(pcCelia, 0);

addToRunTimeSelectionTable
(
    capillarityModel,
    pcCelia,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcCelia::pcCelia
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname,
    const word porousRegion
)
    :
    capillarityModel(mesh, transportProperties.subDict(typeName + "Coeffs"), Sname, porousRegion),
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
        dimensionedScalar(dimless, capillarityProperties_.lookupOrDefault<scalar>("alpha"+porousRegion,0))
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
        dimensionedScalar(dimless, capillarityProperties_.lookupOrDefault<scalar>("beta"+porousRegion,0))
    )
{
    Info << "Celia parameters for capillarity pressure model" << nl << "{" << endl;
    Info <<  "    " << Sname << porousRegion << "min" << " ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    " << Sname << porousRegion << "max" << " ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info << "    alpha" << porousRegion << " ";
    if (alpha_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alpha_).value() << endl;}
    Info << "    beta" << porousRegion << " ";
    if (beta_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(beta_).value() << endl;}
    Info << "} \n" << endl;
}

// ************************************************************************* //
