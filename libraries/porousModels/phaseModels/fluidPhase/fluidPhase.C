/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Redistributed as part of pfasFOAM by:
      Energy & Intelligence Lab (Eilab) - University of Utah
      (Nicholas T. Baker & Salah A. Faroughi), 2026

  Modifications (2026):
      - None to the code of this file; it is unchanged from the upstream
        source (apart from this header).

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

#include "fluidPhase.H"
#include "fixedValueFvPatchFields.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPhase::fluidPhase
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
    :  
    phase(mesh, transportProperties, phaseName),
    U_
    (
        IOobject
        (
            "U" + phaseName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

    phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phi"+phaseName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                Foam::linearInterpolate(U_) & mesh.Sf(),
                phiTypes
            )
        );
    if (phiPtr_->headerOk()) Info << nl << "Reading field phi" << phaseName << endl;
    else Info<< nl << "Computing field phi" << phaseName << " from field U" << phaseName << endl;
}


Foam::autoPtr<Foam::fluidPhase> Foam::fluidPhase::New
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
{
    return autoPtr<fluidPhase>
    (
        new fluidPhase(mesh, transportProperties, phaseName)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPhase::~fluidPhase()
{}

// ************************************************************************* //
