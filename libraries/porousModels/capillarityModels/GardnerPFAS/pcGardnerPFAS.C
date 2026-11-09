/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Implementation of the Gardner capillarity model with Szyszkowski's
      equation to account for surfactant (PFAS) effects on capillarity pressure.

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

#include "pcGardnerPFAS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{
defineTypeNameAndDebug(pcGardnerPFAS, 0);

addToRunTimeSelectionTable
(
    capillarityModel,
    pcGardnerPFAS,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcGardnerPFAS::pcGardnerPFAS
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname,
    const word porousRegion
)
    :
    capillarityModel(mesh, transportProperties.subDict(typeName + "Coeffs"), Sname, porousRegion),
    alpha_G_
    (
        IOobject
        (
            "alpha_G"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,0,0,0,0), capillarityProperties_.lookupOrDefault<scalar>("alpha_G"+porousRegion,1))
    ),
    a_
    (
        IOobject
        (
            "a"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,-3,0,0,1), capillarityProperties_.lookupOrDefault<scalar>("a"+porousRegion,0))
    ),
    b_
    (
        IOobject
        (
            "b"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, capillarityProperties_.lookupOrDefault<scalar>("b"+porousRegion,0))
    )
{
    Info << "Exponential parameters for capillarity pressure model" << nl << "{" << endl;
    Info <<  "    " << Sname << porousRegion << "min" << " ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    " << Sname << porousRegion << "max" << " ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info << "    alpha_G" << porousRegion << " ";
    if (alpha_G_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alpha_G_).value() << endl;}
    Info << "    a" << porousRegion << " ";
    if (a_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(a_).value() << endl;}
    Info << "    b" << porousRegion << " ";
    if (b_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(b_).value() << endl;}
    Info << "} \n" << endl;
}

// ************************************************************************* //
