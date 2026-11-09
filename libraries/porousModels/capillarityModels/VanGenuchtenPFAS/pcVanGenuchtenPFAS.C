/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Implementation of the Van Genuchten capillarity model with terms for PFAS
      modelling.

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

#include "pcVanGenuchtenPFAS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{
defineTypeNameAndDebug(pcVanGenuchtenPFAS, 0);

addToRunTimeSelectionTable
(
    capillarityModel,
    pcVanGenuchtenPFAS,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcVanGenuchtenPFAS::pcVanGenuchtenPFAS
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname,
    const word porousRegion
)
    :
    capillarityModel(mesh, transportProperties.subDict(typeName + "Coeffs"), Sname, porousRegion),
    m_
    (
        IOobject
        (
            "m"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, capillarityProperties_.lookupOrDefault<scalar>("m"+porousRegion, 0))
    ),
    n_(1/(1-m_)),
    alpha_ // necessary for Richards solver
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
        dimensionedScalar(dimless, capillarityProperties_.lookupOrDefault<scalar>("alpha"+porousRegion, 0))
    ),
    pc0_
    (
        IOobject
        (
            "pc0"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1,-1,-2,0,0), capillarityProperties_.lookupOrDefault<scalar>("pc0"+porousRegion,0))
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
    ),
    c_
    (
        IOobject
        (
            "c"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1,-3,0,0,0), capillarityProperties_.lookupOrDefault<scalar>("c"+porousRegion,0))
    )
{
    if (gMin(m_) == 0) FatalErrorIn("Foam::capillarityModels::pcVanGenuchtenPFAS::pcVanGenuchtenPFAS") << "m" << porousRegion << "=0 in pcVanGenuchtenPFAS" << abort(FatalError);
    Info << "Van Genuchten parameters for capillarity pressure model" << nl << "{" << endl;
    Info <<  "    " << Sname << porousRegion << "min" << " ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    " << Sname << porousRegion << "max" << " ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info << "    m" << porousRegion << " ";
    if (m_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(m_).value() << endl;}
    Info << "    pc0" << porousRegion << " ";
    if (pc0_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(pc0_).value() << endl;}
    Info << "    alpha" << porousRegion << " ";
    if (alpha_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alpha_).value() << endl;}
    Info << "    a" << porousRegion << " ";
    if (a_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(a_).value() << endl;}
    Info << "    b" << porousRegion << " ";
    if (b_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(b_).value() << endl;}
    Info << "    c" << porousRegion << " ";
    if (c_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(c_).value() << endl;}
    Info << "} \n" << endl;

}

// ************************************************************************* //
