/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Modifications Copyright (C) 2026
      Energy & Intelligence Lab (Eilab) - University of Utah

  Authors of modifications: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Modifications (2026):
      - Construct lpower_ from field 'l<region>' or transportProperties
        (default 1.0)

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

#include "krVanGenuchten.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krVanGenuchten, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krVanGenuchten,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krVanGenuchten::krVanGenuchten
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& Sname,
    const word porousRegion
)
    :
    relativePermeabilityModel(mesh, transportProperties.subDict(typeName + "Coeffs"), Sname, porousRegion),
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
        dimensionedScalar(dimless, modelProperties_.lookupOrDefault<scalar>("m"+porousRegion,0))
    ),
    kramax_
    (
        IOobject
        (
            "kr"+Sname+"max"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.lookupOrDefault<scalar>("kr"+Sname+"max"+porousRegion,1.0))
    ),
    krbmax_
    (
        IOobject
        (
            "kr"+Sname+"max"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.lookupOrDefault<scalar>("kr"+Sname+"max"+porousRegion,1.0))
    ),
    lpower_
    (
        IOobject
        (
            "l"+porousRegion,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, modelProperties_.lookupOrDefault<scalar>("l"+porousRegion,1.0))
    )
{
    if (gMin(m_) <= 0)
    {
        FatalErrorIn
            (
                "in krVanGenuchten.C"
            )
            << "Relative permeability coefficient m equal or less than 0"
                << nl << "m" + porousRegion << " is required"
                << exit(FatalError);
    }
    Info << "Van Genuchten parameters for relative permeability model" << nl << "{" << endl;
    Info <<  "    " << Sname << porousRegion << "min" << " ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    " << Sname << porousRegion << "max" << " ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info << "    m" << porousRegion << " ";
    if (m_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(m_).value() << endl;}
    Info << "    kramax" << porousRegion << " ";
    if (kramax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(kramax_).value() << endl;}
    Info << "    krbmax" << porousRegion << " ";
    if (krbmax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(krbmax_).value() << endl;}
    Info << "} \n" << endl;
}

// ************************************************************************* //
