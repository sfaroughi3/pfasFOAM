/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________
 		 |_   __ \|_   \  /   _||_   __  |
   		   | |__) | |   \/   |    | |_ \_|
   		   |  ___/  | |\  /| |    |  _|
    		  _| |_    _| |_\/_| |_  _| |_
   		 |_____|  |_____||_____||_____|
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
    flows through porous media.

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

\*---------------------------------------------------------------------------*/

#include "GuoDispersion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
defineTypeNameAndDebug(GuoDispersion, 0);

addToRunTimeSelectionTable
(
    dispersionModel,
    GuoDispersion,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModels::GuoDispersion::GuoDispersion
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
    :
    dispersionModel(name, transportProperties, mesh),
    GuoDispersionCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    alphaL_
    (
        IOobject
        (
            "alphaL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(GuoDispersionCoeffs_.lookup("alphaL"))
    ),
    alphaT_
    (
        IOobject
        (
            "alphaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(GuoDispersionCoeffs_.lookup("alphaT"))
    ),
    theta_s_
    (
        IOobject
        (
            "theta_s",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(GuoDispersionCoeffs_.lookup("theta_s"))
    )
{
    Info << "Parameters for alpha dispersion model" << nl << "{" << endl;
    Info << "    theta_s ";
    if (theta_s_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(theta_s_).value() << endl;}
    Info <<  "    alphaL ";
    if (alphaL_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alphaL_).value() << endl;}
    Info << "    alphaT ";
    if (alphaT_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alphaT_).value() << endl;}
    Info << "}" << endl;
}

// ************************************************************************* //
