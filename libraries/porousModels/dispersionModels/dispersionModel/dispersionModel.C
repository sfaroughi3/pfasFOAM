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

#include "dispersionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dispersionModel, 0);
defineRunTimeSelectionTable(dispersionModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModel::dispersionModel
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
    :
    name_(name),
    transportProperties_(transportProperties),
    mesh_(mesh),
    Dm_(transportProperties.lookup("Dm")),
    Deff_
    (
        IOobject
        (
            name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionSet(0,2,-1,0,0,0,0)
    )
{}

// ************************************************************************* //
