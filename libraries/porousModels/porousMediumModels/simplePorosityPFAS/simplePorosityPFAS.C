/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Implementation of the simple porosity medium model (one kr / one pc law)
      used by the PFAS solver.

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

#include "simplePorosityPFAS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumModels
{
defineTypeNameAndDebug(simplePorosityPFAS, 0);

addToRunTimeSelectionTable
(
    twophasePorousMediumModel,
    simplePorosityPFAS,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumModels::simplePorosityPFAS::simplePorosityPFAS
(
    const word Sname,
    const fvMesh& mesh,
    const IOdictionary& transportProperties,
    const autoPtr<incompressiblePhase>& phase,
    const word porousRegion
)
    :
    twophasePorousMediumModel(Sname, mesh, transportProperties, phase, porousRegion)
{}

void Foam::porousMediumModels::simplePorosityPFAS::rewindTime()
{
    //- nothing for simple porosity model
}

void Foam::porousMediumModels::simplePorosityPFAS::correct()
{
    //- nothing for simple porosity model
}

void Foam::porousMediumModels::simplePorosityPFAS::correct(volScalarField& h, const bool steady, const bool conservative)
{
    //- nothing for simple porosity model
}
// ************************************************************************* //
