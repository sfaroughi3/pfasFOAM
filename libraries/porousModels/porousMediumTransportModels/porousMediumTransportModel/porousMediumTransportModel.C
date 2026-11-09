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

#include "porousMediumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(porousMediumTransportModel, 0);
defineRunTimeSelectionTable(porousMediumTransportModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumTransportModel::porousMediumTransportModel
(
    const word& phaseName,
    const porousMediumModel& pmModel
)
    :
    transportProperties_(pmModel.transportProperties()),
    phaseName_(phaseName),
    speciesNames_(transportProperties_.lookupOrDefault("species", wordList(1, "C"))),
    pmModel_(pmModel),
    composition_(
        transportProperties_,
        speciesNames_,
        pmModel_.mesh(),
        word::null,
        pmModel_.eps(),
        &sourceEventList_,
        "C"
    )
{
    composition_.check_eps();
}

Foam::wordList Foam::porousMediumTransportModel::speciesNames(Foam::word porousRegion)
{
    wordList speciesNamesList = wordList(speciesNames_);
    forAll(speciesNamesList, wordi)
    {
        speciesNamesList[wordi] += porousRegion;
    }
    return speciesNamesList;
}

// ************************************************************************* //
