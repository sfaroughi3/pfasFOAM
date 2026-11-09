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

Namespace
    Foam::porousMediumModels::simplePorosity

Class
    Foam::simplePorosity

Description
    Simple porosity model (one kr / one pc law)

SourceFiles
    simplePorosity.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef simplePorosity_H
#define simplePorosity_H

#include "volFields.H"
#include "twophasePorousMediumModel.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"
#include "incompressiblePhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumModels
{

/*---------------------------------------------------------------------------* \
             Class simplePorosity Declaration
\*---------------------------------------------------------------------------*/
class simplePorosity
:
    public twophasePorousMediumModel
{
protected:

    //- Disallow copy construct
    simplePorosity(const simplePorosity&);

    //- Disallow default bitwise assignment
    void operator=(const simplePorosity&);

public:

    //- Runtime type information
    TypeName("simplePorosity");

    // Constructors

    //- Construct from components
    simplePorosity
    (
        const word Sname,
        const fvMesh& mesh,
        const IOdictionary& transportProperties,
        const autoPtr<incompressiblePhase>& phase,
        const word porousRegion
    );

    //- Destructor
    ~simplePorosity()
        {}

    //- Correct the relative permeabilities
    void rewindTime();
    void correct();
    void correct(volScalarField& h, const bool steady = false, const bool conservative = true);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousMediumModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
