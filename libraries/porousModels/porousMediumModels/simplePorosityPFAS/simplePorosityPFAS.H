/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Declaration of the simple porosity medium model (one kr / one pc law)
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

Namespace
    Foam::porousMediumModels::simplePorosityPFAS

Class
    Foam::simplePorosityPFAS

SourceFiles
    simplePorosityPFAS.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef simplePorosityPFAS_H
#define simplePorosityPFAS_H

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
             Class simplePorosityPFAS Declaration
\*---------------------------------------------------------------------------*/
class simplePorosityPFAS
:
    public twophasePorousMediumModel
{
protected:

    //- Disallow copy construct
    simplePorosityPFAS(const simplePorosityPFAS&);

    //- Disallow default bitwise assignment
    void operator=(const simplePorosityPFAS&);

public:

    //- Runtime type information
    TypeName("simplePorosityPFAS");

    // Constructors

    //- Construct from components
    simplePorosityPFAS
    (
        const word Sname,
        const fvMesh& mesh,
        const IOdictionary& transportProperties,
        const autoPtr<incompressiblePhase>& phase,
        const word porousRegion
    );

    //- Destructor
    ~simplePorosityPFAS()
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
