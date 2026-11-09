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
    Foam::porousMediumTransportModels::simplePorosityTransport

Class
    Foam::simplePorosityTransport

Description
    Simple porosity model for solute transport

SourceFiles
    simplePorosityTransport.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef simplePorosityTransport_H
#define simplePorosityTransport_H

#include "volFields.H"
#include "porousMediumTransportModel.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"
#include "incompressiblePhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumTransportModels
{

/*---------------------------------------------------------------------------* \
             Class simplePorosityTransport Declaration
\*---------------------------------------------------------------------------*/

class simplePorosityTransport
:
    public porousMediumTransportModel
{
protected:

    //- Disallow copy constructn
    simplePorosityTransport(const simplePorosityTransport&);

    //- Disallow default bitwise assignment
    void operator=(const simplePorosityTransport&);

public:

    //- Runtime type information
    TypeName("simplePorosityTransport");

    // Constructors

    //- Construct from components
    simplePorosityTransport
    (
        const word& phaseName,
        const porousMediumModel& pmModel
    );

    //- Destructor
    ~simplePorosityTransport()
        {}

    //- Solve tracer tranport equation
    void solveTransport
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& theta
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousMediumTransportModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
