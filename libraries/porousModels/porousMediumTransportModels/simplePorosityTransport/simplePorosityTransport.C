/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Modifications Copyright (C) 2026
      Energy & Intelligence Lab (Eilab) - University of Utah

  Authors of modifications: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

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

#include "fvCFD.H"
#include "simplePorosityTransport.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumTransportModels
{
defineTypeNameAndDebug(simplePorosityTransport, 0);

addToRunTimeSelectionTable
(
    porousMediumTransportModel,
    simplePorosityTransport,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumTransportModels::simplePorosityTransport::simplePorosityTransport
(
    const word& phaseName,
    const porousMediumModel& pmModel
)
    :
    porousMediumTransportModel(phaseName, pmModel)
{}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

void Foam::porousMediumTransportModels::simplePorosityTransport::solveTransport
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& theta
)
{
    composition_.correct(U, theta);

    dictionary solverDict = pmModel_.mesh().solverDict("C");
    forAll(composition_.Y(), speciesi)
    {
        auto& C = composition_.Y(speciesi);
        const auto& R = composition_.R(speciesi);
        const auto& Deff = composition_.Deff(speciesi);
        const auto& lambda = composition_.lambda(speciesi);
        const auto& sourceTerm_tracer = composition_.sourceTerm(speciesi);
        const auto& xzero = pmModel_.transportProperties().lookup<dimensionedScalar>("xzero");
        fvScalarMatrix CEqn
            (
                R * fvm::ddt(theta,C)
                + fvm::div(phi, C, "div(phi,C)")
                - fvm::laplacian(theta*Deff, C, "laplacian(Deff,C)")
                ==
                - sourceTerm_tracer
                - R * theta * fvm::Sp(lambda,C)
            );

        CEqn.solve(solverDict);
    }
}

// ************************************************************************* //
