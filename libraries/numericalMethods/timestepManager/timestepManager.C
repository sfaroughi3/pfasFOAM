/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / numericalMethods

  Modifications Copyright (C) 2026
      Energy & Intelligence Lab (Eilab) - University of Utah

  Authors of modifications: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Modifications (2026):
      - Vmax_/V2max_/V3max_ now track max|field| over all cells instead of only the
        cell holding the maximum time derivative
      - Track the cell index of the maximum derivative and add a (currently
        commented-out) call to printMaxLocation() for debugging
      - Clamp the argument of the 2nd/3rd-order timestep estimates at zero before
        taking the root, avoiding NaN timesteps

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

#include "timestepManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //
Foam::timestepManager::timestepManager(
    const Time& runTime,
    const volScalarField& vf,
    const scalar& truncationError,
    const labelList* dryCells
)
    :
    runTime_(runTime),
    vf_(vf),
    truncationError_(truncationError),
    dryCells_(dryCells),
    timeScheme_(vf.mesh().ddtScheme("ddt("+vf.name()+")")),
    d3dt3Operator_(vf.mesh(),runTime.deltaTValue()),
    d2dt2Operator_(vf.mesh()),
    dVmax_(0),
    dV2max_(0),
    dV3max_(0),
    Vmax_(0),
    V2max_(0),
    V3max_(0)
{
    //- derivative initialization to keep 1st user-defined time step
    if (timeScheme_ == "backward")
    {
        V3max_ = max(SMALL, gMax(vf_.internalField()));
        dV3max_ = 3*truncationError_*(V3max_+VSMALL)/Foam::pow(runTime_.deltaTValue(),3);
    }
    else if  (timeScheme_ == "CrankNicolson")
    {
        V3max_ = max(SMALL, gMax(vf_.internalField()));
        dV3max_ = 12*truncationError_*(V3max_+VSMALL)/Foam::pow(runTime_.deltaTValue(),3);
    }
    else if (timeScheme_ == "Euler")
    {
        V2max_ = max(SMALL, gMax(vf_.internalField()));
        dV2max_ = 2*truncationError_*(V2max_+VSMALL)/Foam::pow(runTime_.deltaTValue(),2);
    }
    else if (timeScheme_ != "steadyState")
    {
        FatalErrorIn("timestepManager.C") << "time scheme " << timeScheme_
            << " does not work with timestepManager class" << exit(FatalError);
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timestepManager::~timestepManager()
{}

// * * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * //

void timestepManager::update1stOrder()
{
    volScalarField dVdT(vf_-vf_.oldTime());
    if (dryCells_) dryCellsDerivativeUpdate(dVdT);
    dVmax_ = 0;
    Vmax_ = 0;
    label maxCelli = -1;
    forAll(dVdT, celli)
    {
        if (mag(dVdT[celli]) > dVmax_)
        {
            dVmax_ = mag(dVdT[celli]);
            maxCelli = celli;
        }
        Vmax_ = max(Vmax_, mag(vf_[celli]));
    }
    scalar dVmaxLocal = dVmax_;
    reduce(dVmax_, maxOp<scalar>());
    reduce(Vmax_, maxOp<scalar>());
    if (Vmax_ == 0) Vmax_ = SMALL;

    //printMaxLocation("d("+vf_.name()+")/dt", dVmax_, dVmaxLocal, maxCelli);
}

void timestepManager::update2ndOrder()
{
    volScalarField dV2dT2(d2dt2Operator_.fvcD2dt2(vf_));
    if (dryCells_) dryCellsDerivativeUpdate(dV2dT2);
    dV2max_ = 0;
    V2max_ = 0;
    label maxCelli = -1;
    forAll(dV2dT2, celli)
    {
        if (mag(dV2dT2[celli]) > dV2max_)
        {
            dV2max_ = mag(dV2dT2[celli]);
            maxCelli = celli;
        }
        V2max_ = max(V2max_, mag(vf_[celli]));
    }
    scalar dV2maxLocal = dV2max_;
    reduce(dV2max_, maxOp<scalar>());
    reduce(V2max_, maxOp<scalar>());
    if (V2max_ == 0) V2max_ = SMALL;

    //printMaxLocation("d2("+vf_.name()+")/dt2", dV2max_, dV2maxLocal, maxCelli);
}

void timestepManager::update3rdOrder()
{
    volScalarField dV3dT3(d3dt3Operator_.fvcD3dt3(vf_));
    if (dryCells_) dryCellsDerivativeUpdate(dV3dT3);
    dV3max_ = 0;
    V3max_ = 0;
    label maxCelli = -1;
    forAll(dV3dT3, celli)
    {
        if (mag(dV3dT3[celli]) > dV3max_)
        {
            dV3max_ = mag(dV3dT3[celli]);
            maxCelli = celli;
        }
        V3max_ = max(V3max_, mag(vf_[celli]));
    }
    scalar dV3maxLocal = dV3max_;
    reduce(dV3max_, maxOp<scalar>());
    reduce(V3max_, maxOp<scalar>());
    if (V3max_ == 0) V3max_ = SMALL;

    //printMaxLocation("d3("+vf_.name()+")/dt3", dV3max_, dV3maxLocal, maxCelli);
}

//Debugging function that will show location of highest error
void timestepManager::printMaxLocation
(
    const string& derivName,
    const scalar maxValue,
    const scalar localMaxValue,
    const label maxCelli
) const
{
    //- identify the (single) processor owning the global max, then
    //  broadcast its cell location to all processors so that master
    //  (the only one whose Info<< is actually visible) can print it
    label owner =
        (maxCelli != -1 && localMaxValue == maxValue) ? Pstream::myProcNo() : -1;
    reduce(owner, maxOp<label>());
    if (owner == -1) return;

    point maxLoc(vector::zero);
    if (Pstream::myProcNo() == owner) maxLoc = vf_.mesh().C()[maxCelli];
    reduce(maxLoc, sumOp<vector>());

    Info<< "timestepManager: max(" << derivName << ") = " << maxValue
        << " at " << maxLoc << endl;
}

void timestepManager::dryCellsDerivativeUpdate(volScalarField& field)
{
    const labelList& dryCells = *dryCells_;
    forAll(dryCells, pointi) field[dryCells[pointi]] = 0;
}

// * * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * //

void timestepManager::updateDerivatives()
{
    update1stOrder();
    if (timeScheme_ == "Euler") update2ndOrder();
    else update3rdOrder();
}


scalar timestepManager::computeTimestep
(
    const scalar& truncationError
)
{
    scalar dt = GREAT;
    if (timeScheme_ == "Euler")
    {
        dt = Foam::pow(max(scalar(0), 2*truncationError*(V2max_+VSMALL)/(dV2max_+VSMALL)), 1./2.);
    }
    else if (timeScheme_ == "backward")
    {
        d3dt3Operator_.storeDeltaT00(runTime_.deltaT0Value());
        dt = Foam::pow(max(scalar(0), 3*truncationError*(V3max_+VSMALL)/(dV3max_+VSMALL)), 1./3.);
    }
    else if (timeScheme_ == "CrankNicolson")
    {
        d3dt3Operator_.storeDeltaT00(runTime_.deltaT0Value());
        dt = Foam::pow(max(scalar(0), 12*truncationError*(V3max_+VSMALL)/(dV3max_+VSMALL)), 1./3.);
    }
    reduce(dt, minOp<scalar>());
    return dt;
}

scalar timestepManager::computeTimestep()
{
    scalar dt = computeTimestep(truncationError_);
    return dt;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
