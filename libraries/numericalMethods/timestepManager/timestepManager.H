/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / numericalMethods

  Modifications Copyright (C) 2026
      Energy & Intelligence Lab (Eilab) - University of Utah

  Authors of modifications: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Modifications (2026):
      - Declare printMaxLocation(), a debugging helper that prints the location of
        the cell with the global maximum derivative

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

Class
    Foam::fv::timestepManager

Description
    Timestep manager compute timestep based on temporal truncation error.
    Current version works with volScalarField and 3 time schemes:
    Euler / backward / CrankNicolson

SourceFiles
    timestepManager.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef timestepManager_H
#define timestepManager_H

#include "fvMesh.H"
#include "EulerD2dt2Scheme.H"
#include "EulerD3dt3Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Class timestepManager Declaration
\*---------------------------------------------------------------------------*/

class timestepManager
{
    
private:

    //- reference to time
    const Time& runTime_;

    //- reference to field
    const volScalarField& vf_;

    //- User defined truncation error
    const scalar truncationError_;

    //- reference to dry cells list
    const labelList* dryCells_;

    //- time scheme for this field
    const word timeScheme_;

    //- derivative operator
    EulerD3dt3Scheme d3dt3Operator_;
    fv::EulerD2dt2Scheme<scalar> d2dt2Operator_;

    //- maximal derivative/value
    scalar dVmax_, dV2max_, dV3max_;
    scalar Vmax_, V2max_, V3max_;

    void update1stOrder();
    void update2ndOrder();
    void update3rdOrder();
    void dryCellsDerivativeUpdate(volScalarField&);

    //- print (on master) the location of the cell holding the global max
    //  derivative; localMaxValue/maxCelli are this processor's local
    //  candidate, maxValue is the already-reduced global value
    void printMaxLocation
    (
        const string& derivName,
        const scalar maxValue,
        const scalar localMaxValue,
        const label maxCelli
    ) const;
 
public:

    // Constructors

    //- Construct from field
    timestepManager(
        const Time&, const volScalarField&, const scalar&,
        const labelList* dryCells = NULL
    );

    virtual ~timestepManager();

    //- Access functions
    scalar Vmax(){return Vmax_;}
    scalar dVmax(){return dVmax_;}
    scalar V2max(){return V2max_;}
    scalar dV2max(){return dV2max_;}
    scalar V3max(){return V3max_;}
    scalar dV3max(){return dV3max_;}

    //- Access functions
    void updateDerivatives();
    scalar computeTimestep(const scalar&);
    scalar computeTimestep();
    

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
