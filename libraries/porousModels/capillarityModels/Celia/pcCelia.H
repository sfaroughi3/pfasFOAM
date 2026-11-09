/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Declaration of the Celia capillarity model (saturation and capillary
      capacity from head pressure).

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
    Foam::capillarityModels::Celia

SourceFiles
    pcCelia.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef pcCelia_H
#define pcCelia_H

#include "capillarityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
    Class pcCelia Declaration
    \*---------------------------------------------------------------------------*/

class pcCelia
:
    public capillarityModel
{
    // pcVan Genuchten coefficients
    volScalarField alpha_;
    volScalarField beta_;


public:

    //- Runtime type information
    TypeName("Celia");

    // Constructors

    //- Construct from components
    pcCelia
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~pcCelia()
        {}

    // Member Functions

    //- Correct the capillarity
    void correct(const volScalarField& Sb)
        {
            Se_ = (Sb - Smin_)/(Smax_-Smin_);
        }

    //- Compute capillarity capacity from head pressure and return moisture (Richards' model)
    tmp<volScalarField> correctAndSb(volScalarField& h)
        {
            FatalErrorIn("volScalarField correctAndSb(volScalarField h) const")
                << "Not implemented in pcCelia.H " << abort(FatalError);
            return tmp<volScalarField>();
        };

    //- Compute capillarity capacity from head pressure and return moisture (Richards' model)
    scalar correctAndSb(volScalarField& h, const label& celli)
        {
            FatalErrorIn("scalar correctAndSb(const volScalarField& h, const label& celli)")
                << "Not implemented in pcCelia.H " << abort(FatalError);
            return 0;
        };
    //- PFAS version of Sb. To allow compatibility with PFAS solver for testing.
    tmp<volScalarField> correctAndSb(volScalarField& h, volScalarField& c)
        {
            // Update effective water saturation
            h.dimensions().reset(dimless);
            Se_ = alpha_ / (alpha_ + pow(mag(h) * 100, beta_));;
            // Update capillarity capacity
            Ch_.dimensions().reset(dimless);
            Ch_ = beta_ *  Se_ * pow(100, beta_) * (pow(mag(h), beta_-1))  * (Smax_ - Smin_) / (alpha_ + pow(mag(h) * 100, beta_));
            h.dimensions().reset(dimLength);
            Ch_.dimensions().reset(dimless/dimLength);

            return Smin_+(Smax_-Smin_)*Se_;
        }

    //- Compute saturation from given head pressure
    tmp<volScalarField> S(volScalarField& h)
        {
            FatalErrorIn("tmp<volScalarField> S(const volScalarField& h)")
                << "Not implemented in pcCelia.H " << abort(FatalError);
            return tmp<volScalarField>();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
