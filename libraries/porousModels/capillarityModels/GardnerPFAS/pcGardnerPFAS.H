/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Declaration of the Gardner capillarity model with Szyszkowski's equation
      to account for surfactant (PFAS) effects on capillarity pressure.

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
    Foam::capillarityModels::GardnerPFAS

SourceFiles
    pcGardnerPFAS.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef pcGardnerPFAS_H
#define pcGardnerPFAS_H

#include "capillarityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
    Class pcGardnerPFAS Declaration
    \*---------------------------------------------------------------------------*/

class pcGardnerPFAS
:
    public capillarityModel
{
    // pcVan Genuchten coefficients
    volScalarField alpha_G_;
    volScalarField a_, b_;

public:

    //- Runtime type information
    TypeName("GardnerPFAS");

    // Constructors

    //- Construct from components
    pcGardnerPFAS
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~pcGardnerPFAS()
        {}

    // Member Functions

    //- Correct the capillarity
    void correct(const volScalarField& Sb)
        {
            return;
        }

    //- Compute capillarity capacity from head pressure and return moisture (Richards' model)
    tmp<volScalarField> correctAndSb(volScalarField& h)
        {
            FatalErrorIn("volScalarField correctAndSb(volScalarField h) const")
                << "Not implemented in pcGardnerPFAS.H " << abort(FatalError);
            return tmp<volScalarField>();
        };

    //- Compute capillarity capacity from head pressure and return moisture (Richards' model)
    scalar correctAndSb(volScalarField& h, const label& celli)
        {
            FatalErrorIn("scalar correctAndSb(const volScalarField& h, const label& celli)")
                << "Not implemented in pcGardnerPFAS.H " << abort(FatalError);
            return 0;
        };
    //- PFAS equation. Not valid here, thus ignore
    tmp<volScalarField> correctAndSb(volScalarField& h, volScalarField& c)
        {
            volScalarField tempvar = (c / a_) + 1;
            volScalarField conversion_ = alpha_G_/(1-b_*log(tempvar));

            h.dimensions().reset(dimless);
            Se_ =  Foam::exp(conversion_* (min(h, -VSMALL) + VSMALL));
            h.dimensions().reset(dimLength);

            Ch_.dimensions().reset(dimless);
            Ch_ = (Smax_-Smin_) * conversion_ * Se_ ;
            Ch_.dimensions().reset(dimless/dimLength);

            return Smin_+(Smax_-Smin_)*Se_;
        }

    //- Compute saturation from given head pressure
    tmp<volScalarField> S(volScalarField& h)
        {
            FatalErrorIn("tmp<volScalarField> S(const volScalarField& h)")
                << "Not implemented in pcGardnerPFAS.H " << abort(FatalError);
            return tmp<volScalarField>();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
