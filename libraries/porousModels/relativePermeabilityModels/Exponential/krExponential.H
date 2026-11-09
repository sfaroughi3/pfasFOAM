/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Declaration of the exponential relative permeability model,
      kr = exp(alpha_G h).

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
    Foam::relativePermeabilityModels::krExponential

SourceFiles
    krExponential.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef krExponential_H
#define krExponential_H

#include "relativePermeabilityModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
    Class krExponential Declaration
    \*---------------------------------------------------------------------------*/

class krExponential
:
    public relativePermeabilityModel
{
    //- Exponential coefficients
    volScalarField alpha_G_;

    //- cut-off point saturation (computed with entry pressure)

public:

    //- Runtime type //Information
    TypeName("rExponential");

    // Constructors

    //- Construct from components
    krExponential
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~krExponential()
        {}

    // Member Functions

    //- Correct the relative permeabilities
     void correct(const volScalarField& Sb)
        {
            FatalErrorIn("krExponential.H") << " Exponential model can be used only with groundwater solvers" << abort(FatalError);
        }
    void correctkra(const volScalarField& Sb)
        {
            FatalErrorIn("krExponential.H") << " Exponential model can be used only with groundwater solvers" << abort(FatalError);
        }
    void correctkrb(const volScalarField& Sb)
        {
            volScalarField h = Sb;
            dkrbdS_.dimensions().reset(dimless);
            h.dimensions().reset(dimless);
            krb_ = Foam::exp(alpha_G_ * (min(h, -VSMALL) + VSMALL));
            //- face-based kr, evaluated at the face-interpolated head
            const surfaceScalarField hface(fvc::interpolate(h));
            krbhave_ = Foam::exp(fvc::interpolate(alpha_G_) * (min(hface, -VSMALL) + VSMALL));
            dkrbdS_ =  alpha_G_ * dkrbdS_;
            h.dimensions().reset(dimLength);
        }

    void correctkrb(const surfaceScalarField& Sb, const volScalarField& Sb2)
        {
            FatalErrorIn("krExponential.H") << "Exponential model can be used only in test case with groundwatertransportPFAS" << abort(FatalError);
        }

    void correctkrb(const volScalarField& Sb, const label& celli)
        {
            FatalErrorIn("krExponential.H") << "Exponential model can be used only in test case with groundwatertransportPFAS" << abort(FatalError);
        }
    // void correctkrb(const volScalarField& h)
    //     {
    //         krb_ = Ks_ * Aparam_ / (Aparam_ + pow(pos(h), gammapow_));
    //     }
    tmp<volScalarField> kr(const volScalarField& S)
        {
            FatalErrorIn("krExponential.H") << "Exponential model can be used only in test case with groundwatertransportPFAS" << abort(FatalError);
            // volScalarField Se((S-Smin_)/(Smax_-Smin_));
            // return pow(Se,1) * pow(
            //     (1-pow(1-pow(Sc_*Se,1/m_),m_)) /
            //     (1-pow(1-pow(Sc_,1/m_),m_))
            //     ,2);
            return tmp<volScalarField>();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
