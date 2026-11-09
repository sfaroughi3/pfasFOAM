/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Modifications Copyright (C) 2026
      Energy & Intelligence Lab (Eilab) - University of Utah

  Authors of modifications: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Modifications (2026):
      - Add correctkrb(Sb, Sb2): cell kr and dkr/dS from the water content Sb2
        (effective saturation clamped away from 0 and 1) and face-based krbhave_
        from the face-interpolated water content Sb
      - Include fvc.H; keep the commented-out correctkrb(h) placeholder

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
    Foam::relativePermeabilityModels::krIppisch

Description
    Ippisch Vogel Bastian relative permeability model (AWR 2006)

SourceFiles
    krIppisch.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef krIppisch_H
#define krIppisch_H

#include "relativePermeabilityModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
    Class krIppisch Declaration
    \*---------------------------------------------------------------------------*/

class krIppisch
:
    public relativePermeabilityModel
{
    //- Ippisch coefficients
    volScalarField m_,n_,alpha_,tau_,he_;

    //- cut-off point saturation (computed with entry pressure)
    volScalarField Sc_;

public:

    //- Runtime type information
    TypeName("Ippisch");

    // Constructors

    //- Construct from components
    krIppisch
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~krIppisch()
        {}

    // Member Functions

    //- Correct the relative permeabilities
     void correct(const volScalarField& Sb)
        {
            FatalErrorIn("krIppisch.H") << " Ippisch model can be used only with groundwater solvers" << abort(FatalError);
        }
    void correctkra(const volScalarField& Sb)
        {
            FatalErrorIn("krIppisch.H") << " Ippisch model can be used only with groundwater solvers" << abort(FatalError);
        }
    void correctkrb(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            krb_ = pow(Se_,tau_) * pow(
                (1-pow(1-pow(Sc_*Se_,1/m_),m_)) /
                (1-pow(1-pow(Sc_,1/m_),m_))
                ,2);
        }
    //- Cell kr and dkr/dS from the water content Sb2, and face-based kr
    //  (krbhave_) from the face-interpolated water content Sb
    void correctkrb(const surfaceScalarField& Sb, const volScalarField& Sb2)
        {
            const surfaceScalarField Sminf(fvc::interpolate(Smin_));
            const surfaceScalarField Smaxf(fvc::interpolate(Smax_));
            const surfaceScalarField Sef
            (
                min(scalar(1), max(scalar(0), (Sb - Sminf)/(Smaxf - Sminf)))
            );
            const surfaceScalarField mf(fvc::interpolate(m_));
            const surfaceScalarField Scf(fvc::interpolate(Sc_));
            const surfaceScalarField Fcf(1 - pow(1 - pow(Scf,1/mf), mf));
            const surfaceScalarField Fxf(1 - pow(1 - pow(Scf*Sef,1/mf), mf));
            krbhave_ = pow(Sef, fvc::interpolate(tau_)) * pow(Fxf/Fcf, 2);

            Se_ = min(scalar(1) - SMALL, max(SMALL, (Sb2 - Smin_)/(Smax_ - Smin_)));
            const volScalarField ScSe(Sc_*Se_);
            const volScalarField Fx(1 - pow(1 - pow(ScSe,1/m_), m_));
            const volScalarField Fc(1 - pow(1 - pow(Sc_,1/m_), m_));
            //- derivative of F(x) = 1-(1-x^(1/m))^m evaluated at Sc*Se
            const volScalarField dFx(pow(1 - pow(ScSe,1/m_) + VSMALL, m_-1)*pow(ScSe,1/m_-1));
            krb_ = pow(Se_,tau_) * pow(Fx/Fc, 2);
            dkrbdS_ =
                (
                    tau_*pow(Se_,tau_-1)*pow(Fx/Fc, 2)
                  + 2*pow(Se_,tau_)*(Fx/Fc)*Sc_*dFx/Fc
                )/(Smax_ - Smin_);
        }
    void correctkrb(const volScalarField& Sb, const label& celli)
        {
            scalar Se = (Sb[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            krb_[celli] = pow(Se_[celli],tau_[celli]) * pow(
                (1-pow(1-pow(Sc_[celli]*Se,1/m_[celli]),m_[celli])) /
                (1-pow(1-pow(Sc_[celli],1/m_[celli]),m_[celli]))
                ,2);
        }

    // void correctkrb(const volScalarField& h)
    //     {
    //         FatalErrorIn("krIppisch.H") << "This is only for test cases in Celia" << abort(FatalError);
    //     }

    tmp<volScalarField> kr(const volScalarField& S)
        {
            volScalarField Se((S-Smin_)/(Smax_-Smin_));
            return pow(Se,tau_) * pow(
                (1-pow(1-pow(Sc_*Se,1/m_),m_)) /
                (1-pow(1-pow(Sc_,1/m_),m_))
                ,2);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
