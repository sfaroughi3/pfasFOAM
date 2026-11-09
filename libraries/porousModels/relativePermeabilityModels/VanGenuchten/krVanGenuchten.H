/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Modifications Copyright (C) 2026
      Energy & Intelligence Lab (Eilab) - University of Utah

  Authors of modifications: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Modifications (2026):
      - Add lpower_ so the Mualem exponent is configurable (was fixed at 0.5)
      - Add correctkrb(Sb, Sb2) computing cell kr and dkr/dS (zero at Se >= 1)
        and the face-based krbhave_
      - Include fvc.H; keep alternative formulations as commented-out code

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
    Foam::relativePermeabilityModels::krVanGenuchten

Description
    Standard Van Genuchten relative permeability model.

SourceFiles
    krVanGenuchten.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef krVanGenuchten_H
#define krVanGenuchten_H

#include "relativePermeabilityModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
    Class krVanGenuchten Declaration
    \*---------------------------------------------------------------------------*/

class krVanGenuchten
:
    public relativePermeabilityModel
{
    //- Van Genuchten coefficients
    volScalarField m_;
    volScalarField lpower_;

    //- end points
    volScalarField kramax_;
    volScalarField krbmax_;

public:

    //- Runtime type information
    TypeName("VanGenuchten");

    // Constructors

    //- Construct from components
    krVanGenuchten
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~krVanGenuchten()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            kra_ = kramax_ * pow(1-Se_,0.5) * pow(1-pow(Se_,1/m_),2*m_);
            //krb_ = krbmax_ * pow(Se_,0.5) * pow(1-pow(1-pow(Se_,1/m_),m_),2);

            dkradS_ = - pow((1-pow(Se_,1/m_)+VSMALL),2*m_-1) * (-5*pow(Se_,1/m_+1)+4*pow(Se_,1/m_)+Se_);
            dkradS_ *= 1/(2*pow((1-Se_),0.5)*Se_);
            dkradS_ *=  1/(Smax_ - Smin_);
            dkradS_ *= kramax_;

            dkrbdS_ = 0.5 * (1-pow((1-pow(Se_,1/m_)),m_));
            dkrbdS_ *= ( 4 * pow(Se_,1/m_-0.5) * pow( (1-pow(Se_,1/m_)+VSMALL) , m_-1)) - ( pow((1-pow(Se_,1/m_)),m_) -1) / pow(Se_,0.5);
            dkrbdS_ *= 1/(Smax_ - Smin_);
            dkrbdS_ *= krbmax_;
        }
    void correctkra(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            kra_ = kramax_ * pow(1-Se_,0.5) * pow(1-pow(Se_,1/m_),2*m_);
            dkradS_ = - pow((1-pow(Se_,1/m_)+VSMALL),2*m_-1) * (-5*pow(Se_,1/m_+1)+4*pow(Se_,1/m_)+Se_);
            dkradS_ *= 1/(2*pow((1-Se_),0.5)*Se_);
            dkradS_ *=  1/(Smax_ - Smin_);
            dkradS_ *= kramax_;
        }
    void correctkrb(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            //surfaceScalarField sbinterp ("sbinterp", fvc::interpolate(Sb_));
            //krb_ = pow(sbinterp,lpower_) * pow(1-pow(1-pow(sbinterp,1/m_),m_),2);
            krb_ = krbmax_ * pow(Se_,lpower_) * pow(1-pow(1-pow(Se_,1/m_),m_),2);
            dkrbdS_ = lpower_ * (1-pow((1-pow(Se_,1/m_)),m_));
            dkrbdS_ *= ( 4 * pow(Se_,1/m_-0.5) * pow( (1-pow(Se_,1/m_)+VSMALL) , m_-1)) - ( pow((1-pow(Se_,1/m_)),m_) -1) * pow(Se_,lpower_-1);
            dkrbdS_ *= 1/(Smax_ - Smin_);
            dkrbdS_ *= krbmax_;

            // Se_= (Sb-Smin_)/(Smax_-Smin_);
            // krb_ = krbmax_ * pow(Se_,0.5) * pow(1-pow(1-pow(Se_,1/m_),m_),2);
            // dkrbdS_ = 0.5 * (1-pow((1-pow(Se_,1/m_)),m_));
            // dkrbdS_ *= ( 4 * pow(Se_,1/m_-0.5) * pow( (1-pow(Se_,1/m_)+VSMALL) , m_-1)) - ( pow((1-pow(Se_,1/m_)),m_) -1) / pow(Se_,0.5);
            // dkrbdS_ *= 1/(Smax_ - Smin_);
            // dkrbdS_ *= krbmax_;
        }

    void correctkrb(const surfaceScalarField& Sb, const volScalarField& Sb2)
        {
            //- face-based kr, evaluated at the face-interpolated water content Sb
            const surfaceScalarField Sminf(fvc::interpolate(Smin_));
            const surfaceScalarField Smaxf(fvc::interpolate(Smax_));
            const surfaceScalarField mf(fvc::interpolate(m_));
            const surfaceScalarField Sef
            (
                min(scalar(1) - SMALL, max(SMALL, (Sb - Sminf)/(Smaxf - Sminf)))
            );
            krbhave_ =
                fvc::interpolate(krbmax_)
              * pow(Sef, fvc::interpolate(lpower_))
              * pow(1 - pow(1 - pow(Sef,1/mf), mf), 2);

            Se_ = min(scalar(1) - SMALL, max(SMALL, (Sb2- Smin_)/(Smax_ - Smin_)));
            krb_ = krbmax_ * pow(Se_,lpower_) * pow(1-pow(1-pow(Se_,1/m_),m_),2);
            forAll(Se_, theta_i)
            {
                if (Se_[theta_i] >= 1)
                {
                    dkrbdS_[theta_i] = 0;
                }
                else
                {
                    dkrbdS_[theta_i] = lpower_[theta_i] * (1-pow((1-pow(Se_[theta_i],1/m_[theta_i])),m_[theta_i]));
                    dkrbdS_[theta_i] *= ( 4 * pow(Se_[theta_i],1/m_[theta_i]-0.5) * pow( (1-pow(Se_[theta_i],1/m_[theta_i])+VSMALL) , m_[theta_i]-1))
                                        - ( pow((1-pow(Se_[theta_i],1/m_[theta_i])),m_[theta_i]) -1) * pow(Se_[theta_i],lpower_[theta_i]-1);
                    dkrbdS_[theta_i] *= 1/(Smax_[theta_i] - Smin_[theta_i]);
                    dkrbdS_[theta_i] *= krbmax_[theta_i];
                }
            }
            // dkrbdS_ = lpower_ * (1-pow((1-pow(Se_,1/m_)),m_));
            // dkrbdS_ *= ( 4 * pow(Se_,1/m_-0.5) * pow( (1-pow(Se_,1/m_)+VSMALL) , m_-1)) - ( pow((1-pow(Se_,1/m_)),m_) -1) * pow(Se_,lpower_-1);
            // dkrbdS_ *= 1/(Smax_ - Smin_);
            // dkrbdS_ *= krbmax_;
            //dkrbdS_ *= krbmax_;
        }
    void correctkrb(const volScalarField& Sb, const label& celli)
        {
            scalar Se = (Sb[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            krb_[celli] = krbmax_[celli] * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_[celli]),m_[celli]),2);
        }

    // void correctkrb(const volScalarField& h)
    //     {
    //         FatalErrorIn("VanGenuchten.H") << "This is only for test cases in Celia" << abort(FatalError);
    //     }

    tmp<volScalarField> kr(const volScalarField& S)
        {
            volScalarField Se((S-Smin_)/(Smax_-Smin_));
            return krbmax_ * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_),m_),2);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
