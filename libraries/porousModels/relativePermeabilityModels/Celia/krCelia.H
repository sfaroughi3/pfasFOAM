/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Declaration of the Celia relative permeability model,
      kr = A / (A + |100 h|^gamma).

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
    Foam::relativePermeabilityModels::krCelia

SourceFiles
    krCelia.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef krCelia_H
#define krCelia_H

#include "relativePermeabilityModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
    Class krCelia Declaration
    \*---------------------------------------------------------------------------*/

class krCelia
:
    public relativePermeabilityModel
{
    //- Celia coefficients
    volScalarField gammapow_;
    volScalarField Aparam_;
    volScalarField alpha_;
    volScalarField beta_;

    //- cut-off point saturation (computed with entry pressure)

public:

    //- Runtime type //Information
    TypeName("rCelia");

    // Constructors

    //- Construct from components
    krCelia
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~krCelia()
        {}

    // Member Functions

    //- Correct the relative permeabilities
     void correct(const volScalarField& Sb)
        {
            FatalErrorIn("krCelia.H") << " Celia model can be used only with groundwater solvers" << abort(FatalError);
        }
    void correctkra(const volScalarField& Sb)
        {
            FatalErrorIn("krCelia.H") << " Celia model can be used only with groundwater solvers" << abort(FatalError);
        }
    void correctkrb(const volScalarField& Sb)
        {
            // Se_= (Sb-Smin_)/(Smax_-Smin_);
            volScalarField h = Sb;
            scalar Aparamval = Aparam_.internalField()[0];
            scalar gammapowval = gammapow_.internalField()[0];
            dkrbdS_.dimensions().reset(dimless);
            h.dimensions().reset(dimless);
            surfaceScalarField hface = fvc::interpolate(h);
            krb_ = Aparam_ / (Aparam_ + pow(mag(h*100), gammapow_));
            krbhave_ = Aparamval / (Aparamval + pow(mag(hface*100), gammapowval));

            //There is an issue with K(h), as the derivative with respect to theta is very expensive and annoying
            // As such, since this is really only for testing purposes, lets divide by Ch so that in the Jacobian
            // it cancels out. I know this is not great form, but to let everything be solved by one solver, this
            // is a sacrifice we have to make.

            volScalarField Setemp = alpha_ / (alpha_ + pow(mag(h) * 100, beta_));;
            volScalarField Chtemp = beta_ *  Setemp * pow(100, beta_) * (pow(mag(h), beta_-1))  * (Smax_ - Smin_) / (alpha_ + pow(mag(h) * 100, beta_));

            //Now actually calculate dkrbdS/Ch
            dkrbdS_ =  pow(100, gammapow_) * gammapow_ * pow(mag(h), gammapow_ - 1) * Aparam_ / pow(Aparam_ + pow(mag(h*100), gammapow_), 2);
            dkrbdS_ /= Chtemp;
            h.dimensions().reset(dimLength);


            //Info << "hi" << endl;
            // volScalarField hvalue = mag( pow(alpha_ * (1 / (Se_ + VSMALL) - 1), 1/beta_));
            // Info << hvalue.internalField()[10] << endl;
            // krb_ = Aparam_ / (Aparam_ + pow(mag(hvalue), gammapow_));
            // dkrbdS_ = Aparam_;
            // dkrbdS_ *= gammapow_ * alpha_ / beta_;
            // dkrbdS_ *= pow(hvalue, gammapow_ - 1.0) * mag( pow(alpha_ * (1 / (Se_ + VSMALL) - 1), 1/beta_ - 1)) ;
            // dkrbdS_ /= pow(Aparam_ + pow(hvalue, gammapow_), 2);
            // dkrbdS_ *= (Smax_-Smin_);
            // dkrbdS_ /= pow(Sb-Smin_ + VSMALL, 2);
            // Info << dkrbdS_.internalField()[10] << endl;
        }

    void correctkrb(const surfaceScalarField& Sb, const volScalarField& Sb2)
        {
            FatalErrorIn("krCelia.H") << "Celia model can be used only in test case with groundwatertransportPFAS" << abort(FatalError);
        }

    void correctkrb(const volScalarField& Sb, const label& celli)
        {
            FatalErrorIn("krCelia.H") << "Celia model can be used only in test case with groundwatertransportPFAS" << abort(FatalError);
        }
    // void correctkrb(const volScalarField& h)
    //     {
    //         krb_ = Ks_ * Aparam_ / (Aparam_ + pow(pos(h), gammapow_));
    //     }
    tmp<volScalarField> kr(const volScalarField& S)
        {
            FatalErrorIn("krCelia.H") << "Celia model can be used only in test case with groundwatertransportPFAS" << abort(FatalError);
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
