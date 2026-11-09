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
    Foam::relativePermeabilityModels::krBrooksAndCorey

Description
    Standard Brooks and Corey relative permeability model.

SourceFiles
    krBrooksAndCorey.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef krBrooksAndCorey_H
#define krBrooksAndCorey_H

#include "relativePermeabilityModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
    Class krBrooksAndCorey Declaration
\*---------------------------------------------------------------------------*/

class krBrooksAndCorey
:
    public relativePermeabilityModel
{
    volScalarField n_;

    //- end points
    volScalarField kramax_;
    volScalarField krbmax_;

public:

    //- Runtime type information
    TypeName("BrooksAndCorey");

    // Constructors

    //- Construct from components
    krBrooksAndCorey
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~krBrooksAndCorey()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            kra_ = kramax_* pow((scalar(1)-Se_),n_);
            krb_ = krbmax_ * pow(Se_,n_);
            dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
            dkrbdS_ = krbmax_*n_*pow(Se_,n_-1)/(Smax_- Smin_);
        }
    void correctkra(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            kra_ = kramax_* pow((scalar(1)-Se_),n_);
            dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
        }
    void correctkrb(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            krb_ = krbmax_ * pow(Se_,n_);
            dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
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
            krbhave_ = fvc::interpolate(krbmax_) * pow(Sef, fvc::interpolate(n_));

            Se_ = min(scalar(1) - SMALL, max(SMALL, (Sb2 - Smin_)/(Smax_ - Smin_)));
            krb_ = krbmax_ * pow(Se_,n_);
            dkrbdS_ = krbmax_*n_*pow(Se_,n_-1)/(Smax_ - Smin_);
        }
    void correctkrb(const volScalarField& Sb, const label& celli)
        {
            scalar Se = (Sb[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            krb_[celli] = krbmax_[celli] * pow(Se,n_[celli]);
        }

    // void correctkrb(const volScalarField& h)
    //     {
    //         FatalErrorIn("BrooksAndCorey.H") << "This is only for test cases in Celia" << abort(FatalError);
    //     }

    tmp<volScalarField> kr(const volScalarField& S)
        {
            volScalarField Se((S-Smin_)/(Smax_-Smin_));
            return krbmax_ * pow(Se,n_);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
