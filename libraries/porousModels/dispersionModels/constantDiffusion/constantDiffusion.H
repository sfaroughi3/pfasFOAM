/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / porousModels

  Copyright (C) 2026 Energy & Intelligence Lab (Eilab) - University of Utah

  Authors: Nicholas T. Baker & Salah A. Faroughi
  Affiliation: Eilab - University of Utah

  Description:
      Declaration of a dispersion model with constant diffusion coefficients
      Dx, Dy and Dz read from the constantDiffusionCoeffs dictionary.

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
    Foam::dispersionModels::constantDiffusion

SourceFiles
    constantDiffusion.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef constantDiffusion_H
#define constantDiffusion_H

#include "dispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{

/*---------------------------------------------------------------------------*\
                  Class constantDiffusion Declaration
\*---------------------------------------------------------------------------*/

class constantDiffusion
:
    public dispersionModel
{
    //- Private data

    //- Dictionary for model coefficients
    dictionary constantDiffusionCoeffs_;

    //- model coefficients
    volScalarField Dx_,Dy_,Dz_;

public:

    //- Runtime type information
    TypeName("constantDiffusion");

    // Constructors

    //- Construct from components
    constantDiffusion
    (
        const word& name,
        const dictionary& transportProperties,
        const fvMesh& mesh
    );

    //- Destructor
    ~constantDiffusion()
        {}

    // Member Functions

    //- Correct the dispersion giving water content or porosity field
    void correct(const volScalarField& C, const volVectorField& U, const volScalarField& theta)
        {
            dimensionedScalar USMALL("",dimLength/dimTime,1e-36);

            forAll(Deff_.internalField(),patchi)
            {
                Deff_.ref()[patchi] = tensor(Dx_.internalField()[0], 0, 0, 0, Dy_.internalField()[0], 0, 0, 0, Dz_.internalField()[0]);
            }
            forAll(Deff_.boundaryField(), patchi)
            {
                Deff_.boundaryFieldRef()[patchi] = tensor(Dx_.internalField()[0], 0, 0, 0, Dy_.internalField()[0], 0, 0, 0, Dz_.internalField()[0]);
            }
        }

    //- Correct the dispersion giving scalar porosity
    void correct(const volScalarField& C, const volVectorField& U, const dimensionedScalar& eps)
        {
            dimensionedScalar USMALL("",dimLength/dimTime,1e-36);

            forAll(Deff_.internalField(),patchi)
            {
                Deff_.ref()[patchi] = tensor(Dx_.internalField()[0], 0, 0, 0, Dy_.internalField()[0], 0, 0, 0, Dz_.internalField()[0]);
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dispersionModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
