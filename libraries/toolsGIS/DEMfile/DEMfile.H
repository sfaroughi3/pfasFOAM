/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / toolsGIS

  Redistributed as part of pfasFOAM by:
      Energy & Intelligence Lab (Eilab) - University of Utah
      (Nicholas T. Baker & Salah A. Faroughi), 2026

  Modifications (2026):
      - None to the code of this file; it is unchanged from the upstream
        source (apart from this header).

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
    Foam::DEMfile

Description
    DEM file handler/reader with bilinear interpolation for the
    porousMultiphaseFoam toolbox

SourceFiles
    DEMfile.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef DEMfile_H
#define DEMfile_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class fluidDEMfile Declaration
    \*---------------------------------------------------------------------------*/

class DEMfile
{

private:

    word name_;

    scalar x0_,y0_;
    scalar dx_,dy_;
    label nx_,ny_;
    scalarList zvalues_;

public:

    // Constructors

    //- Copy constructor
    DEMfile
    (
        const DEMfile& fileToCopy
    );

    //- Construct from Istream
    DEMfile
    (
        const word& fileName
    );
    
    // Selectors

    //- Destructor
    virtual ~DEMfile();

    // Member Functions
    
    //- Access function
    const word& name() const {return name_;}
    const scalar& x0() const {return x0_;}
    const scalar& y0() const {return y0_;}
    const scalar& dx() const {return dx_;}
    const scalar& dy() const {return dy_;}
    const scalarList& z() const {return zvalues_;}

    //- Bilinear interpolation using given point
    scalar interpolate(const point& position);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
