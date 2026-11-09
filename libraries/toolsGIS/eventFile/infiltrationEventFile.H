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
    Foam::infiltrationEventFile

Description
    Event file handler/reader for the porousMultiphaseFoam toolbox which
    contain 'M' dates with uniform or non uniform infiltration values where
    infiltration is expressed in m/s. The class reads the uniform file :

    date time1
    infiltration1
    date time2
    infiltration2
    date timeM
    infiltrationM

    or the non-uniform file :

    date time1
    infiltrationCell11 infiltrationCell21 ... infiltrationCellN1
    date time2
    infiltrationCell12 infiltrationCell22 ... infiltrationCellN2
    date timeM
    infiltrationCell1M infiltrationCell2M ... infiltrationCellNM

SourceFiles
    infiltrationEventFile.C

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#ifndef infiltrationEventFile_H
#define infiltrationEventFile_H

#include "fvMesh.H"
#include "eventFile.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class fluidinfiltrationEventFile Declaration
    \*---------------------------------------------------------------------------*/

class infiltrationEventFile
:
    public eventFile
{
    //- True if infiltration is uniform
    bool uniform_;

public:

    // Constructors    
    //- Copy constructor
    infiltrationEventFile
    (
        const infiltrationEventFile& fileToCopy
    );

    //- Construct from Istream
    infiltrationEventFile
    (
        const word& fileName
    );
    
    // Selectors

    //- Destructor
    ~infiltrationEventFile();

    // Member Functions
    
    //- Access function
    const bool& isUniform() {return uniform_;}

    //- reset infiltration as uniform if size does not correspond to mesh size
    void setInfiltrationAsUniform() {uniform_ = true;}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
