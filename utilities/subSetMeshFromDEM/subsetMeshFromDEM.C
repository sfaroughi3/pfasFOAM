/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / utilities/subSetMeshFromDEM

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

Application
    subsetMeshFromDEM

Description
    Generate automatically a subsetMesh using two DEM files (.xyz)
    Boundary faces created from first file are added to the patch "top"
    Boundary faces created from second file are added to the patch "bottom"

Usage
    subsetMeshFromDEM topDEMfile bottomDEMfile

  Previous license:
      Copyright (C) Toulouse INP, Pierre Horgue
      This file is part of porousMultiphaseFoam, an extension of OpenFOAM
      developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
      flows through porous media.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "DEMfile.H"
#include "cellSet.H"
#include "fvMeshSubset.H"
#include "fvMeshTools.H"

int main(int argc, char *argv[])
{

    argList::validArgs.append("DEMtop"),
    argList::validArgs.append("DEMbottom");
    
    Foam::argList args(argc,argv); 

    #include "createTime.H"
    #include "createMesh.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- select using DEM-top file
    Info << nl << "Read DEM-top file";
    DEMfile topFile(args[1]);

    cellSet topCells(mesh, "topCells",  IOobject::NO_READ);
    forAll(mesh.C(),celli)
    {
        scalar ztop = topFile.interpolate(mesh.C()[celli]);
        if (mesh.C()[celli].z() < ztop) topCells.insert(celli);
    } 

    fvMeshSubset topMesh(mesh);
    topMesh.setCellSubset(topCells, -1, true);
    const fvMesh& interMesh = topMesh.subMesh();
    
    //- select using DEM-bottom file
    Info << nl << "Read DEM-bottom file";
    DEMfile bottomFile(args[2]);

    cellSet bottomCells(interMesh, "bottomCells",  IOobject::NO_READ);
    forAll(interMesh.C(),celli)
    {
        scalar zbottom = bottomFile.interpolate(interMesh.C()[celli]);
        if (interMesh.C()[celli].z() > zbottom) bottomCells.insert(celli);
    } 

    fvMeshSubset finalMesh(topMesh.subMesh());
    finalMesh.setCellSubset(bottomCells, -1, true);

    //- write final mesh
    fvMesh& outMesh = finalMesh.subMesh();

    outMesh.setInstance(runTime.constant());
    outMesh.write();
        
    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
