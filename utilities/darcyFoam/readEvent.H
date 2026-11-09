/*---------------------------------------------------------------------------*\
  Project: pfasFOAM / utilities/darcyFoam

  Redistributed as part of pfasFOAM by:
      Energy & Intelligence Lab (Eilab) - University of Utah
      (Nicholas T. Baker & Salah A. Faroughi), 2026

  Modifications (2026):
      - None to the code of this file; it is unchanged from the upstream
        source (apart from this header).

  License:
      SPDX-License-Identifier: GPL-3.0-or-later

      This file is part of pfasFOAM.

  Previous license:
      The upstream porousMultiphaseFoam file carries no header of its own.
      It is part of porousMultiphaseFoam, an extension of OpenFOAM
      Copyright (C) Toulouse INP, Pierre Horgue (as stated in the headers of
      the other porousMultiphaseFoam source files), distributed under the GNU
      General Public License, version 3 or (at your option) any later version
      (see COPYING_OPENFOAM).
\*---------------------------------------------------------------------------*/

bool sourceEventIsPresent = transportProperties.found("eventFileWaterSource");
word sourceEventFileName = transportProperties.lookupOrDefault<word>("eventFileWaterSource","");
sourceEventFile sourceEvent(sourceEventFileName);
if (sourceEventIsPresent)
{
    sourceEvent.updateIndex(runTime.startTime().value());
    sourceEvent.updateValue(runTime);
    sourceEvent.storeOldValues();
}
