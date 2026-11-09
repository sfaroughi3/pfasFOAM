{
    //- Set h equal to old-time value
    h = h.oldTime();
    porousModel->rewindTime();
    C_normal = C_normal.oldTime();
    C_airwater = C_airwater.oldTime();
    C_sorption = C_sorption.oldTime();
    //- Rewind time
    runTime.setTime(runTime.value()-runTime.deltaTValue(),runTime.timeIndex());
    //- Recompute time step
    iterPFAS += 1;
    iterPicard += 1;
    iterNewton += 1;
    #include "setDeltaT.H"
    //- Update new time
    runTime.setTime(runTime.value()+runTime.deltaTValue(),runTime.timeIndex());
    //- Update properties
    #include "updateProperties.H"
}
