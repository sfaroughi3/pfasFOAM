if(adjustTimeStep)
{
    //- in case of non convergence of the Picard or Newton  algorithm
    if (iterPicard == maxIterPicard+1 || iterNewton == maxIterNewton+1 || iterPFAS == maxIterPFAS+1)
    {
        runTime.setDeltaT(max(minDeltaT, dTFactDecrease*runTime.deltaTValue()));
    }
    else
    {
        //- h-variation control
        scalar dtForh = dtManager.computeTimestep();

        //- C-variation control
        scalar dtForC = dtManagerC.computeTimestep();
        //- set minimal delta t
        runTime.setDeltaT
            (
                max(
                    minDeltaT,
                    min(
                        min(dtForh, dtForC),
                        min(1.2*runTime.deltaTValue(), maxDeltaT)
                    )
                )
            );
        if (min(dtForh,dtForC ) == dtForh)
        {
            Info << "h is controlling time step" << endl;
        } else
        {
            Info << "C is controlling time step" << endl;
        }
    }

    //- Patch event starting (stiff boundary condition change)
    //  => strong time step reduction
    forAll(patchEventList,patchEventi)
    {
        if (patchEventList[patchEventi]->iterator() != patchEventIteratorOld[patchEventi])
        {
            patchEventIteratorOld[patchEventi] = patchEventList[patchEventi]->iterator();
            Info << "Patch event " << patchEventList[patchEventi]->name()
                << " starting => reducing time step by factor dTFactEvent" << endl;
            runTime.setDeltaT(max(minDeltaT, dTFactEvent*runTime.deltaTValue()));
        }
    }

    //-Adjust time step to explicitly compute (water) event time
    scalar timeOfNextEvent = GREAT;
    if (eventTimeTracking)
    {
        if (eventIsPresent_water) timeOfNextEvent = min(timeOfNextEvent,event_water.currentEventEndTime());

        scalar timeToNextEvent = timeOfNextEvent-runTime.value();
        scalar nSteps =  timeToNextEvent/runTime.deltaTValue();
        if ((nSteps < labelMax) && (nSteps != 0))
        {
            const label nStepsToNextEvent = label(max(nSteps, 1) + 0.99);
            runTime.setDeltaTNoAdjust(timeToNextEvent/nStepsToNextEvent);
        }

        //- To handle close event times (inferior to current timestep)
        if (nSteps == 0)
        {
            scalar timeToCloseEvent = GREAT;
            if (eventIsPresent_water)
            {
                if (event_water.currentEventEndTime() != runTime.value())
                {
                    timeToCloseEvent = min(timeToCloseEvent,event_water.currentEventEndTime()-runTime.value());
                }
            }
            forAll(patchEventList,patchEventi)
            {
                if (patchEventList[patchEventi]->currentEventEndTime() != runTime.value())
                {
                    timeToCloseEvent = min(timeToCloseEvent,patchEventList[patchEventi]->currentEventEndTime()-runTime.value());
                }
            }
            runTime.setDeltaTNoAdjust(min(runTime.deltaTValue(),timeToCloseEvent));
        }
    }

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;

}

// ************************************************************************* //
