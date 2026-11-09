//- update patch event values before solving h equation so that
//  event-based boundary conditions (eventInfiltration/eventFlux) are
//  synchronized with the current time step
forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);

//- compute water source term
if (eventIsPresent_water)
{
    event_water.updateValue(runTime);
    sourceTerm = event_water.dtValuesAsField();
}
