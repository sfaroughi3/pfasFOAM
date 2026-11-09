const scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

const scalar minDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 0);

const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

const scalar truncationError_global =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError",0.01);

const scalar truncationError_h =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError_h",truncationError_global);

const scalar dTFactDecrease =
    runTime.controlDict().lookupOrDefault<scalar>("dTFactDecrease",0.8);

//- time step reduction factor applied when a patch event starts
const scalar dTFactEvent =
    runTime.controlDict().lookupOrDefault<scalar>("dTFactEvent",0.1);

timestepManager dtManager(runTime, theta, truncationError_h);

const scalar truncationError_C =
    runTime.controlDict().lookupOrDefault<scalar>("truncationError_C",truncationError_global);

timestepManager dtManagerC(runTime, C_normal, truncationError_C);

//- Derivative initialisation to keep 1st user-defined time step
if (adjustTimeStep)
{
    Info << nl << "Timestepping for h equation is based on time-scheme truncation error with :"
        << nl << "{"
        << nl << "    truncationError = " << truncationError_h
        << nl << "    dTFactDecrease  = " << dTFactDecrease
        << nl << "}"
        << endl;

    Info << nl << "Timestepping for C equation is based on time-scheme truncation error with :"
        << nl << "{"
        << nl << "    truncationError = " << truncationError_C
        << nl << "}"
        << endl;
}

const bool eventTimeTracking =
     runTime.controlDict().lookupOrDefault("eventTimeTracking", false);

// ************************************************************************* //
