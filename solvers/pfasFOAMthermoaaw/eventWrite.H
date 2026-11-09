//- write using event file
if (outputEventIsPresent)
{
    if (outputEvent.currentEventEndTime() <= runTime.value())
    {
        //- compute and write time interpolated fields
        volScalarField h_interpolated = outputEvent.timeInterpolate(h, runTime);
        volScalarField theta_interpolated = outputEvent.timeInterpolate(theta, runTime);
        volVectorField Utheta_interpolated = outputEvent.timeInterpolate(Utheta, runTime);
        surfaceScalarField phi_interpolated = outputEvent.timeInterpolate(phi, runTime);

        //- water CSV output
        if (CSVoutput)
        {
            waterMassBalanceCSV << outputEvent.currentEventEndTime() << " ";
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    waterMassBalanceCSV << " " << gSum(phi_interpolated.boundaryField()[patchi]);
                }
            }
            waterMassBalanceCSV << endl;
        }

        //- tracer outputs
        volScalarField C_interpolated = outputEvent.timeInterpolate(C_normal, runTime);
        scalar totalCmass_interpolated = outputEvent.timeInterpolate(totalCmassOld, totalCmass, runTime);

        //- CSV output at given times only
        if (CSVoutput)
        {
            CmassBalanceCSV << outputEvent.currentEventEndTime() << " " << totalCmass_interpolated;
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    CmassBalanceCSV << " " << gSum(phi_interpolated.boundaryField()[patchi]*C_interpolated.boundaryField()[patchi]);
                }
            }
            CmassBalanceCSV << endl;
        }

        outputEvent.updateIndex(runTime.value());
    }
}
else
{
    //- write fields using usual openfoam rules
    runTime.write();

    //- write CSV at all times
    if (CSVoutput)
    {
            waterMassBalanceCSV << runTime.timeName() << " ";
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    waterMassBalanceCSV << " " << gSum(phi.boundaryField()[patchi]);
                }
            }
            waterMassBalanceCSV << endl;

            CmassBalanceCSV << runTime.timeName() << " " << totalCmass;
            forAll(mesh.boundaryMesh(),patchi)
            {
                if (mesh.boundaryMesh()[patchi].type() == "patch")
                {
                    CmassBalanceCSV << " " << gSum(phi.boundaryField()[patchi]*C_normal.boundaryField()[patchi]);
                }
            }
            CmassBalanceCSV << endl;

    }
}
