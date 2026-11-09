{
    //- water mass balance terminal display
    Info << "Water mass balance (m3/s) : sourceTerm = " << fvc::domainIntegrate(sourceTerm).value() << " ; ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  gSum(phi.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;

    //- tracer mass balance terminal display
    Kaw = sigma_0 * b_PFAS / (R_gas * T_amb  * (a_PFAS + max(verytiny, C_normal.prevIter()) - verytiny));
    Info<< "Concentration: Min (C) = " << gMin(C_normal.internalField())
        << " Max(C_aq) = " << gMax(C_normal.internalField())
        << endl;
    scalar caqint = fvc::domainIntegrate(C_normal*theta).value();
    scalar cs2int = fvc::domainIntegrate(rhob*C_sorption).value();
    scalar csint = fvc::domainIntegrate(rhob * pow(C_normal,N_power) * Fsorption * Kfit).value() + cs2int;
    scalar caw2int = fvc::domainIntegrate(C_airwater).value();
    scalar cawint = fvc::domainIntegrate(C_normal * Aaw * Kaw * Faw).value() + caw2int;

    totalCmassOld = totalCmass;
    totalCmass = caqint + csint + cawint;

    Info<< " mass(C_aq) = " << caqint
        << " mass(C_s,2) = " << cs2int
        << " mass(C_s) = " << csint
        << " mass(C_aw,2) = " << caw2int
        << " mass(C_aw) = " << cawint
        << " mass(total) = " << totalCmass
        << " dCmax = " << dtManagerC.dVmax()*runTime.deltaTValue()
        << endl;

    Info << "C mass balance (kg/s) : ";
    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            scalar advectionterm = runTime.deltaTValue() * gSum(phi.boundaryField()[patchi]*C_normal.boundaryField()[patchi]);

            const scalarField theta_f = theta.boundaryField()[patchi];
            const tensorField D_f     = dispersioncoeff.boundaryField()[patchi];
            const scalarField snGradC = C_normal.boundaryField()[patchi].snGrad();
            const vectorField nf      = mesh.boundary()[patchi].nf()();
            const vectorField Sf      = mesh.Sf().boundaryField()[patchi];
            scalar diffusionterm = runTime.deltaTValue()
                * gSum(((theta_f * (D_f & nf)) & Sf) * snGradC);
            totalCinout += -advectionterm + diffusionterm;
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  -advectionterm + diffusionterm << " ; ";
        }
    }
    Info << "Total flux in and out over time = " << totalCinout << endl;
}
