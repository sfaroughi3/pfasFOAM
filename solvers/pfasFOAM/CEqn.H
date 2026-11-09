//- solve the PFAS transport equation
C_normal.storePrevIter();

Sw = theta/thetamax;
forAll(Sw, i)
{
    if(Sw[i] >= 1.0)
    {
        Sw[i] = 1.0;
    }
}
Aaw = Sw * Sw * xtwo + Sw * xone + xzero;
Kaw = sigma_0 * b_PFAS / (R_gas * T_amb  * (a_PFAS + max(verytiny, C_normal) - verytiny));
cpow = pow( max(verytiny, C_normal.prevIter() - verytiny), N_power-1);
cpowcurr = pow(max(verytiny, C_normal) - verytiny, N_power);

dismodel->correct(C_normal, Utheta, theta);

fvScalarMatrix C_sorpddt
    (
        fvm::ddt(C_sorption)
        + fvm::Sp(alpha_s, C_sorption)
        ==
        alpha_s * ((1 - Fsorption) * Kfit * cpowcurr)
    );

C_sorpddt.solve();

fvScalarMatrix C_airwaterddt
    (
        fvm::ddt(C_airwater)
        + fvm::Sp(alpha_aw, C_airwater)
        ==
        alpha_aw * ((1 - Faw) * Aaw * Kaw * C_normal)
    );
C_airwaterddt.solve();

fvScalarMatrix CEqn
    (
        fvm::ddt(theta,C_normal)
        + Faw * fvm::ddt(Aaw, Kaw, C_normal)
        //+ rhob * fvm::Sp(Fsorption * Kfit * cpow/ runTime.deltaT(), C_normal)
        + rhob * Fsorption*Kfit * fvm::ddt(cpow, C_normal) // shouldn't this be cpowcurr??
        + fvm::Sp(rhob * alpha_s * (1 - Fsorption) * Kfit * cpow, C_normal)
        + fvm::Sp(alpha_aw *(1 - Faw)*Aaw*Kaw, C_normal)
        + fvm::div(phi, C_normal, "div(phi,C_normal)")
        - fvm::laplacian(theta*dispersioncoeff, C_normal, "laplacian(Deff,C_normal)")
        ==
        rhob * alpha_s * C_sorption
        + alpha_aw  * C_airwater
        //+ rhob * Fsorption* Kfit * cpowold / runTime.deltaT()
    );

CEqn.relax();
CEqn.solve();

volScalarField deltaC = C_normal - C_normal.prevIter();

PFASResidual = gMax(mag(deltaC.internalField())()); /// gMax(mag(C_normal.primitiveField())());
Info << "Residual PFAS = " << PFASResidual << endl;

C_normal = max(C_normal, verytiny) - verytiny;
