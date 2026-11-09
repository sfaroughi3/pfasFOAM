{
    const scalarField& thetaI  = theta.primitiveField();
    scalarField&       AawI = Aaw.primitiveFieldRef();

    forAll(AawI, celli)
    {
        scalar Se = (thetaI[celli] - theta_r)/(theta_s - theta_r);
        Se = min(max(Se, SeMin), 1.0);

        const scalar x = (Se - SeMin)/dSeTab;
        const label  k = min(label(x), nAwTab - 1);
        const scalar w = x - k;
        const scalar hcInt = (1.0 - w)*hcIntTab[k] + w*hcIntTab[k + 1];

        AawI[celli] = porosity*rhoW*gMag*(theta_s - theta_r)*hcInt/sigmaAW;   // [1/m]
    }
    Aaw.correctBoundaryConditions();
}
