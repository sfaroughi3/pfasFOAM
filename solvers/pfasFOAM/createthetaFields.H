//- relative permeability (kr)
surfaceScalarField thetasurf = fvc::interpolate(theta);

if (!hbasedkf)
{
    krModel->correctkrb(thetasurf, theta);

}else
{
    krModel->correctkrb(h);
}

const volScalarField& krtheta = krModel->krb();
surfaceScalarField krthetaf ("krthetaf",fvc::interpolate(krtheta,"krtheta"));

if (averagekf)
{
        krthetaf = surfaceScalarField ("krthetaf",fvc::interpolate(krtheta,"krtheta"));
} else
{
        krthetaf = krModel->krbave();
}
#include "inletKr.H"

//- mobility and fractional flow
surfaceScalarField Lf ("Lf",rhotheta*Kf*krthetaf/mutheta);
surfaceScalarField Mf ("Mf",mag(g)*Lf);

//- fluxes depending on saturation
surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());
Utheta.correctBoundaryConditions();
forAll(mesh.boundary(),patchi)
{
    if (isA< fixedValueFvPatchField<vector> >(Utheta.boundaryField()[patchi]))
    {
        phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    }
}

//- Test if gravity is present
if (mag(g).value() == 0)
{
    FatalErrorIn("createthetaFields.H")
        << " Magnitude of gravity mag(g) equal to zero " << abort(FatalError);
}

//- null Pc field for darcy velocity boundary conditions
surfaceScalarField phiPcNull("phiPc",0*phiG);

//- h equation residual
volScalarField ResiduN
    (
        Ss*pcModel->Se() * fvc::ddt(h)
           + fvc::ddt(theta)
           - fvc::laplacian(Mf,h)
           + fvc::div(phiG)
           + sourceTerm
    );
if (writeResiduals)
{
    ResiduN.writeOpt()=IOobject::AUTO_WRITE;
    ResiduN.rename("hEqnResidual");
}

// Integrate the inital fields and values for the mass balance statistics.
scalar caqint = fvc::domainIntegrate(C_normal*theta).value();
scalar csint = fvc::domainIntegrate(rhob * pow(C_normal,N_power) * Fsorption * Kfit).value() + fvc::domainIntegrate(rhob*C_sorption).value();
scalar cawint = fvc::domainIntegrate(C_normal * Aaw * Kaw * Faw).value() + fvc::domainIntegrate(C_airwater).value();
scalar totalCinout = caqint + csint + cawint;

//- total PFAS mass (aqueous + solid + air-water) at the current and previous time step, for the CSV output
scalar totalCmass = caqint + csint + cawint;
scalar totalCmassOld = totalCmass;
Info<< " mass(C_aq) = " << caqint
    << " mass(C_s) = " << csint
    << " mass(C_aw) = " << cawint
    << " mass(total) = " << caqint + csint + cawint
    << endl;
