//- update saturation
theta = pcModel->correctAndSb(h, C_normal);

//- relative permeability computation
thetasurf = fvc::interpolate(theta);

// If correlation in Celia is desired or exponential models, need h for input
if (!hbasedkf)
{
    krModel->correctkrb(thetasurf, theta);
}else
{
    krModel->correctkrb(h);
}

if (averagekf)
{
        krthetaf = surfaceScalarField ("krthetaf",fvc::interpolate(krtheta,"krtheta"));
} else
{
        krthetaf = krModel->krbave();
}
#include "inletKr.H"
//- mobility computation
Lf = rhotheta*Kf*krthetaf/mutheta;
Mf = mag(g)*Lf;

//- compute fluxes
phiG = (Lf * g) & mesh.Sf();
phi = phiG-(Mf*fvc::snGrad(h))*mesh.magSf();
Utheta = fvc::reconstruct(phi);
Utheta.correctBoundaryConditions();
forAll(mesh.boundary(),patchi)
{
    if (isA< fixedValueFvPatchField<vector> >(Utheta.boundaryField()[patchi]))
    {
        phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    }
}
h.correctBoundaryConditions();

//- refresh saturation so boundary theta is consistent with the updated h_b
theta = pcModel->correctAndSb(h, C_normal);

//- Compute h equation residual
ResiduN = - fvc::laplacian(Mf,h) + fvc::div(phiG) + porousModel->exchangeTerm() + sourceTerm;
if (!steady)
{
    ResiduN += Ss*pcModel->Se() * fvc::ddt(h);
    if (massConservative) ResiduN += fvc::ddt(theta);
    else ResiduN += pcModel->Ch() * fvc::ddt(h);
}

forAll(fixedPotentialIDList, pointi) ResiduN[fixedPotentialIDList[pointi]]=0;
hEqnRes = gMax(mag(ResiduN.internalField())());
hEqnResidual = gMax(mag(deltah.internalField())());
