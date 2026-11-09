// ---- Build once (e.g. at the end of createFields.H) -----------------------
IOdictionary awiDict
(
    IOobject("awiProperties", runTime.constant(), mesh,
             IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE)
);

const scalar porosity = awiDict.lookup<scalar>("porosity");  // "n" in Eq. 3
const scalar sigmaAW  = awiDict.lookup<scalar>("sigma");     // [N/m]
const scalar alphaVG  = awiDict.lookup<scalar>("alphaVG");   // [1/m]
const scalar nVG      = awiDict.lookup<scalar>("nVG");       // [-]
const scalar theta_s  = awiDict.lookupOrDefault<scalar>("theta_s", 0.0);
const scalar theta_r  = awiDict.lookupOrDefault<scalar>("theta_r", 0.0);
const scalar rhoW     = awiDict.lookupOrDefault<scalar>("rhoW", 1000.0);
const scalar gMag     = awiDict.lookupOrDefault<scalar>("gMag", 9.81);

const scalar mVG     = 1.0 - 1.0/nVG;
const label  nAwTab  = 20000;
const scalar SeMin   = awiDict.lookupOrDefault<scalar>("SeMin", 1e-4);
const scalar dSeTab  = (1.0 - SeMin)/nAwTab;

scalarField hcIntTab(nAwTab + 1, 0.0);   // integral from Se_k to 1 of h_c dSe [m]
for (label k = nAwTab - 1; k >= 0; --k)
{
    const scalar SeMid = SeMin + (k + 0.5)*dSeTab;
    const scalar hc =
        Foam::pow(Foam::pow(SeMid, -1.0/mVG) - 1.0, 1.0/nVG)/alphaVG;
    hcIntTab[k] = hcIntTab[k + 1] + hc*dSeTab;
}

volScalarField Aaw
(
    IOobject("Aaw", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimless/dimLength, 0)
);

#include "Updateaw.H"
Aaw.oldTime();
