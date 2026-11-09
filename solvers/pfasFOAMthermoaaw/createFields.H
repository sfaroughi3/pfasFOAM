//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

Info << "Reading Pressure head h \n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating deltah field (for Newton loop)" << endl;
volScalarField deltah
(
    IOobject
    (
        "deltah",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    h
);
deltah == dimensionedScalar("",dimLength,0);

//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info<< "Reading transportProperties" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//- list that receives event files of event-based boundary conditions
List<patchEventFile*> patchEventList;
eventInfiltration::setEventFileRegistry(&patchEventList, h.name());
eventFlux::setEventFileRegistry(&patchEventList, "C_aq");

//////////////////////////////////////////////////////////////////
////////////////////// Concentration Fields //////////////////////
//////////////////////////////////////////////////////////////////

const dimensionedScalar N_power ( transportProperties.lookup("Npower") );
Info << "Reading Concentration of Aqueous phase C_aq (C) \n" << endl;
volScalarField C_normal
(
    IOobject
    (
        "C_aq",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading Concentration of time based Solid phase C_s2 \n" << endl;
volScalarField C_sorption
(
    IOobject
    (
        "C_s2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading Concentration of time based air water interface C_aw2 \n" << endl;
volScalarField C_airwater
(
    IOobject
    (
        "C_aw2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Load all PFAS constants so they are easy to work with
const dimensionedScalar& xzero = transportProperties.lookup("xzero");
const dimensionedScalar& xone = transportProperties.lookup("xone");
const dimensionedScalar& xtwo = transportProperties.lookup("xtwo");
const dimensionedScalar& Fsorption = transportProperties.lookup("Fs");
const dimensionedScalar& Kfit = transportProperties.lookup("Kf");
const dimensionedScalar& rhob = transportProperties.lookup("rhob");
const dimensionedScalar& alpha_s = transportProperties.lookup("alpha_s");
const dimensionedScalar& alpha_aw = transportProperties.lookup("alpha_aw");
const dimensionedScalar& Faw = transportProperties.lookup("Faw");
const dimensionedScalar& a_PFAS = transportProperties.lookup("a");
const dimensionedScalar& b_PFAS = transportProperties.lookup("b");
const dimensionedScalar& sigma_0 = transportProperties.lookup("sigma_0");
const dimensionedScalar& T_amb = transportProperties.lookup("T_amb");
const dimensionedScalar& R_gas = transportProperties.lookup("R_gas");
const dimensionedScalar& thetamax = transportProperties.lookup("thetamax");

// Load Dispersion model
autoPtr<dispersionModel> dismodel = dispersionModel::New("DeffModel",transportProperties,mesh);
const volTensorField& dispersioncoeff = dismodel->Deff();
dimensionedScalar verytiny("small", dimensionSet(0,-3,0,0,1,0,0), 1e-15);
dimensionedScalar verytinyvalue("small", dimensionSet(0,0,0,0,0,0,0), 1e-15);

// Pre load several terms for calculation
C_normal.storePrevIter();
volScalarField Kaw = sigma_0 * b_PFAS / (R_gas * T_amb  * (a_PFAS + max(verytiny, C_normal) - verytiny));
volScalarField cpow = pow( max(verytiny, C_normal.prevIter() - verytiny), N_power-1);
volScalarField cpowcurr = pow(max(verytiny, C_normal) - verytiny, N_power);

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////

autoPtr<incompressiblePhase> phasetheta = incompressiblePhase::New(mesh,transportProperties,"theta");
volVectorField& Utheta = phasetheta->U();
const dimensionedScalar& rhotheta = phasetheta->rho();
const dimensionedScalar& mutheta = phasetheta->mu();
phasetheta->phi().writeOpt()=IOobject::NO_WRITE;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Utheta) & mesh.Sf()
);

/////////////////////////////////////////////////////////////////////////////
//////////////// RELATIVE PERMABILITY - CAPILLARY MODEL /////////////////////
/////////////////////////////////////////////////////////////////////////////

// Relative permeability / Capillary pressure models
autoPtr<twophasePorousMediumModel> porousModel = twophasePorousMediumModel::New("theta", mesh, transportProperties, phasetheta);
autoPtr<relativePermeabilityModel>& krModel = porousModel->krModel();
autoPtr<capillarityModel>& pcModel = porousModel->pcModel();

// Water content computation
Info<< "Computing field theta \n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimless,
    calculatedFvPatchScalarField::typeName
);
theta = pcModel->correctAndSb(h, C_normal);
theta.write();

#include "Aawtable.H"
/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Intrinsic permeability
porousModel->check_K();
const volScalarField& K = porousModel->K();

// permeability interpolation (harmonic by default)
surfaceScalarField Kf(fvc::interpolate(K,"K"));

// specific storage
dimensionedScalar Ss(transportProperties.lookupOrDefault<dimensionedScalar>("Ss",dimensionedScalar("Ss",dimless/dimLength,0.)));
Info << nl << "Reading speficic storage : Ss = " << Ss.value() << endl;

///////////////////////////////////////////////////////////////////
////////////////////////// FORCING TERMS //////////////////////////
///////////////////////////////////////////////////////////////////

//- For constant or event source injection/extraction
volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimless/dimTime,0)
);

////////////////////////////////////////////////////
//////////////////// OUTPUT CSV ////////////////////
////////////////////////////////////////////////////

bool CSVoutput=runTime.controlDict().lookupOrDefault<bool>("CSVoutput",true);
OFstream waterMassBalanceCSV("waterMassBalance.csv");
if (CSVoutput)
{
    waterMassBalanceCSV << "#Time ";
    forAll(mesh.boundaryMesh(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            waterMassBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
        }
    }
    waterMassBalanceCSV << endl;
}

OFstream CmassBalanceCSV("CmassBalance.csv");
if (CSVoutput)
{
    CmassBalanceCSV << "#Time TotalMass(kg)";
    forAll(mesh.boundaryMesh(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            CmassBalanceCSV << " flux(" << phi.boundaryField()[patchi].patch().name() << ")";
        }
    }
    CmassBalanceCSV << endl;
}

bool writeResiduals=false;
