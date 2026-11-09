const dictionary& PicardDict = mesh.solutionDict().subOrEmptyDict("Picard");
const scalar tolerancePicard = PicardDict.lookupOrDefault<scalar>("tolerance", GREAT);
label maxIterPicard = PicardDict.lookupOrDefault<label>("maxIter",10);
const dictionary& NewtonDict = mesh.solutionDict().subOrEmptyDict("Newton");
const scalar toleranceNewton = NewtonDict.lookupOrDefault<scalar>("tolerance", GREAT);
const label maxIterNewton = NewtonDict.lookupOrDefault<label>("maxIter",10);
const dictionary& PFASCDict = mesh.solutionDict().subOrEmptyDict("PFAS");
const scalar tolerancePFAS = PFASCDict.lookupOrDefault<scalar>("tolerance", GREAT);
const label maxIterPFAS = PFASCDict.lookupOrDefault<label>("maxIter",10);

Info << nl << "Picard loop control" << nl << "{"
<< nl << "    tolerance = " << tolerancePicard
<< "    maximum number of iteration = " << maxIterPicard
<< nl << "}" << endl;
Info << nl << "Newton loop control" << nl << "{"
<< nl << "    tolerance = " << toleranceNewton
<< "    maximum number of iteration = " << maxIterNewton
<< nl << "}" << endl;
Info << nl << "PFAS loop control" << nl << "{"
<< nl << "    tolerance = " << tolerancePFAS
<< "    maximum number of iteration = " << maxIterPFAS
<< nl << "}" << endl;

if (!PicardDict.found("tolerance") && !NewtonDict.found("tolerance"))
{
    FatalErrorIn("readPicardNewtonControls.H")
        << "Missing keyword tolerance in system/fvSolution.Picard or system/fvSolution.Newton"
            << abort(FatalError);
}

//- to switch between conservative mixed formulation (by default) and head pressure (non conservative)
bool massConservative = transportProperties.lookupOrDefault<bool>("massConservative",true);
bool averagekf = transportProperties.lookupOrDefault<bool>("averageKf",true);
bool hbasedkf = transportProperties.lookupOrDefault<bool>("hbasedkf",false);

Info << "Mass conervative " << massConservative << nl
<<  "Solve Kf then average " << averagekf << nl
<< "An H based relative permeability " << hbasedkf << endl;
