Info << "Reading field p_rgh\n"
     << endl;
volScalarField p_rgh(
    IOobject(
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

Info << "Reading field U\n"
     << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

Info << "Reading field variance\n"
     << endl;
volScalarField var(
    IOobject(
        "var",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

#include "createPhi.H"

Info << "Reading transportProperties\n"
     << endl;
incompressibleTwoPhaseMixture mixture(U, phi);

volScalarField &alpha1(mixture.alpha1());
volScalarField &alpha2(mixture.alpha2());

const dimensionedScalar &rho1 = mixture.rho1();
const dimensionedScalar &rho2 = mixture.rho2();

dimensionedScalar Dab("Dab", dimViscosity, mixture);

// Read the reciprocal of the turbulent Schmidt number
dimensionedScalar alphatab("alphatab", dimless, mixture);

// Need to store rho for ddt(rho, U)
volScalarField rho("rho", alpha1 *rho1 + alpha2 * rho2);
rho.oldTime();

// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi(
    IOobject(
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    rho1 *phi);

// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence(
    incompressible::momentumTransportModel::New(U, phi, mixture));

// Accessing the delta filter
volScalarField LESdelta(
    IOobject(
        "LESdelta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("LESdelta", dimLength, SMALL));

if (mesh.objectRegistry::foundObject<volScalarField>("delta"))
{
    const volScalarField &delta =
        mesh.objectRegistry::lookupObject<volScalarField>("delta");
    
    LESdelta = delta;
}

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p(
    IOobject(
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    p_rgh + rho * gh);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue);

if (p_rgh.needReference())
{
    p += dimensionedScalar(
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell));
    p_rgh = p - rho * gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

Info << "Reading chemicalReactionDictionary" << endl
     << endl;

Cphi Cphi_(mesh, turbulence);
chemReac betaPDF(mesh, alpha1, var, phi, turbulence, alphatab, Dab);

// Defining neutralization
betaPDF.lookUpTableB();