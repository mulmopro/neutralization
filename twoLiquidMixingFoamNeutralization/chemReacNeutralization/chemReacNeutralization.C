/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::chemReac

Description

SourceFiles
    chemReac.C

\*---------------------------------------------------------------------------*/

#include "chemReacNeutralization.H"
#include "cellSet.H"

#include "CEqn.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chemReac::chemReac(
    const fvMesh &mesh,
    const volScalarField &alpha1,
    const volScalarField &var,
    const surfaceScalarField &phi,
    const incompressible::momentumTransportModel &turbulence,
    const dimensionedScalar &alphatab,
    const dimensionedScalar &Dab)

    : IOdictionary(

          IOobject(
              "chemicalReactionProperties",
              mesh.time().constant(),
              mesh,
              IOobject::MUST_READ_IF_MODIFIED,
              IOobject::NO_WRITE)),

      mesh_(mesh),

      alpha1_(alpha1),

      var_(var),

      phi_(phi),

      turbulence_(turbulence),

      alphatab_(alphatab),

      D_(Dab),

      HCl_(readScalar(subDict("Concentrations").lookup("A"))),

      NaOH_(readScalar(subDict("Concentrations").lookup("B"))),

      rhow_(readScalar(lookup("waterDensity"))),

      r_(readScalar(lookup("stoichiometricRatio"))),

      zeta_s(r_ * HCl_ / (r_ * HCl_ + NaOH_)),

      OH(
          IOobject(
              "OH",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh)
{
    Foam::chemReac::mixtureFraction();
    Foam::chemReac::segregationMatrix();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::chemReac::~chemReac()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

double Foam::chemReac::incbeta(double a, double b, double x)
{
    if (x < 0.0 || x > 1.0)
        return 1.0 / 0.0;

    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    if (x > (a + 1.0) / (a + b + 2.0))
    {
        return (1.0 - incbeta(b, a, 1.0 - x)); /*Use the fact that beta is symmetrical.*/
    }

    /*Find the first part before the continued fraction.*/
    const double lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
    const double front = exp(log(x) * a + log(1.0 - x) * b - lbeta_ab) / a;

    /*Use Lentz's algorithm to evaluate the continued fraction.*/
    double f = 1.0, c = 1.0, d = 0.0;

    int i, m;
    for (i = 0; i <= 200; ++i)
    {
        m = i / 2;

        double numerator;
        if (i == 0)
        {
            numerator = 1.0; /*First numerator is 1.0.*/
        }
        else if (i % 2 == 0)
        {
            numerator = (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m)); /*Even term.*/
        }
        else
        {
            numerator = -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1)); /*Odd term.*/
        }

        /*Do an iteration of Lentz's algorithm.*/
        d = 1.0 + numerator * d;
        if (fabs(d) < 1e-30)
            d = 1e-30;
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if (fabs(c) < 1e-30)
            c = 1e-30;

        const double cd = c * d;
        f *= cd;

        /*Check for stop.*/
        if (fabs(1.0 - cd) < 1.0e-8)
        {
            return front * (f - 1.0);
        }
    }

    return 1.0 / 0.0; /*Needed more loops, did not converge.*/
}

void Foam::chemReac::mixtureFraction()
{
    std::cout.precision(4);
    Info << "mixture fraction method" << endl
         << "-----------------------" << endl;

    double a, b;
    a = 0.0;
    b = 1.0;
    double h;
    h = (b - a) / (N_ - 1);

    for (int i = 0; i <= N_ - 1; i++)
    {
        mmf[i] = a + i * h; // Mixture fraction array
        Info << mmf[i] << " ";
    }
    Info << endl
         << endl;
}

void Foam::chemReac::segregationMatrix()
{
    Info << "Segregation Matrix method" << endl
         << "-------------------------" << endl;

    for (int i = 0; i < N_; i++)
    {
        for (int j = 0; j < N_; j++)
        {
            Is[i][j] = mmf[i];
            Info << Is[i][j] << " ";
        }
        Info << endl;
    }
    Info << endl
         << endl;
}

void Foam::chemReac::lookUpTableB()
{
    for (int j = 1; j < N_ - 1; j++)

    {
        double mean_mf = mmf[j]; // mean mixture fraction
        for (int i = 1; i < N_ - 1; i++)
        {
            double I_s = Is[i][j]; // segregation coefficient
            double ni = mean_mf * (1 / I_s - 1);
            double w = (1 - mean_mf) * (1 / I_s - 1);
            double beta = Foam::tgamma(ni) * Foam::tgamma(w) / Foam::tgamma(ni + w);
            double beta1 = Foam::tgamma(ni + 1) * Foam::tgamma(w) / Foam::tgamma(ni + w + 1);
            // Analytical solution
            OH_[i][j] = NaOH_ / (beta * (1 - zeta_s)) * (beta1 * (1 - Foam::chemReac::incbeta(ni + 1, w, zeta_s)) - zeta_s * beta * (1 - Foam::chemReac::incbeta(ni, w, zeta_s)));
        }
    }

    // implementing known analytical solutions
    for (int j = 0; j <= N_ - 1; j++)
    {
        if (mmf[j] >= 0 && mmf[j] <= zeta_s)
        {
            OH_[0][j] = 0;
        }
        else if (mmf[j] >= zeta_s && mmf[j] <= 1)
        {
            OH_[0][j] = NaOH_ * ((mmf[j] - zeta_s) / (1 - zeta_s));
        }
    }

    // NaOH

    for (int i = 0; i <= N_ - 1; i++)
    {
        OH_[i][0] = 0;
    }

    for (int j = 0; j <= N_ - 1; j++)
    {
        OH_[N_ - 1][j] = NaOH_ * mmf[j];
    }

    for (int i = 0; i <= N_ - 1; i++)
    {
        OH_[i][N_ - 1] = NaOH_;
    }

    // Print basic component look-up table
    if (r_ == 1)
    {
        Info << "NaOH look-up table" << endl
             << "---------------------" << endl;
    }
    else if (r_ == 2)
    {
        Info << "Na+ look-up table" << endl
             << "---------------------" << endl;
    }
    for (int i = 0; i <= N_ - 1; i++)
    {
        for (int j = 0; j <= N_ - 1; j++)
        {
            Info << OH_[i][j] << " ";
        }
        Info << endl;
    }
    Info << "\n\n";
}

void Foam::chemReac::update()
{
    double OF_mmf = 0;
    double OF_Is = 0;

    forAll(alpha1_, Celli)
    {

        if (alpha1_[Celli] >= 0.0 && alpha1_[Celli] <= 1.0)
        {
            OF_mmf = alpha1_[Celli];
        }
        else if (alpha1_[Celli] < 0.0)
        {
            OF_mmf = 0.0;
        }
        else if (alpha1_[Celli] > 1.0)
        {
            OF_mmf = 1.0;
        }

        if (var_[Celli] >= 0.0 && var_[Celli] <= OF_mmf * (1 - OF_mmf))
        {
            OF_Is = var_[Celli] / var_max;
        }
        if (var_[Celli] > OF_mmf * (1 - OF_mmf))
        {
            OF_Is = 1.0;
        }

        int a0, a1, b0, b1;
        a0 = 0;
        a1 = 0;

        for (int i = 0; i <= N_ - 2; i++)
        {
            if (OF_mmf <= mmf[i + 1] && OF_mmf >= mmf[i])
            {
                a0 = i;
                a1 = i + 1;
            }
        }

        b0 = 0;
        b1 = 0;
        for (int i = 0; i < N_ - 1; i++)
        {
            if (OF_Is >= Is[i][a0] && OF_Is <= Is[i + 1][a0])
            {
                b0 = i;
                b1 = i + 1;
            }
        }

        // Bilinear interpolation for NaOH
        double Na1 = (OH_[b1][a0] - OH_[b0][a0]) / (Is[b1][a0] - Is[b0][a0]) * (OF_Is - Is[b0][a0]) + OH_[b0][a0];
        double Na2 = (OH_[b1][a1] - OH_[b0][a1]) / (Is[b1][a1] - Is[b0][a1]) * (OF_Is - Is[b0][a1]) + OH_[b0][a1];
        OH[Celli] = ((Na2 - Na1) / (mmf[a1] - mmf[a0]) * (OF_mmf - mmf[a0]) + Na1);
    }
}

// ************************************************************************* //