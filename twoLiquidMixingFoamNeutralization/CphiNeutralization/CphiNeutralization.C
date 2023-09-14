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
    Foam::Cphi

Description

SourceFiles
    CphiNeutralization.C

\*---------------------------------------------------------------------------*/

#include "CphiNeutralization.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Cphi::Cphi(
    const fvMesh &mesh,
    const incompressible::momentumTransportModel &turbulence)

    : IOdictionary(

          IOobject(
              "chemicalReactionProperties",
              mesh.time().constant(),
              mesh,
              IOobject::MUST_READ_IF_MODIFIED,
              IOobject::NO_WRITE)),

      mesh_(mesh),

      turbulence_(turbulence),

      Cf_(readScalar(lookup("Cf"))),

      gamma_(
          IOobject(
              "Cphi",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("gamma", dimless, 0.0)),

      nu_(
          IOobject(
              "nu",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("nu_", dimViscosity, 1e-06))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Cphi::~Cphi()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::Cphi::calculate()
{

    // forAll(mesh_.C(), i)
    // {
    //     Re_l = turbulence_.k()()[i] / Foam::pow(turbulence_.epsilon()()[i] * nu_[i], 0.5);
    //     gamma_[i] = 0.0;

    //     for (int j = 0; j < 7; j++)
    //     {
    //         gamma_[i] = (gamma_[i] + a[j] * Foam::pow(Foam::log10(Re_l), j));
    //     }
        
    //     gamma_[i] =  Cf_ * gamma_[i];
    // }
    gamma_ = Cf_;
}

// ************************************************************************* //