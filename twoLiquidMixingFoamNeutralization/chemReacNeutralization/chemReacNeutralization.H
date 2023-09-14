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

#ifndef chemReac_H
#define chemReac_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "geometricZeroField.H"
#include "kinematicMomentumTransportModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class chemReac Declaration
\*---------------------------------------------------------------------------*/

    class chemReac : public IOdictionary
    {
        // Private data

        //- Reference to mesh
        const fvMesh &mesh_;

        //- Reference to alpha1
        const volScalarField &alpha1_;

        //- Reference to var
        const volScalarField &var_;

        //- Reference to the flux at faces
        const surfaceScalarField &phi_;

        //- Reference to the turbulence
        const incompressible::momentumTransportModel &turbulence_;

        // Read the reciprocal of the turbulent Schmidt number
        const dimensionedScalar &alphatab_;

        // Read the diffusivity of ion
        const dimensionedScalar &D_;

        //-Number of nodes for look-up table construction
        static const label N_ = 11;

        //-Max variance value
        const scalar var_max = 0.25+1e-6;

        //- HCl concentration
        scalar HCl_;

        //- NaOH concentration
        scalar NaOH_;

        //- water density
        scalar rhow_;

        //- Stoichiometric ratio (A+rB)
        scalar r_;

        //- Stoichiometric mixture fraction
        scalar zeta_s;

        //- Mixture fraction array
        scalar mmf[N_];

        //- Segregation index 2D array
        scalar Is[N_][N_];

        //- Look-up table of B component (OH-)
        scalar OH_[N_][N_];

        //- OH concentration field
        volScalarField OH;

        //- Dimension of concentration
        const dimensionSet dimConc = dimMoles / dimVol;

    public:
        // Constructors

        //- Construct from fvMesh
        chemReac(const fvMesh &mesh, const volScalarField &alpha1, const volScalarField &var, const surfaceScalarField &phi, const incompressible::momentumTransportModel &turbulence, const dimensionedScalar &alphatab, const dimensionedScalar &Dab);

        //- Destructor
        virtual ~chemReac();

        // Member Functions

        //- Find the beta function value
        double incbeta(double, double, double);

        //- Generate the mixture fraction array
        void mixtureFraction();

        //- Generate the segregation coefficient 2Darray
        void segregationMatrix();

        //- Generate the concentration matrix of component B (Na+)
        void lookUpTableB();

        //- Update the volScalarFields knowing alpha and variance
        void update();

        // Access function

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif