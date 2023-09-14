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
    Cphi.C

\*---------------------------------------------------------------------------*/

#ifndef CphiNeutralization_H
#define CphiNeutralization_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "geometricZeroField.H"
#include "kinematicMomentumTransportModel.H"

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                        Class Cphi Declaration
    \*---------------------------------------------------------------------------*/

    class Cphi : public IOdictionary
    {
        // Private data

        //- Reference to mesh
        const fvMesh &mesh_;

        //- Reference to the turbulence
        const incompressible::momentumTransportModel &turbulence_;

        //- DimensionedField Cf
        scalar Cf_;

        //- DimensionedField Cphi
        volScalarField gamma_;

        //- DimensionedField nu_
        DimensionedField<scalar, volMesh> nu_;

        //- Coefficients for Cphi calcultion
        scalar a[7] = {0.4093, 0.6015, 0.5851, 0.09472, -0.3909, 0.1461, -0.01604};

        //- Local Reynolds number
        scalar Re_l;

    public:
        // Constructors

        //- Construct from fvMesh
        Cphi(const fvMesh &mesh, const incompressible::momentumTransportModel &turbulence);

        //- Destructor
        virtual ~Cphi();

        // Member Functions
        void calculate();

        // Access function

        //- Return the Cphi
        inline const volScalarField &CphiValue() const;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

inline const Foam::volScalarField &
Foam::Cphi::CphiValue() const
{
    return gamma_;
}

#endif