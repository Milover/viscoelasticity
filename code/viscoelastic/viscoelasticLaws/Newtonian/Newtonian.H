/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Newtonian

SourceFiles
    Newtonian.C

\*---------------------------------------------------------------------------*/

#ifndef Newtonian_H
#define Newtonian_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Newtonian Declaration
\*---------------------------------------------------------------------------*/

class Newtonian
:
    public viscoelasticLaw
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar eta_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Newtonian(const Newtonian&);

        //- Disallow default bitwise assignment
        void operator=(const Newtonian&);

public:

    //- Runtime type information
    TypeName("Newtonian");

    // Constructors

        //- Construct from components
        Newtonian
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~Newtonian()
        {}


    // Member Functions

		//- Return density
		virtual dimensionedScalar rho() const
		{
			return rho_;
		}

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const
        {
            return tau_;
        }

        //- Return the viscosity
        virtual tmp<volScalarField> eta() const;

        //- Return the solvent viscosity
        virtual tmp<volScalarField> etaS() const
		{
			return eta();
		}

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(const volVectorField& U) const;

        //- Correct the viscoelastic stress
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
