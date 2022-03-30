/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::pBar

Group
    grpFieldFunctionObjects

Description
    Computes the Peclet number as a \c surfaceScalarField.

    Operands:
    \table
      Operand        | Type               | Location
      input          | surfaceScalarField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -                  | -
      output field   | surfaceScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    pBar1
    {
        // Mandatory entries (unmodifiable)
        type            pBar;
        libs            (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        rho             rho;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: pBar                | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      rho          | Name of density field              | word |  no   | rho
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func pBar
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression

SourceFiles
    pBar.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pBar_H
#define functionObjects_pBar_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class pBar Declaration
\*---------------------------------------------------------------------------*/

class pBar
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the dimensionless pressure field and return true
		//  if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("pBar");


    // Constructors

        //- Construct from Time and dictionary
        pBar
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        pBar(const pBar&) = delete;

        //- No copy assignment
        void operator=(const pBar&) = delete;


    //- Destructor
    virtual ~pBar() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
