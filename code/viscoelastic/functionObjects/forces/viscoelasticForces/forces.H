/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::viscoelasticForces

Group
    grpForcesFunctionObjects

Description
    Calculates the viscoelastic forces and moments by integrating the pressure
	and skin-friction forces over a given list of patches.

    Forces and moments are calculated in a global Cartesian coordinate system
    by default, or using a user-defined system.

    Results are written to multiple files as a function of time in the
    postProcessing/\<functionObjectName\> directory:
    - force.dat          : viscoelasticForces
    - moment.dat         : moments

Usage
    Example of function object specification:
    \verbatim
    viscoelasticForces1
    {
        type        viscoelasticForces;
        libs        ("$(FOAM_USER_LIBBIN)/libviscoelasticTransportModels.so");
        ...
        log         yes;
        writeFields yes;
        patches     (walls);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: viscoelasticForces       | yes         |
        log          | Write force data to standard output | no | no
        writeFields  | Write the force and moment fields | no | no
        patches      | Patches included in the viscoelasticForces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        CofR         | Centre of rotation (see below) | no   |
    \endtable

Note
  - The centre of rotation (CofR) for moment calculations can either be
    specified by an \c CofR entry, or be taken from origin of the local
    coordinate system.  For example,
    \verbatim
        CofR        (0 0 0);
    \endverbatim
    or
    \verbatim
        origin  (0 0 0);
        e1      (0 1 0);
        e3      (0 0 1);
    \endverbatim
    or
    \verbatim
        coordinateSystem
        {
            origin  (0 0 0);
            rotation
            {
                type    axes;
                e3      (0 0 1);
                e1      (1 0 0);
            }
        }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    viscoelasticForces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_viscoelasticForces_H
#define functionObjects_viscoelasticForces_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "cartesianCS.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class viscoelasticForces Declaration
\*---------------------------------------------------------------------------*/

class viscoelasticForces
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Pressure, viscous and porous force per bin
        List<Field<vector>> force_;

        //- Pressure, viscous and porous moment per bin
        List<Field<vector>> moment_;

        // File streams

            //- Forces
            autoPtr<OFstream> forceFilePtr_;

            //- Moments
            autoPtr<OFstream> momentFilePtr_;


        // Read from dictionary

            //- Patches to integrate viscoelasticForces over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Reference pressure
            scalar pRef_;

            //- Coordinate system used when evaluating forces/moments
            coordSystem::cartesian coordSys_;

            //- Write fields flag
            bool writeFields_;

            //- Initialised flag
            bool initialised_;


    // Protected Member Functions

        //- Create a field name
        word fieldName(const word& name) const;

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

        //- Set the co-ordinate system from dictionary and axes names
        void setCoordinateSystem
        (
            const dictionary& dict,
            const word& e3Name = word::null,
            const word& e1Name = word::null
        );

        //- Initialise the fields
        void initialise();

        //- Reset the fields prior to accumulation of force/moments
        void resetFields();

        //- Return the viscoelastic stress
        tmp<volSymmTensorField> tau() const;

        //- Return density
        dimensionedScalar rho() const;

        //- Add patch contributions to force and moment fields
        void addToFields
        (
            const label patchi,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT
        );

        //- Add cell contributions to force and moment fields
        void addToFields
        (
            const labelList& cellIDs,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT
        );

        //- Helper function to write integrated forces and moments
        void writeIntegratedForceMoment
        (
            const string& descriptor,
            const vectorField& fm0,
            const vectorField& fm1,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write force data
        void writeForces();

        //- No copy construct
        viscoelasticForces(const viscoelasticForces&) = delete;

        //- No copy assignment
        void operator=(const viscoelasticForces&) = delete;


public:

    //- Runtime type information
    TypeName("viscoelasticForces");


    // Constructors

        //- Construct from Time and dictionary
        viscoelasticForces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        viscoelasticForces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~viscoelasticForces() = default;


    // Member Functions

        //- Read the viscoelasticForces data
        virtual bool read(const dictionary&);

        //- Calculate the viscoelasticForces and moments
        virtual void calcForcesMoment();

        //- Return the total force
        virtual vector forceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the viscoelasticForces
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
