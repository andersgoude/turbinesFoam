/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

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
    Foam::fv::crossFlowTurbineADSource

Description
    Cell based momentum source that represents a cross-flow turbine

SourceFiles
    crossFlowTurbineADSource.C

\*---------------------------------------------------------------------------*/

#ifndef crossFlowTurbineADSource_H
#define crossFlowTurbineADSource_H

#include "crossFlowTurbineALSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class crossFlowTurbineALSource Declaration
\*---------------------------------------------------------------------------*/

class crossFlowTurbineADSource
:
    public crossFlowTurbineALSource
{

protected:

    // Protected data

    // Protected Member Functions

        //- Special rotation function for actuator disc
        void rotateAD(bool updateOnly = false);

        //- Rotate the turbine a specified angle about its axis
        virtual void rotate(scalar radians);


public:

    //- Runtime type information
    TypeName("crossFlowTurbineADSource");


    // Constructors

        //- Construct from components
        crossFlowTurbineADSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~crossFlowTurbineADSource();


    // Member Functions

        // Routines for actuatorModelBase for parallel efficiency

            //- Alloction part of the initial setup
            virtual void allocateAL();

            //- Work to be done first time for calculateALData
            virtual void initializeAL();

            //- Set up positions by calling findCells on all
            //  actuator line elements
            virtual void setupPositions(bool includeRing);

            //- Calculate the force once velocity is sampled
            virtual void calculateForces();

        // Source term addition

            virtual void addForce
            (
                volVectorField& forceField,
                scalar scale,
                bool compressible
            );

            //- Add source term to turbulence model equation
            virtual void addTurbulence
            (
                fvMatrix<scalar>& eqn,
                const word fieldName
            );

        // I-O

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
