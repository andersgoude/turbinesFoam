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

        //- number of divisions per revolution
        label divisions_;
        
        //- Run the force calculation loop additional times to allow
        // dynamic stall models with hysteresis develop the solution
        // 1 additional loop should be sufficient
        label dynStallLoop_;
        
        //- Multiply the force with this value (basically allows to set up
        // one bladed turbines that act like 3 bladed turbines)
        scalar bladeMultiplier_;

        //- Compress the interaction field to improve memory acces
        bool compactField_;

        //- Cache intearction point list to avoid unnecessary iterations
        bool cacheInteractions_;
        
        //-first time step indicator, used to determine when to build cells
        bool firstUse_;

        //- Cached positions for more compact disc operations
        labelList localToGlobal_;

        //- Cached positions for more compact disc operations
        vectorField activePositions_;

        //- ForceField corresponding to the compact field
        vectorField activeForceField_;

        //- ForceField after relaxation is applied
        vectorField filteredForceField_;

        //- Apply relaxation to force field
        bool relaxForceField_;

        //- Starting iteration for relaxation
        label relaxStartTime_;

        //- Relaxation value (note equation is 1/(1+relaxValue_))
        scalar relaxValue_;

        //- Relaxation growth factor (relaxation += relaxGrowthValue_)
        scalar relaxGrowthValue_;

        //- Only grow if relative difference is larger than threshold
        scalar relaxGrowthThreshold_;

        //- Peak value for relaxation
        scalar relaxMaxValue_;
        
        
    // Protected Member Functions

        //- Special rotation function for actuator disc
        void rotateAD(bool updateOnly = false);

        //- Apply the forcefield with relaxation
        void updateForceField(volVectorField &forceField);

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

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
