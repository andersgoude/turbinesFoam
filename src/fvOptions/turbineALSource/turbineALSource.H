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
    Foam::fv::turbineALSource

Description
    Cell based momentum source, which is a collection of actuatorLineSources

SourceFiles
    turbineALSource.C

\*---------------------------------------------------------------------------*/

#ifndef turbineALSource_H
#define turbineALSource_H

#include "cellSetOption.H"
//#include "NamedEnum.H
#include "actuatorModelBase.H"
#include "actuatorLineSource.H"
#include "volFieldsFwd.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class turbineALSource Declaration
\*---------------------------------------------------------------------------*/

class turbineALSource
:
    public actuatorModelBase
{

protected:

    // Protected data

        //- Runtime pointer
        const Time& time_;

        //- Time value used to determine whether to rotate
        scalar lastRotationTime_;

        //- Origin of turbine coordinate system
        vector origin_;

        //- Turbine axis of rotation
        vector axis_;

        //- Reference density for incompressible case
        scalar rhoRef_;

        //- Rotational speed in rad/s
        //  Positive anti-clockwise when looking along axis direction
        scalar omega_;

        //- Azimuthal angle in degrees (does not reset each rotation)
        scalar baseAngleDeg_;

        //- List of angles for each position
        List<scalar> angleDeg_;

        //- internal time for actuator disc (for use by dynamic stall model)
        //- used to perform a loop within the openfoam-step
        scalar baseCustomTime_;

        //- List of custom time values for actuator disc
        List<scalar> customTime_;

        //- deltaT value associated with customTime
        scalar customDeltaT_;

        // current position for the actuator disc
        label azimuthIndex_;

        //- Number of blades
        label nBlades_;

        //- Inlet velocity for specified inflow
        vector freeStreamVelocity_;

        //- Free stream velocity direction
        vector freeStreamDirection_;

        //- (Initial) radial direction
        vector radialDirection_;

        //- Rotor tip speed ratio
        scalar tipSpeedRatio_;

        //- Rotor radius for calculating tip speed ratio
        scalar rotorRadius_;

        //- List of actuator lines as blades
        PtrList<actuatorLineSource> blades_;

        //- Blade definitions dictionary
        dictionary bladesDict_;

        //- List of blade names
        wordList bladeNames_;

        //- Dictionary for profile coefficient data
        dictionary profileData_;

        //- Force vector
        vector force_;

        //- Torque about the axis
        scalar torque_;

        //- Frontal area in m^2
        scalar frontalArea_;

        //- Nondimensional power coefficient
        scalar powerCoefficient_;

        //- Nondimensional drag coefficient
        scalar dragCoefficient_;

        //- Nondimensional torque coefficient
        scalar torqueCoefficient_;

        // Actuator disc specific versions
        //- Nondimensional power coefficient
        scalar meanPowerCoefficient_;

        //- Nondimensional drag coefficient
        scalar meanDragCoefficient_;

        //- Nondimensional torque coefficient
        scalar meanTorqueCoefficient_;

        //- Output file stream
        OFstream* outputFile_;

        //- Dynamic stall dictionary
        dictionary dynamicStallDict_;

        //- Amplitude of prescribed TSR oscillation
        scalar tsrAmplitude_;

        //- Angle of first peak in TSR oscillation (rad)
        scalar tsrPhase_;

        //- Mean tip speed ratio
        scalar meanTSR_;

        //- Individual blade moments about turbine origin
        List<vector> bladeMoments_;

    // variables used to automatically detect active region for
    // actuator line simulations
    
        //- Maximum possible drag coefficient for cell list creation
        scalar maxDragCoefficient_;

        //- lower axial point of cylinder (relative to center)
        scalar cylMin_;

        //- upper axial point of cylinder (relative to center)
        scalar cylMax_;

        //- cylinder ratiud
        scalar cylRadius_;

        //- Maximum chord length of any element
        scalar chordMax_;

        //- largest value of chord length * chordFactor for any element
        scalar epsilonLiftMax_;

        //- largest value of chord length * dragFactor for any element
        scalar dragMax_;

        //- largest value of meshFactor for any element
        scalar meshFactorMax_;


    // Actuator disc variables

        //- number of divisions per revolution
        label divisions_;
        
        //- Run the force calculation loop additional times to allow
        // dynamic stall models with hysteresis develop the solution
        // 1 additional loop should be sufficient
        label dynStallLoop_;
        
        //- Multiply the force with this value (basically allows to set up
        // one bladed turbines that act like 3 bladed turbines)
        label bladeMultiplier_;

        //- Compress the interaction field to improve memory acces
        bool compactField_;

        //- Cache intearction point list to avoid unnecessary iterations
        bool cacheInteractions_;

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

        //- Rotate a vector
        void rotateVector
        (
            vector& vectorToRotate,
            vector rotationPoint,
            vector axis,
            scalar radians
        );

        //- Calculate the cone angle from elementdata
        scalar calculateCone
        (
            const List<List<scalar>>& elementData,
            const label j
        );

        //- Create actuator lines for blades
        virtual void createBlades();

        //- Create the co-ordinate system
        virtual void createCoordinateSystem();

        //- Update tip speed ratio and angular velocity
        virtual void updateTSROmega();

        //- Rotate the turbine according to its angular velocity and time step
        virtual void rotate();

        //- Rotate the turbine a specified angle about its axis
        virtual void rotate(scalar radians);

        //- Print performance
        virtual void printPerf();

        //- Determine minimum cylinder that fits all actuator line elements
        void determineBoundingCylinder(PtrList<actuatorLineSource> &lines);

        //- Determine maximum cell volume within the cylinder
        scalar maxCellVolumeInCylinder();

        //- Create the list of influence cells for one angle
        void addPointsInCylinder
        (
            scalar maxVol,
            labelList& globalToLocal,
            label& nActive
        );

        //- Apply the forcefield for actuator line
        void updateForceFieldAL(volVectorField &forceField);

        //- Apply the forcefield for actuator disc with relaxation
        void updateForceFieldAD(volVectorField &forceField);


public:

    //- Runtime type information
    TypeName("turbineALSource");


    // Constructors

        //- Construct from components
        turbineALSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~turbineALSource();


    // Member Functions

        // Access

            //- Return const access to runTime
            inline const Time& time() const;

            //- Return the reference density for incompressible case
            inline scalar rhoRef() const;


        // Source term addition


        // I-O

            //- Create the turbine output file
            virtual void createOutputFile();

            //- Write the turbine performance to file
            virtual void writePerf();

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

            //- Print dictionary values
            virtual void printCoeffs() const;

            // return power cofficinet
            scalar powerCoefficient() const;
            
            // return drag coefficient
            scalar dragCoefficient() const;

            // return power cofficinet
            scalar meanPowerCoefficient() const;
            
            // return drag coefficient
            scalar meanDragCoefficient() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
