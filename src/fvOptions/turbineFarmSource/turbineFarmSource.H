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
    Foam::fv::actuatorModelBase

Description
    Base class for all actuator models, with focus on handling parallel
    computing efficiently

SourceFiles
    actuatorModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef turbineFarmSource_H
#define turbineFarmSource_H

#include "actuatorModelBase.H"
#include "axialFlowTurbineALSource.H"
#include "crossFlowTurbineALSource.H"
//#include "NamedEnum.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "interpolationCellPoint.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class turbineFarmSource  Declaration
\*---------------------------------------------------------------------------*/

class turbineFarmSource
:
    public actuatorModelBase
{

protected:

    // Protected data

        //- Axial flow turbines
        PtrList<axialFlowTurbineALSource> axialFlowTurbines_;

        //- Cross flow turbines
        PtrList<crossFlowTurbineALSource> crossFlowTurbines_;

        //- Axial flow turbine definitions dictionary
        dictionary axialFlowTurbinesDict_;

        //- Cross flow turbine definitions dictionary
        dictionary crossFlowTurbinesDict_;

        //- List of turbine names
        wordList axialFlowTurbineNames_;

        //- List of turbine names
        wordList crossFlowTurbineNames_;

        //- Number of axial flow turbines
        label nAxialFlowTurbines_;

        //- Number of axial flow turbines
        label nCrossFlowTurbines_;

        //- Individual axial flow turbine power
        List<scalar> axialFlowTurbinePowerCoefficients_;

        //- Individual axial flow turbine thrust coefficients
        List<scalar> axialFlowTurbineDragCoefficients_;

        //- Individual cross flow turbine power
        List<scalar> crossFlowTurbinePowerCoefficients_;

        //- Individual cross flow turbine thrust coefficients
        List<scalar> crossFlowTurbineDragCoefficients_;

        //- Output file stream
        OFstream* outputFile_;

    // Protected Member Functions

        //- Create axial flow turbines
        void createAxialFlowTurbines();

        //- Create cross flow turbines
        void createCrossFlowTurbines();

        //- Collect power and thrust data from turbines
        void getTurbineData();

public:

    //- Runtime type information
    TypeName("turbineFarmSource");


    // Constructors

        //- Construct from components
        turbineFarmSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~turbineFarmSource();


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

            //- Virtual function that each class implements to create
            //  the field for itself + all child objects
            virtual void createForceFieldForChildren(const bool compressible);

        // Access


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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
