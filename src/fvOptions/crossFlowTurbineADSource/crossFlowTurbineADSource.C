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

\*---------------------------------------------------------------------------*/

#include "crossFlowTurbineADSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(crossFlowTurbineADSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        crossFlowTurbineADSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineADSource::crossFlowTurbineADSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    crossFlowTurbineALSource(name, modelType, dict, mesh)
{
    baseCustomTime_ = mesh.time().value();
    rotateAD(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineADSource::~crossFlowTurbineADSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::crossFlowTurbineADSource::rotateAD(bool updateOnly)
{
    scalar radians = 2*mathematical::pi/divisions_;
    customDeltaT_ = radians/omega_;

    //updateOnly is intended for first step only to set custom speed
    if (updateOnly == false)
    {
        customTime_[azimuthIndex_] = baseCustomTime_;
        baseCustomTime_ += customDeltaT_;
        rotate(radians);
        angleDeg_[azimuthIndex_] = baseAngleDeg_;
        baseAngleDeg_ += radToDeg(radians);
        //lastRotationTime_ = time_.value();
    }
    updateTSROmega();
    
    forAll(actuatorLines_, i)
    {
        actuatorLines_[i]->setCustomTime(baseCustomTime_, customDeltaT_);
    }
}

void Foam::fv::crossFlowTurbineADSource::rotate(scalar radians)
{
    if (debug)
    {
        Info<< "Rotating " << name_ << " " << radians << " radians"
            << endl << endl;
    }

    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, axis_, radians);
        blades_[i].setSpeed(origin_, axis_, omega_);
    }

    if (hasStruts_)
    {
        forAll(struts_, i)
        {
            struts_[i].rotate(origin_, axis_, radians);
            struts_[i].setSpeed(origin_, axis_, omega_);
        }
    }

    if (hasShaft_)
    {
        shaft_->rotate(origin_, axis_, radians);
        shaft_->setSpeed(origin_, axis_, omega_);
    }
}

void Foam::fv::crossFlowTurbineADSource::initializeAL()
{
    // if compactField is false, we do not need influenceCells
    if (compactField_ == false)
    {
        return;
    }

    forAll(actuatorLines_, i)
    {
        for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
        {
            forAll(actuatorLines_[i]->elements(), j)
            {
                actuatorLines_[i]->setAzimuthIndex(azimuthIndex_, false);
                actuatorLines_[i]->elements()[j].calcInfluenceEpsilon
                (
                    maxDragCoefficient_
                );
            }
        }
    }
    distributeEpsilon();

    labelList globalToLocal(mesh_.nCells(), -1);

    label nActive = 0;
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        // Add scalar source term from blades
        forAll(actuatorLines_, i)
        {
            actuatorLines_[i]->setAzimuthIndex(azimuthIndex_, false);
            actuatorLines_[i]->constructInfluenceCellList
            (
                azimuthIndex_,
                globalToLocal,
                nActive
            );
        }
        rotateAD();
    }
    label localCells = mesh_.nCells();
    label nCellsGlobal = localCells;
    reduce(nCellsGlobal, sumOp<label>());

    label nActiveGlobal = nActive;
    reduce(nActiveGlobal, sumOp<label>());

    // Print only once
    if (Pstream::master())
    {
        Info<< "Active cells participating in the force field: "
            << nActiveGlobal << " of " << nCellsGlobal << endl;
    }

    activePositions_.setSize(nActive);
    activeForceField_.setSize(nActive, Zero);
    localToGlobal_.setSize(nActive);

    const vectorField& C = mesh_.C();
    forAll(globalToLocal, globalI)
    {
        label localI = globalToLocal[globalI];

        if (localI != -1)
        {
            activePositions_[localI] = C[globalI];
            localToGlobal_[localI] = globalI;
        }
    }

    forAll(actuatorLines_, i)
    {
        actuatorLines_[i]->setCompactFields
        (
            activePositions_,
            activeForceField_
        );
    }
    
    baseAngleDeg_ = 0;
    angleDeg_[0] = 0;
}

void Foam::fv::crossFlowTurbineADSource::setupPositions(bool includeRing)
{
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        forAll(actuatorLines_, i)
        {
            actuatorLines_[i]->setAzimuthIndex(azimuthIndex_, false);
            actuatorLines_[i]->findCells(includeRing);
        }
        rotateAD();
    }
}

void Foam::fv::crossFlowTurbineADSource::calculateForces()
{
    // code can run extra revolutions to make dynamic stall converge
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
        {
            forAll(actuatorLines_, i)
            {
                actuatorLines_[i]->setAzimuthIndex(azimuthIndex_, false);
                actuatorLines_[i]->setCustomTime
                (
                    customTime_[azimuthIndex_],
                    customDeltaT_
                );
                actuatorLines_[i]->calculateElementForces();
            }
        }
    }
}

void Foam::fv::crossFlowTurbineADSource::addForce
(
    volVectorField &forceField,
    scalar scale,
    bool compressible
)
{
    // forceField should be the average during one revolution here
    if (compactField_)
    {
        if (activeForceField_.size() > 0)
        {
            activeForceField_ = vector::zero;
        }
    }

    meanPowerCoefficient_ = 0;
    meanDragCoefficient_ = 0;
    meanTorqueCoefficient_ = 0;
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        // Zero out force vector and field
        force_ *= 0;

        // Create local moment vector
        vector moment(vector::zero);

        // Add source for blade actuator lines
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_, true);
            blades_[i].setCustomTime // Not needed in current implementation
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            blades_[i].addForceFromChild
            (
                forceField,
                static_cast<scalar>(bladeMultiplier_)/divisions_,
                compressible
            );

            bladeMoments_[i] = blades_[i].moment(origin_);
            
            // when using nBlades == 1 with bladeMultiplier
            // emulate 3 blades with even spacing
            if (bladeMultiplier_ > 1 && nBlades_ == 1)
            {
                for (label k = 0; k < bladeMultiplier_; k++)
                {
                    label newazimuthIndex =
                        (azimuthIndex_ + divisions_/bladeMultiplier_*k)
                        % divisions_;
                    blades_[i].setAzimuthIndex(newazimuthIndex, false);
                    force_ += blades_[i].force();
                    moment += blades_[i].moment(origin_);
                }
                blades_[i].setAzimuthIndex(azimuthIndex_, false);
            }
            else
            {
                force_ += bladeMultiplier_*blades_[i].force();
                moment += bladeMultiplier_*bladeMoments_[i];
            }
        }

        if (hasStruts_)
        {
            // Add source for strut actuator lines
            forAll(struts_, i)
            {
                struts_[i].setAzimuthIndex(azimuthIndex_, true);
                struts_[i].setCustomTime
                (
                    customTime_[azimuthIndex_],
                    customDeltaT_
                );
                struts_[i].addForceFromChild
                (
                    forceField,
                    static_cast<scalar>(bladeMultiplier_)/divisions_,
                    compressible
                );

                // when using nBlades == 1 with bladeMultiplier
                // emulate 3 blades with even spacing
                // Note: code will assume that if you only give one blade
                // you also only gave the struts for one blade
                // checking nStruts == 1 does not work as one blade can have
                // multiple struts
                if (bladeMultiplier_ > 1 && nBlades_ == 1)
                {
                    for (label k = 0; k < bladeMultiplier_; k++)
                    {
                        label newazimuthIndex =
                            (azimuthIndex_ + divisions_/bladeMultiplier_*k)
                            % divisions_;
                        struts_[i].setAzimuthIndex(newazimuthIndex, false);
                        force_ += struts_[i].force();
                        moment += struts_[i].moment(origin_);
                    }
                    struts_[i].setAzimuthIndex(azimuthIndex_, false);
                }
                else
                {
                    force_ += bladeMultiplier_*struts_[i].force();
                    moment += bladeMultiplier_*struts_[i].moment(origin_);
                }
            }
        }

        if (hasShaft_)
        {
            // Add source for shaft actuator line
            shaft_->setAzimuthIndex(azimuthIndex_, true);
            shaft_->setCustomTime
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            shaft_->addForceFromChild
            (
                forceField,
                1.0/divisions_,
                compressible
            );
            force_ += shaft_->force();
            moment += shaft_->moment(origin_);
        }

        // Torque is the projection of the moment from
        // all blades on the axis
        torque_ = moment & axis_;

        torqueCoefficient_ =
            torque_/(0.5*frontalArea_*rotorRadius_
            * magSqr(freeStreamVelocity_));
        powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
        dragCoefficient_ =
            force_ & freeStreamDirection_
            / (0.5*frontalArea_*magSqr(freeStreamVelocity_));


        meanPowerCoefficient_ += powerCoefficient_;
        meanDragCoefficient_ += dragCoefficient_;
        meanTorqueCoefficient_ += torqueCoefficient_;
        // Print performance to terminal
        if (printPerf_)
        {
            printPerf();
        }

        // Write performance data
        // Note this will write multiples if there are
        // multiple PIMPLE loops
        if (Pstream::master())
        {
            writePerf();
        }
    }
    meanPowerCoefficient_ /= divisions_;
    meanDragCoefficient_ /= divisions_;
    meanTorqueCoefficient_ /= divisions_;
    
    // When using compressed fields, restore them to the original forceField
    updateForceFieldAD(forceField);
}

void Foam::fv::crossFlowTurbineADSource::addTurbulence
(
    fvMatrix<scalar>& eqn,
    const word fieldName
)
{
    // kField should be the average during one revolution here
    fvMatrix<scalar> kField(eqn.psi(), eqn.dimensions());
    kField *= dimensionedScalar("zero", eqn.dimensions(), 0.0);
    fvMatrix<scalar> kFieldShaft(eqn.psi(), eqn.dimensions());
    kFieldShaft *=
        dimensionedScalar("zero", eqn.dimensions(), 0.0);
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        // Add scalar source term from blades
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_, false);
            blades_[i].setCustomTime // Not needed in current implementation
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            blades_[i].addTurbulence(kField, fieldName);
        }

        if (hasStruts_)
        {
            // Add source for strut actuator lines
            forAll(struts_, i)
            {
                struts_[i].setAzimuthIndex(azimuthIndex_, false);
                struts_[i].setCustomTime
                (
                    customTime_[azimuthIndex_],
                    customDeltaT_
                );
                struts_[i].addTurbulence(kField, fieldName);
            }
        }

        if (hasShaft_)
        {
            // Add source for shaft actuator line
            shaft_->setAzimuthIndex(azimuthIndex_, false);
            shaft_->setCustomTime
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            shaft_->addTurbulence(kFieldShaft, fieldName);
        }
    }
    eqn += (static_cast<scalar>(bladeMultiplier_)/divisions_)*kField
            + (1.0/divisions_)*kFieldShaft;
}

void Foam::fv::crossFlowTurbineADSource::allocateAL()
{
    angleDeg_.setSize(divisions_);
    customTime_.setSize(divisions_);
    forAll(actuatorLines_, i)
    {
        actuatorLines_[i]->allocateInfluenceCells
        (
            divisions_,
            cacheInteractions_
        );
    }
    label nEpsilon = 0;
    forAll(actuatorLines_, i)
    {
        forAll(actuatorLines_[i]->elements(), j)
        {
            nEpsilon += actuatorLines_[i]->elements()[j].epsilonCount();
        }
    }
    // if-statement should not be necessary as this runs before
    // actuatorModelBase changes its size for the farm case.
    if (nEpsilon > epsilon_.size())
    {
        // When running turbineFarmSource, epsilon_ is not allocated
        // in actuatorModelBase for this class
        epsilon_.resize(nEpsilon);
    }
}

// ************************************************************************* //
