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

#include "axialFlowTurbineADSource.H"
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
    defineTypeNameAndDebug(axialFlowTurbineADSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        axialFlowTurbineADSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::axialFlowTurbineADSource::axialFlowTurbineADSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    axialFlowTurbineALSource(name, modelType, dict, mesh)
{
    customTime_ = mesh.time().value();
    rotateAD(true);

    // override the nBlades value for the end effects calculation
    // if bladeMultiplier is used
    effectiveNBlades_ = bladeMultiplier_*nBlades_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::fv::axialFlowTurbineADSource::~axialFlowTurbineADSource()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::axialFlowTurbineADSource::rotateAD(bool updateOnly)
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

    //Info << "rotateAD called for time " << time_.value()
    //     << " custom time: " << customTime_ << endl;
    forAll(actuatorLines_, i)
    {
        actuatorLines_[i]->setCustomTime(baseCustomTime_, customDeltaT_);
    }
}

void Foam::fv::axialFlowTurbineADSource::rotate(scalar radians)
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

    if (hasHub_)
    {
        hub_->rotate(origin_, axis_, radians);
        hub_->setSpeed(origin_, axis_, omega_);
    }
}

void Foam::fv::axialFlowTurbineADSource::allocateAL()
{
    angleDeg_.setSize(divisions_);
    customTime_.setSize(divisions_);
    forAll(blades_, i)
    {
        blades_[i].allocateInfluenceCells(divisions_, cacheInteractions_);
    }

    if (hasHub_)
    {
        // Add source for hub actuator line
        hub_->allocateInfluenceCells(divisions_, cacheInteractions_);
    }

    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->allocateInfluenceCells(1, cacheInteractions_);
    }

    if (hasNacelle_)
    {
        // Add source for nacelle actuator line
        nacelle_->allocateInfluenceCells(1, cacheInteractions_);
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

void Foam::fv::axialFlowTurbineADSource::initializeAL()
{
    // if compactField is false, we do not need influenceCells
    if (compactField_ == false)
    {
        return;
    }

    if (hasTower_)
    {
        tower_->setAzimuthIndex(0, false); // not really needed, remove later
        tower_->calcInfluenceEpsilon(maxDragCoefficient_);
    }

    if (hasNacelle_)
    {
        nacelle_->setAzimuthIndex(0, false); // not really needed, remove later
        nacelle_->calcInfluenceEpsilon(maxDragCoefficient_);
    }
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_, false);
            blades_[i].calcInfluenceEpsilon(maxDragCoefficient_);
        }

        if (hasHub_)
        {
            hub_->setAzimuthIndex(azimuthIndex_, false);
            hub_->calcInfluenceEpsilon(maxDragCoefficient_);
        }
        rotateAD();
    }
    distributeEpsilon();

    labelList globalToLocal(mesh_.nCells(), -1);

    label nActive = 0;
    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->constructInfluenceCellList
        (
            0,
            globalToLocal,
            nActive
        );
    }

    if (hasNacelle_)
    {
        // Add source for nacelle actuator line
        nacelle_->constructInfluenceCellList
        (
            0,
            globalToLocal,
            nActive
        );
    }
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        // Add scalar source term from blades
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_, false);
            blades_[i].constructInfluenceCellList
            (
                azimuthIndex_,
                globalToLocal,
                nActive
            );
        }

        if (hasHub_)
        {
            // Add source for hub actuator line
            hub_->setAzimuthIndex(azimuthIndex_, false);
            hub_->constructInfluenceCellList
            (
                azimuthIndex_,
                globalToLocal,
                nActive
            );
        }
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
    filteredForceField_.setSize(0);
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

void Foam::fv::axialFlowTurbineADSource::setupPositions(bool includeRing)
{
    // A bit more complicated than other turbine types,
    // as tower and nacelle are not rotating
    if (hasTower_)
    {
        tower_->setAzimuthIndex(0, false); // not really needed, remove later
        tower_->findCells(includeRing);
    }

    if (hasNacelle_)
    {
        nacelle_->setAzimuthIndex(0, false); // not really needed, remove later
        nacelle_->findCells(includeRing);
    }
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_, false);
            blades_[i].findCells(includeRing);
        }

        if (hasHub_)
        {
            hub_->setAzimuthIndex(azimuthIndex_, false);
            hub_->findCells(includeRing);
        }
        rotateAD();
    }
}

void Foam::fv::axialFlowTurbineADSource::calculateForces()
{
    // A bit more complicated than other turbine types,
    // as tower and nacelle are not rotating
    if (hasTower_)
    {
        tower_->setAzimuthIndex(0, false); // not really needed, remove later
        tower_->calculateElementForces();
    }

    if (hasNacelle_)
    {
        nacelle_->setAzimuthIndex(0, false); // not really needed, remove later
        nacelle_->calculateElementForces();
    }
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
        {
            forAll(blades_, i)
            {
                blades_[i].setAzimuthIndex(azimuthIndex_, false);
                blades_[i].setCustomTime
                (
                    customTime_[azimuthIndex_],
                    customDeltaT_
                );
            }

            if (hasHub_)
            {
                hub_->setAzimuthIndex(azimuthIndex_, false);
                hub_->setCustomTime
                (
                    customTime_[azimuthIndex_],
                    customDeltaT_
                );
            }

            if (endEffectsActive_ and endEffectsModel_ != "liftingLine")
            {
                // Calculate end effects based on current velocity field
                calcEndEffects();
            }
            
            forAll(blades_, i)
            {
                blades_[i].calculateElementForces();
            }

            if (hasHub_)
            {
                hub_->calculateElementForces();
            }
        }
    }
}

void Foam::fv::axialFlowTurbineADSource::addForce
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

    // tower and nacelle are not rotating,
    // so we only need to calculate the force field once
    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->setAzimuthIndex(0, true);
        tower_->addForceFromChild(forceField, 1.0, compressible);
    }

    if (hasNacelle_)
    {
        // Add source for tower actuator line
        nacelle_->setAzimuthIndex(0, true);
        nacelle_->addForceFromChild(forceField, 1.0, compressible);
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

        // Should not be needed, but for future safety
        forAll(actuatorLines_, i)
        {
            actuatorLines_[i]->setCustomTime
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
        }

        // Add source for blade actuator lines
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_, true);
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

        if (hasHub_)
        {
            // Add source for hub actuator line
            hub_->setAzimuthIndex(azimuthIndex_, true);
            hub_->addForceFromChild
            (
                forceField,
                1.0/divisions_,
                compressible
            );

            force_ += hub_->force();
            moment += hub_->moment(origin_);
        }

        // tower and nacell are both added outside look as they do not move
        if (hasTower_)
        {
            if (includeTowerDrag_)
            {
                force_ += tower_->force();
            }
        }

        if (hasNacelle_)
        {
            if (includeNacelleDrag_)
            {
                force_ += nacelle_->force();
            }
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

void Foam::fv::axialFlowTurbineADSource::addTurbulence
(
    fvMatrix<scalar>& eqn,
    const word fieldName
)
{
    // kField should be the average during one revolution here
    fvMatrix<scalar> kField(eqn.psi(), eqn.dimensions());
    kField *= dimensionedScalar("zero", eqn.dimensions(), 0.0);
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        // Add scalar source term from blades
        forAll(actuatorLines_, i)
        {
            actuatorLines_[i]->setAzimuthIndex(azimuthIndex_, false);
            actuatorLines_[i]->setCustomTime
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            actuatorLines_[i]->addTurbulence(kField, fieldName);
        }
    }
    eqn += (static_cast<scalar>(bladeMultiplier_)/divisions_)*kField;
}

// ************************************************************************* //
