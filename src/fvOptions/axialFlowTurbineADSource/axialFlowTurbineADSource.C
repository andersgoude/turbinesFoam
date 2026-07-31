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
    axialFlowTurbineALSource(name, modelType, dict, mesh),
    firstUse_(true)
{
    read(dict);
    customTime_ = mesh.time().value();
    rotateAD(true);

    // override the nBlades value for the end effects calculation
    // if bladeMultiplier is used
    effectiveNBlades_ = bladeMultiplier_*nBlades_;
    forAll(blades_, i)
    {
        blades_[i].setApplyForce(false);
    }

    if (hasHub_)
    {
        hub_->setApplyForce(false);
    }

    if (hasTower_)
    {
        tower_->setApplyForce(false);
    }

    if (hasNacelle_)
    {
        nacelle_->setApplyForce(false);
    }
    //buildInfluenceCells();
    // reset these after buildInfluenceCells
    customTime_ = mesh.time().value();
    baseAngleDeg_ = 0;
    angleDeg_[0] = 0;
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
        tower_->setAzimuthIndex(0); // not really needed, remove later
        tower_->calcInfluenceEpsilon();
    }

    if (hasNacelle_)
    {
        nacelle_->setAzimuthIndex(0); // not really needed, remove later
        nacelle_->calcInfluenceEpsilon();
    }
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_);
            blades_[i].calcInfluenceEpsilon();
        }

        if (hasHub_)
        {
            hub_->setAzimuthIndex(azimuthIndex_);
            hub_->calcInfluenceEpsilon();
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
            blades_[i].setAzimuthIndex(azimuthIndex_);
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
            hub_->setAzimuthIndex(azimuthIndex_);
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
    firstUse_ = false;
}

void Foam::fv::axialFlowTurbineADSource::setupPositions(bool includeRing)
{
    // A bit more complicated than other turbine types,
    // as tower and nacelle are not rotating
    if (hasTower_)
    {
        tower_->setAzimuthIndex(0); // not really needed, remove later
        tower_->findCells(includeRing);
    }

    if (hasNacelle_)
    {
        nacelle_->setAzimuthIndex(0); // not really needed, remove later
        nacelle_->findCells(includeRing);
    }
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_);
            blades_[i].findCells(includeRing);
        }

        if (hasHub_)
        {
            hub_->setAzimuthIndex(azimuthIndex_);
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
        tower_->setAzimuthIndex(0); // not really needed, remove later
        tower_->calculateElementForces();
    }

    if (hasNacelle_)
    {
        nacelle_->setAzimuthIndex(0); // not really needed, remove later
        nacelle_->calculateElementForces();
    }
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
        {
            forAll(blades_, i)
            {
                blades_[i].setAzimuthIndex(azimuthIndex_);
                blades_[i].setCustomTime
                (
                    customTime_[azimuthIndex_],
                    customDeltaT_
                );
            }

            if (hasHub_)
            {
                hub_->setAzimuthIndex(azimuthIndex_);
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
    fvMatrix<vector> &eqn,
    volVectorField &forceField,
    scalar scale,
    bool compressible
)
{
    // forceField_ should be the average during one revolution here
    if (compactField_)
    {
        if (activeForceField_.size() > 0)
        {
            activeForceField_ = vector::zero;
        }
    }
    else
    {
        forceField_.primitiveFieldRef() = vector::zero;
    }
    forceField_.correctBoundaryConditions();

    // Check dimensions of force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }
    
    // tower and nacelle are not rotating,
    // so we only need to calculate the force field once
    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->setAzimuthIndex(0);
        tower_->addForce(eqn, forceField_, 1.0, compressible);
    }

    if (hasNacelle_)
    {
        // Add source for tower actuator line
        nacelle_->setAzimuthIndex(0);
        nacelle_->addForce(eqn, forceField_, 1.0, compressible);
    }
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
            blades_[i].setAzimuthIndex(azimuthIndex_);
            blades_[i].addForce
            (
                eqn,
                forceField_,
                bladeMultiplier_/divisions_,
                compressible
            );
            //forceField_ +=
            //    (bladeMultiplier_/divisions_)*blades_[i].forceField();
            //Info<< "Added blade" << endl;
            force_ += bladeMultiplier_*blades_[i].force();
            bladeMoments_[i] = blades_[i].moment(origin_);
            moment += bladeMultiplier_*bladeMoments_[i];
        }

        if (hasHub_)
        {
            // Add source for hub actuator line
            hub_->setAzimuthIndex(azimuthIndex_);
            hub_->addForce
            (
                eqn,
                forceField_,
                1.0/divisions_,
                compressible
            );
        //    forceField_ += (1.0/divisions_)*hub_->forceField();
            force_ += hub_->force();
            moment += hub_->moment(origin_);
        }

        if (hasTower_)
        {
            // Add source for tower actuator line
            //tower_->addSup(eqn, fieldI);
            //forceField_ += (1.0/divisions_)*tower_->forceField();
            if (includeTowerDrag_)
            {
                force_ += tower_->force();
            }
        }

        if (hasNacelle_)
        {
            // Add source for tower actuator line
            //nacelle_->addSup(eqn, fieldI);
            //forceField_ += (1.0/divisions_)*nacelle_->forceField();
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


        // Print performance to terminal
        printPerf();

        // Write performance data
        // Note this will write multiples if there are
        // multiple PIMPLE loops
        if (Pstream::master())
        {
            writePerf();
        }
    }

    // When using compressed fields, restore them to the original forceField
    updateForceField();
}

void Foam::fv::axialFlowTurbineADSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    calculateALData(fieldI);
    
    addForce(eqn, forceField_, 1.0, false);

    eqn += forceField_;
}


void Foam::fv::axialFlowTurbineADSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    calculateALData(fieldI);
    
    addForce(eqn, forceField_, 1.0, true);

    // multiply with local density
    forceField_ *= rho;

    eqn += forceField_;
}


void Foam::fv::axialFlowTurbineADSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    calculateALData(fieldI);

    // forceField_ should be the average during one revolution here
    fvMatrix<scalar> kField(eqn.psi(), eqn.dimensions());
    kField *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        if (endEffectsActive_ and endEffectsModel_ != "liftingLine")
        {
            // Calculate end effects based on current velocity field
            calcEndEffects();
        }
        
        // Add scalar source term from blades
        forAll(actuatorLines_, i)
        {
            actuatorLines_[i]->setAzimuthIndex(azimuthIndex_);
            actuatorLines_[i]->setCustomTime
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            actuatorLines_[i]->addSup(kField, fieldI);
        }
    }
    eqn += (bladeMultiplier_/divisions_)*kField;
}


void Foam::fv::axialFlowTurbineADSource::updateForceField()
{
    if (relaxForceField_)
    {
        if (filteredForceField_.size() == 0)
        {
            filteredForceField_ = activeForceField_;
            forAll(localToGlobal_, forceIndex)
            {
                forceField_[localToGlobal_[forceIndex]] =
                    filteredForceField_[forceIndex];
            }
        }
        else
        {
            if (mesh_.time().value() < relaxStartTime_)
            {
                filteredForceField_ = activeForceField_;
                forAll(localToGlobal_, forceIndex)
                {
                    forceField_[localToGlobal_[forceIndex]] =
                        filteredForceField_[forceIndex];
                }
            }
            else
            {
                scalar alpha = 1.0/(1.0 + relaxValue_);
                forAll(localToGlobal_, forceIndex)
                {
                    filteredForceField_[forceIndex] =
                        (1 - alpha)*filteredForceField_[forceIndex]
                        + alpha*activeForceField_[forceIndex];
                    forceField_[localToGlobal_[forceIndex]] =
                        filteredForceField_[forceIndex];
                }
                if (relaxGrowthValue_ > 0)
                {
                    if (relaxGrowthThreshold_ > 0)
                    {
                        // Difference field
                        vectorField diff =
                            activeForceField_ - filteredForceField_;

                        // L2 norm of the difference
                        scalar diffL2 = Foam::sqrt(Foam::sum(magSqr(diff)));

                        // L2 norm of the reference field
                        scalar refL2 =
                            Foam::sqrt(Foam::sum(magSqr(filteredForceField_)));

                        scalar relativeChange = diffL2 / (refL2 + SMALL);
                        if (relativeChange > relaxGrowthThreshold_)
                        {
                            relaxValue_ += relaxGrowthValue_;
                        }
                    }
                    else
                    {
                        relaxValue_ += relaxGrowthValue_;
                    }
                }
                if (relaxValue_ < relaxMaxValue_)
                {
                    relaxValue_ = relaxMaxValue_;
                }
            }
        }
    }
    else
    {
        if (activeForceField_.size() > 0)
        {
            forAll(localToGlobal_, forceIndex)
            {
                forceField_[localToGlobal_[forceIndex]] =
                    activeForceField_[forceIndex];
            }
        }
    }
}

bool Foam::fv::axialFlowTurbineADSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        //crossFlowTurbineALSource::read(dict);

        // Get number of divisions
        divisions_ = coeffs_.lookupOrDefault("divisions", 180);
        
        // Get number of divisions
        dynStallLoop_ = coeffs_.lookupOrDefault("dynStallLoop", 1);
        
        // Get blade multiplier
        bladeMultiplier_ = coeffs_.lookupOrDefault("bladeMultiplier", 1);

        // Get compact field
        compactField_ = coeffs_.lookupOrDefault("compactField", true);

        // Get cache interactions field
        cacheInteractions_ = coeffs_.lookupOrDefault("cacheInteractions", true);

        // For simplicity, ensure that cacheInteractions cannot be true when
        // compactField is false, to avoid having to implement this path
        // as the interctions cache is the memory consuming part
        if (compactField_ == false)
        {
            cacheInteractions_ = false;
        }

        // Get if we should apply relaxation to the force field
        relaxForceField_ = coeffs_.lookupOrDefault("relaxForceField", false);
        relaxStartTime_ = coeffs_.lookupOrDefault("relaxStartTime", 0);
        relaxValue_ = coeffs_.lookupOrDefault("relaxStartValue", 0.3);
        relaxGrowthValue_ = coeffs_.lookupOrDefault("relaxGrowthValue", 0.02);
        relaxGrowthThreshold_ =
            coeffs_.lookupOrDefault("relaxGrowthThreshold", 0.0);
        relaxMaxValue_ = coeffs_.lookupOrDefault("relaxMaxValue", 0.02);

        if (debug)
        {
            Info << "relaxForceField_ " << relaxForceField_
                 << " relaxStartTime_ " << relaxStartTime_
                 << " relaxValue_ " << relaxValue_
                 << " relaxGrowthValue_ " << relaxGrowthValue_
                 << " relaxGrowthThreshold_ " << relaxGrowthThreshold_
                 << " relaxMaxValue_ " << relaxMaxValue_
                 << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
