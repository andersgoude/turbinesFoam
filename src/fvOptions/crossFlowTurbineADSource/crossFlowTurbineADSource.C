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
    crossFlowTurbineALSource(name, modelType, dict, mesh),
    firstUse_(true)
{
    read(dict);
    baseCustomTime_ = mesh.time().value();
    rotateAD(true);

    forAll(blades_, i)
    {
        blades_[i].setApplyForce(false);
    }

    if (hasStruts_)
    {
        forAll(struts_, i)
        {
            struts_[i].setApplyForce(false);
        }
    }

    if (hasShaft_)
    {
        shaft_->setApplyForce(false);
    }
    //buildInfluenceCells();
    // reset these after buildInfluenceCells
    baseCustomTime_ = mesh.time().value();
    baseAngleDeg_ = 0;
    angleDeg_[0] = 0;
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
                actuatorLines_[i]->setAzimuthIndex(azimuthIndex_);
                actuatorLines_[i]->elements()[j].calcInfluenceEpsilon();
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
            actuatorLines_[i]->setAzimuthIndex(azimuthIndex_);
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
    firstUse_ = false;
}

void Foam::fv::crossFlowTurbineADSource::setupPositions(bool includeRing)
{
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        forAll(actuatorLines_, i)
        {
            actuatorLines_[i]->setAzimuthIndex(azimuthIndex_);
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
                actuatorLines_[i]->setAzimuthIndex(azimuthIndex_);
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
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        // Zero out force vector and field
        force_ *= 0;

        // Create local moment vector
        vector moment(vector::zero);

        // Add source for blade actuator lines
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_);
            blades_[i].setCustomTime // Not needed in current implementation
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            blades_[i].addForce
            (
                eqn,
                forceField_,
                bladeMultiplier_/divisions_,
                compressible
            );
            force_ += bladeMultiplier_*blades_[i].force();
            bladeMoments_[i] = blades_[i].moment(origin_);
            moment += bladeMultiplier_*bladeMoments_[i];
        }

        if (hasStruts_)
        {
            // Add source for strut actuator lines
            forAll(struts_, i)
            {
                struts_[i].setAzimuthIndex(azimuthIndex_);
                struts_[i].setCustomTime
                (
                    customTime_[azimuthIndex_],
                    customDeltaT_
                );
                struts_[i].addForce
                (
                    eqn,
                    forceField_,
                    bladeMultiplier_/divisions_,
                    compressible
                );
                force_ += bladeMultiplier_*struts_[i].force();
                moment += bladeMultiplier_*struts_[i].moment(origin_);
            }
        }

        if (hasShaft_)
        {
            // Add source for shaft actuator line
            shaft_->setAzimuthIndex(azimuthIndex_);
            shaft_->setCustomTime
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            shaft_->addForce
            (
                eqn,
                forceField_,
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

    // In case forceField isn't the same, add the local field to the global one
    if (&forceField != &forceField_)
    {
        forceField += forceField_;
    }

}

void Foam::fv::crossFlowTurbineADSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    calculateALData(fieldI);
    addForce(eqn, forceField_, 1.0, false);
    
    eqn += forceField_;
}


void Foam::fv::crossFlowTurbineADSource::addSup
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


void Foam::fv::crossFlowTurbineADSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    calculateALData(fieldI);
    // forceField_ should be the average during one revolution here
    fvMatrix<scalar> kField(eqn.psi(), eqn.dimensions());
    kField *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    fvMatrix<scalar> kFieldShaft(eqn.psi(), eqn.dimensions());
    kFieldShaft *=
        dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    for (azimuthIndex_ = 0; azimuthIndex_ < divisions_; azimuthIndex_++)
    {
        // Add scalar source term from blades
        forAll(blades_, i)
        {
            blades_[i].setAzimuthIndex(azimuthIndex_);
            blades_[i].setCustomTime // Not needed in current implementation
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            blades_[i].addSup(kField, fieldI);
        }

        if (hasStruts_)
        {
            // Add source for strut actuator lines
            forAll(struts_, i)
            {
                struts_[i].setAzimuthIndex(azimuthIndex_);
                struts_[i].setCustomTime
                (
                    customTime_[azimuthIndex_],
                    customDeltaT_
                );
                struts_[i].addSup(kField, fieldI);
            }
        }

        if (hasShaft_)
        {
            // Add source for shaft actuator line
            shaft_->setAzimuthIndex(azimuthIndex_);
            shaft_->setCustomTime
            (
                customTime_[azimuthIndex_],
                customDeltaT_
            );
            shaft_->addSup(kFieldShaft, fieldI);
        }
    }
    eqn += (bladeMultiplier_/divisions_)*kField
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
}


void Foam::fv::crossFlowTurbineADSource::updateForceField()
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


bool Foam::fv::crossFlowTurbineADSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        //crossFlowTurbineALSource::read(dict);

        // Get number of divisions
        divisions_ = coeffs_.lookupOrDefault("divisions", 180);
        
        // Get number of divisions
        dynStallLoop_ = coeffs_.lookupOrDefault("dynStallLoop", 1);
        
        // Get blade multiplier
        bladeMultiplier_ = coeffs_.lookupOrDefault("bladeMultiplier", 1.0);

        // Get compact field
        compactField_ = coeffs_.lookupOrDefault("compactField", true);

        // Get compact field
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
