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
    customTime_ = mesh.time().value();
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
    customTime_ = mesh.time().value();
    angleDeg_ = 0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineADSource::~crossFlowTurbineADSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::crossFlowTurbineADSource::rotateAD(bool updateOnly)
{
    scalar radians = 2*mathematical::pi/divisions_;
    scalar deltaT = radians/omega_;

    //updateOnly is intended for first step only to set custom speed
    if (updateOnly == false)
    {
        customTime_ += deltaT;
        rotate(radians);
        angleDeg_ += radToDeg(radians);
        //lastRotationTime_ = time_.value();
    }
    updateTSROmega();
    
    // Info << "rotateAD called for time " << time_.value()
    // << " custom time: " << customTime_ << endl;
    forAll(blades_, i)
    {
        blades_[i].setCustomTime(customTime_, deltaT);
    }

    if (hasStruts_)
    {
        forAll(struts_, i)
        {
            struts_[i].setCustomTime(customTime_, deltaT);
        }
    }

    if (hasShaft_)
    {
        shaft_->setCustomTime(customTime_, deltaT);
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


void Foam::fv::crossFlowTurbineADSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Generate UInterp object to be used for all velocity interpolations
    const volVectorField& Uin(eqn.psi());
    interpolationCellPoint<vector> UInterp(Uin);

    if (firstUse_)
    {
        buildInfluenceCells();
        angleDeg_ = 0;
        firstUse_ = false;
    }
    // code can run extra revolutions to make dynamic stall converge
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        // forceField_ should be the average during one revolution here
        //forceField_ *= dimensionedScalar("zero",forceField_.dimensions(),0.0)
        if (activeForceField_.size() > 0)
        {
            activeForceField_ = vector::zero;
        }
        else
        {
            // forceField_ should be the average during one revolution here
            forceField_.primitiveFieldRef() = vector::zero;
            forceField_.correctBoundaryConditions();

            // Check dimensions of force field and correct if necessary
            if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
            {
                forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
            }
        }
        for (int innerStep = 0; innerStep < divisions_; innerStep++)
        {
            // Zero out force vector and field
            force_ *= 0;

            // Create local moment vector
            vector moment(vector::zero);

            // Add source for blade actuator lines
            forAll(blades_, i)
            {
                blades_[i].setAzimuthIndex(innerStep);
                blades_[i].addForce
                (
                    eqn,
                    UInterp,
                    forceField_,
                    fieldI,
                    bladeMultiplier_/divisions_
                );
                //forceField_ +=
                //    (bladeMultiplier_/divisions_)*blades_[i].forceField();
                //Info<< "Added blade" << endl;
                force_ += bladeMultiplier_*blades_[i].force();
                bladeMoments_[i] = blades_[i].moment(origin_);
                moment += bladeMultiplier_*bladeMoments_[i];
            }

            if (hasStruts_)
            {
                // Add source for strut actuator lines
                forAll(struts_, i)
                {
                    struts_[i].setAzimuthIndex(innerStep);
                    struts_[i].addForce
                    (
                        eqn,
                        UInterp,
                        forceField_,
                        fieldI,
                        bladeMultiplier_/divisions_
                    );
                    //forceField_ +=
                    //    (bladeMultiplier_/divisions_)*struts_[i].forceField();
                    force_ += bladeMultiplier_*struts_[i].force();
                    moment += bladeMultiplier_*struts_[i].moment(origin_);
                }
            }

            if (hasShaft_)
            {
                // Add source for shaft actuator line
                shaft_->setAzimuthIndex(innerStep);
                shaft_->addForce
                (
                    eqn,
                    UInterp,
                    forceField_,
                    fieldI,
                    1.0/divisions_
                );
                //forceField_ += (1.0/divisions_)*shaft_->forceField();
                force_ += shaft_->force();
                moment += shaft_->moment(origin_);
            }

            // only needed to to do calculations for the actual loop,
            // not the dummy loops that are used to make dynamic stall converge
            if (currentLoop == dynStallLoop_ - 1)
            {
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
            rotateAD();
        }
    }

    // When using compressed fields, restore them to the original forceField
    if (activeForceField_.size() > 0)
    {
        forAll(localToGlobal_, forceIndex)
        {
            forceField_[localToGlobal_[forceIndex]] =
                activeForceField_[forceIndex];
        }
    }
    eqn += forceField_;
}


void Foam::fv::crossFlowTurbineADSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Generate UInterp object to be used for all velocity interpolations
    const volVectorField& Uin(eqn.psi());
    interpolationCellPoint<vector> UInterp(Uin);

    if (firstUse_)
    {
        buildInfluenceCells();
        angleDeg_ = 0;
        firstUse_ = false;
    }
    // code can run extra revolutions to make dynamic stall converge
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        if (activeForceField_.size() > 0)
        {
            activeForceField_ = vector::zero;
        }
        else
        {
            // forceField_ should be the average during one revolution here
            forceField_.primitiveFieldRef() = vector::zero;
            forceField_.correctBoundaryConditions();

            // Check dimensions of force field and correct if necessary
            if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
            {
                forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
            }
        }
        for (int innerStep = 0; innerStep < divisions_; innerStep++)
        {
            // Zero out force vector and field
            force_ *= 0;

            // Create local moment vector
            vector moment(vector::zero);

            // Add source for blade actuator lines
            forAll(blades_, i)
            {
                blades_[i].setAzimuthIndex(innerStep);
                blades_[i].addForce
                (
                    rho,
                    eqn,
                    UInterp,
                    forceField_,
                    fieldI,
                    bladeMultiplier_/divisions_
                );
                //forceField_ +=
                //    (bladeMultiplier_/divisions_)*blades_[i].forceField();
                force_ += bladeMultiplier_*blades_[i].force();
                bladeMoments_[i] = blades_[i].moment(origin_);
                moment += bladeMultiplier_*bladeMoments_[i];
            }

            if (hasStruts_)
            {
                // Add source for strut actuator lines
                forAll(struts_, i)
                {
                    struts_[i].setAzimuthIndex(innerStep);
                    struts_[i].addForce
                    (
                        rho,
                        eqn,
                        UInterp,
                        forceField_,
                        fieldI,
                        bladeMultiplier_/divisions_
                    );
                    //forceField_ +=
                    //  (bladeMultiplier_/divisions_)*struts_[i].forceField();
                    force_ += bladeMultiplier_*struts_[i].force();
                    moment += bladeMultiplier_*struts_[i].moment(origin_);
                }
            }

            if (hasShaft_)
            {
                // Add source for shaft actuator line
                shaft_->setAzimuthIndex(innerStep);
                shaft_->addForce
                (
                    rho,
                    eqn,
                    UInterp,
                    forceField_,
                    fieldI,
                    1.0/divisions_
                );
                //forceField_ += (1.0/divisions_)*shaft_->forceField();
                force_ += shaft_->force();
                moment += shaft_->moment(origin_);
            }
            
            if (currentLoop == dynStallLoop_ - 1)
            {
                // Torque is the projection of the moment from
                // all blades on the axis
                torque_ = moment & axis_;

                scalar rhoRef;
                coeffs_.lookup("rhoRef") >> rhoRef;
                torqueCoefficient_ =
                    torque_/(0.5*rhoRef*frontalArea_*rotorRadius_
                    * magSqr(freeStreamVelocity_));
                powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
                dragCoefficient_ =
                    force_ & freeStreamDirection_
                    / (0.5*rhoRef*frontalArea_*magSqr(freeStreamVelocity_));

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
            rotateAD();
        }
    }

    // When using compressed fields, restore them to the original forceField
    if (activeForceField_.size() > 0)
    {
        forAll(localToGlobal_, forceIndex)
        {
            forceField_[localToGlobal_[forceIndex]] =
                activeForceField_[forceIndex];
        }
    }

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
    if (firstUse_)
    {
        buildInfluenceCells();
        angleDeg_ = 0;
        firstUse_ = false;
    }
    // code can run extra revolutions to make dynamic stall converge
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        // forceField_ should be the average during one revolution here
        fvMatrix<scalar> kField(eqn.psi(), eqn.dimensions());
        kField *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
        fvMatrix<scalar> kFieldShaft(eqn.psi(), eqn.dimensions());
        kFieldShaft *=
            dimensionedScalar("zero", forceField_.dimensions(), 0.0);
        for (int innerStep = 0; innerStep < divisions_; innerStep++)
        {
            // Add scalar source term from blades
            forAll(blades_, i)
            {
                blades_[i].setAzimuthIndex(innerStep);
                blades_[i].addSup(kField, fieldI);
            }

            if (hasStruts_)
            {
                // Add source for strut actuator lines
                forAll(struts_, i)
                {
                    struts_[i].setAzimuthIndex(innerStep);
                    struts_[i].addSup(kField, fieldI);
                }
            }

            if (hasShaft_)
            {
                // Add source for shaft actuator line
                shaft_->setAzimuthIndex(innerStep);
                shaft_->addSup(kFieldShaft, fieldI);
            }
            rotateAD();
        }
        eqn += (bladeMultiplier_/divisions_)*kField
               + (1.0/divisions_)*kFieldShaft;
    }
}

void Foam::fv::crossFlowTurbineADSource::buildInfluenceCells()
{
    forAll(blades_, i)
    {
        blades_[i].allocateInfluenceCells(divisions_, cacheInteractions_);
    }

    if (hasStruts_)
    {
        forAll(struts_, i)
        {
            struts_[i].allocateInfluenceCells(divisions_, cacheInteractions_);
        }
    }

    if (hasShaft_)
    {
        // Add source for tower actuator line
        shaft_->allocateInfluenceCells(divisions_, cacheInteractions_);
    }

    //- Allocations are still needed for velocity lookup, but if compactField
    // is false, we do not need influenceCells
    if (compactField_ == false)
    {
        return;
    }

    labelList globalToLocal(mesh_.nCells(), -1);

    label nActive = 0;
    for (label innerStep = 0; innerStep < divisions_; innerStep++)
    {
        // Add scalar source term from blades
        forAll(blades_, i)
        {
            blades_[i].constructInfluenceCellList
            (
                innerStep,
                globalToLocal,
                nActive
            );
        }

        if (hasStruts_)
        {
            forAll(struts_, i)
            {
                struts_[i].constructInfluenceCellList
                (
                    innerStep,
                    globalToLocal,
                    nActive
                );
            }
        }

        if (hasShaft_)
        {
            // Add source for tower actuator line
            shaft_->constructInfluenceCellList
            (
                innerStep,
                globalToLocal,
                nActive
            );
        }
        rotateAD();
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

    forAll(blades_, i)
    {
        blades_[i].setCompactFields(activePositions_, activeForceField_);
    }
    if (hasStruts_)
    {
        forAll(struts_, i)
        {
            struts_[i].setCompactFields
            (
                activePositions_,
                activeForceField_
            );
        }
    }

    if (hasShaft_)
    {
        // Add source for tower actuator line
        shaft_->setCompactFields
        (
            activePositions_,
            activeForceField_
        );
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
        
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
