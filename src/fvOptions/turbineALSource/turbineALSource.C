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

#include "turbineALSource.H"
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
    defineTypeNameAndDebug(turbineALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        turbineALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::turbineALSource::rotateVector
(
    vector& vectorToRotate,
    vector rotationPoint,
    vector axis,
    scalar radians
)
{
    // Declare and define the rotation matrix (from SOWFA)
    tensor RM;
    scalar angle = radians;
    RM.xx() = Foam::sqr(axis.x())
            + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle);
    RM.xz() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y())
            + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z())
            + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated.
    vectorToRotate -= rotationPoint;

    // Perform the rotation.
    vectorToRotate = RM & vectorToRotate;

    // Return the rotated point to its new location relative to the rotation
    // point
    vectorToRotate += rotationPoint;
}


scalar Foam::fv::turbineALSource::calculateCone
(
    const List<List<scalar>>& elementData,
    const label j
)
{
    // Use provided cone value if available
    if (elementData[j].size() > 6)
    {
        return degToRad(elementData[j][6]);
    }

    scalar dr;
    scalar dx;

    // Forward difference at first point
    if (j == 0)
    {
        dr = elementData[j+1][1] - elementData[j][1];
        dx = elementData[j+1][0] - elementData[j][0];
    }
    // Backward difference at last point
    else if (j == elementData.size() - 1)
    {
        dr = elementData[j][1] - elementData[j-1][1];
        dx = elementData[j][0] - elementData[j-1][0];
    }
    // Central difference for interior points
    else
    {
        dr = elementData[j+1][1] - elementData[j-1][1];
        dx = elementData[j+1][0] - elementData[j-1][0];
    }

    // Protect against division by zero
    if (Foam::mag(dx) < SMALL)
    {
        if (dr > 0)
        {
            return constant::mathematical::pi/2.0;
        }
        else if (dr < 0)
        {
            return -constant::mathematical::pi/2.0;
        }
        else
        {
            return 0.0;
        }
    }

    return Foam::atan(dr/dx);
}


void Foam::fv::turbineALSource::createCoordinateSystem()
{
    // Should be unique for each type of turbine
}


void Foam::fv::turbineALSource::createBlades()
{
    // Should be unique for each type of turbine
}


void Foam::fv::turbineALSource::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = time_.path()/"../postProcessing/turbines"
            / time_.timeName();
    }
    else
    {
        dir = time_.path()/"postProcessing/turbines"
            / time_.timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_<< "time,angle_deg,tsr,cp,cd,ct";

    forAll(blades_, i)
    {
        *outputFile_<< ",cd_" << bladeNames_[i];
        *outputFile_<< ",ct_" << bladeNames_[i];
    }

    *outputFile_<< endl;
}


void Foam::fv::turbineALSource::updateTSROmega()
{
    // Update tip speed ratio and omega
    scalar theta = degToRad(angleDeg_[azimuthIndex_]);
    tipSpeedRatio_ = meanTSR_ + tsrAmplitude_*cos(nBlades_*(theta - tsrPhase_));
    omega_ = tipSpeedRatio_*mag(freeStreamVelocity_)/rotorRadius_;
}


void Foam::fv::turbineALSource::rotate()
{
    scalar deltaT = time_.deltaT().value();
    scalar radians = omega_*deltaT;
    rotate(radians);
    baseAngleDeg_ += radToDeg(radians);
    angleDeg_[azimuthIndex_] = baseAngleDeg_;
    lastRotationTime_ = time_.value();
    updateTSROmega();
}


void Foam::fv::turbineALSource::rotate(scalar radians)
{
    // Should be defined for each turbine type
}


void Foam::fv::turbineALSource::printPerf()
{
    Info<< "Azimuthal angle (degrees) of " << name_ << ": "
        << angleDeg_[azimuthIndex_] << endl;
    Info<< "Tip speed ratio of " << name_ << ": " << tipSpeedRatio_ << endl;
    Info<< "Power coefficient from " << name_ << ": " << powerCoefficient_
        << endl;
    Info<< "Rotor drag coefficient from " << name_ << ": " << dragCoefficient_
        << endl << endl;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::turbineALSource::turbineALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    actuatorModelBase(name, modelType, dict, mesh),
    time_(mesh.time()),
    lastRotationTime_(time_.value()),
    rhoRef_(1.0),
    omega_(0.0),
    baseAngleDeg_(0.0),
    angleDeg_(1, 0.0),
    baseCustomTime_(0.0),
    customTime_(1, 0.0),
    customDeltaT_(0.0),
    azimuthIndex_(0),
    nBlades_(0),
    freeStreamVelocity_(vector::zero),
    frontalArea_(0.0),
    powerCoefficient_(0.0),
    dragCoefficient_(0.0),
    torqueCoefficient_(0.0),
    cylMin_(GREAT),
    cylMax_(-GREAT),
    cylRadius_(0.0),
    chordMax_(0.0),
    epsilonLiftMax_(0.0),
    dragMax_(0.0),
    meshFactorMax_(0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::turbineALSource::~turbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::turbineALSource::createForceFieldForChildren
(
    const bool compressible
)
{
    // Create for myself (only if writeForceField_ is true)
    createForceField(false, compressible);
    forAll(actuatorLines_, i)
    {
        actuatorLines_[i]->createForceFieldForChildren(compressible);
    }
}

void Foam::fv::turbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


void Foam::fv::turbineALSource::writePerf()
{
    *outputFile_<< time_.value() << "," << angleDeg_[azimuthIndex_] << ","
                << tipSpeedRatio_ << "," << powerCoefficient_ << ","
                << dragCoefficient_ << "," << torqueCoefficient_;

    // Write power, drag, and torque coefficients for each blade
    forAll(blades_, i)
    {
        // Write drag (thrust) coefficient contribution from blade
        scalar bladeCd = blades_[i].force() & freeStreamDirection_
            / (0.5*frontalArea_*magSqr(freeStreamVelocity_));
        *outputFile_<< "," << bladeCd;
        // Write torque coefficient contribution from blade
        scalar bladeTorque = bladeMoments_[i] & axis_;
        scalar bladeCt = bladeTorque
            / (0.5*frontalArea_*rotorRadius_* magSqr(freeStreamVelocity_));
        *outputFile_<< "," << bladeCt;
    }

    *outputFile_<< endl;
}


void Foam::fv::turbineALSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::turbineALSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Read coordinate system/geometry invariant properties
        coeffs_.lookup("origin") >> origin_;
        coeffs_.lookup("axis") >> axis_;
        axis_ /= mag(axis_);
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        coeffs_.lookup("tipSpeedRatio") >> meanTSR_;
        coeffs_.lookup("rotorRadius") >> rotorRadius_;
        tsrAmplitude_ = coeffs_.lookupOrDefault("tsrAmplitude", 0.0);
        tsrPhase_ = coeffs_.lookupOrDefault("tsrPhase", 0.0);

        // Get blade information
        bladesDict_ = coeffs_.subDict("blades");
        nBlades_ = bladesDict_.keys().size();
        bladeNames_ = bladesDict_.toc();
        bladeMoments_.setSize(nBlades_);

        // Set tip speed ratio and omega
        updateTSROmega();

        // Get dynamic stall subdict
        dynamicStallDict_ = coeffs_.subOrEmptyDict("dynamicStall");

        // Get profiles information
        profileData_ = coeffs_.subDict("profileData");

        // For automatically determine cells that can interaction with the
        // element, use this drag coefficient to determin drag based sphere
        // radius drag should usually not be more than 4, but it is best that
        // the user can put a higher value if needed
        maxDragCoefficient_ = coeffs_.lookupOrDefault
        (
            "maxDragCoefficient",
            4.0
        );


        // actuator disc variables below

        // Get number of divisions
        divisions_ = coeffs_.lookupOrDefault("divisions", 180);
        
        // Get number of divisions
        dynStallLoop_ = coeffs_.lookupOrDefault("dynStallLoop", 1);
        
        // Get blade multiplier
        bladeMultiplier_ = coeffs_.lookupOrDefault("bladeMultiplier", 1);

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

        if (divisions_ % bladeMultiplier_ != 0 && nBlades_ == 1)
        {
            FatalErrorIn("void turbineALSource::read()")
                    << "divisions must be a multiple of bladeMultiplier"
                    << "current values are: divisions = " << divisions_
                    << " bladeMultiplier = " << bladeMultiplier_
                    << abort(FatalError);
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
        Info << "bladeMultiplier = " << bladeMultiplier_ << " divisions = " << divisions_ << endl;

        return true;
    }
    else
    {
        Info << "turbineALSource.read() failed " << endl;
        return false;
    }
}


scalar Foam::fv::turbineALSource::powerCoefficient() const
{
    // Return power coefficient
    return powerCoefficient_;
}

scalar Foam::fv::turbineALSource::dragCoefficient() const
{
    // Return drag coefficient
    return dragCoefficient_;
}

scalar Foam::fv::turbineALSource::meanPowerCoefficient() const
{
    // Return power coefficient
    return meanPowerCoefficient_;
}

scalar Foam::fv::turbineALSource::meanDragCoefficient() const
{
    // Return drag coefficient
    return meanDragCoefficient_;
}

void Foam::fv::turbineALSource::determineBoundingCylinder
(
    PtrList<actuatorLineSource> &lines
)
{
    scalar radSqr = cylRadius_*cylRadius_;

    forAll(lines, i)
    {
        forAll(lines[i].elements(), j)
        {
            const point& p = lines[i].elements()[j].position();
            const scalar chordLength = lines[i].elements()[j].chordLength();
            const scalar chordValue =
                chordLength * lines[i].elements()[j].chordFactor();
            const scalar dragValue =
                chordLength * lines[i].elements()[j].dragFactor();
            const scalar meshValue = lines[i].elements()[j].meshFactor();

            // vector from centre to the point
            const vector d = p - origin_;

            // axial coordinate relative to centre
            const scalar s = d & axis_;

            cylMin_ = Foam::min(cylMin_, s);
            cylMax_ = Foam::max(cylMax_, s);
            chordMax_ = Foam::max(chordMax_, chordLength);
            epsilonLiftMax_ = Foam::max(epsilonLiftMax_, chordValue);
            dragMax_ = Foam::max(dragMax_, dragValue);
            meshFactorMax_ = Foam::max(meshFactorMax_, meshValue);

            // radial distance squared from the axis
            const vector radial = d - s*axis_;
            radSqr = Foam::max(radSqr, magSqr(radial));
        }
    }

    // Final radius
    const scalar radius = Foam::sqrt(radSqr);
    cylRadius_ = Foam::max(cylRadius_, radius);
}

scalar Foam::fv::turbineALSource::maxCellVolumeInCylinder()
{
    scalar maxVol  = -GREAT;
    scalar addition = 0.05*cylRadius_; //make cylinder slightly larger

    // if we do not find any cell, increase addition and try again
    while (maxVol < 0 && addition < cylRadius_)
    {
        const point  p1 = origin_ + (cylMin_ - addition)*axis_; // lower end
        const point  p2 = origin_ + (cylMax_ + addition)*axis_; // upper end
        const vector axisVec = p2 - p1;
        const scalar magAxis2 = magSqr(axisVec);
        const scalar invMaxAxis2 = 1.0/magAxis2;
        const scalar rad2     = sqr(cylRadius_ + addition);

        const scalarField& V = mesh_.V();
        const vectorField& C = mesh_.C();

        forAll(C, cellI)
        {
            const vector d = C[cellI] - p1;
            const scalar magD = d & axisVec; // projection onto axis

            // axial bounds
            if (magD > 0 && magD < magAxis2)
            {
                // radial distance squared
                const scalar d2 = magSqr(d) - sqr(magD)*invMaxAxis2;

                if (d2 < rad2) // inside the cylinder
                {
                    if (V[cellI] > maxVol)
                    {
                        maxVol = V[cellI];
                    }
                }
            }
        }

        reduce(maxVol, maxOp<scalar>());
        addition *= 2;
    }
    return maxVol;
}

void Foam::fv::turbineALSource::addPointsInCylinder
(
    scalar maxVol,
    labelList& globalToLocal,
    label& nActive
)
{
    // Epsilon based on drag/momentum thickness
    scalar epsilonDrag = maxDragCoefficient_*dragMax_/2.0;

    scalar epsilonMesh = 2.0*Foam::cbrt(maxVol)*meshFactorMax_;

    // Threshold is based on lift or drag, whichever is larger
    scalar epsilonThreshold = Foam::max(epsilonLiftMax_, epsilonDrag);
    scalar epsilon = Foam::max(epsilonThreshold, epsilonMesh);

    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));

    // Force can be applied within this sphere radius
    scalar sphereRadius = chordMax_ + projectionRadius;

    if (Pstream::master())
    {
        Info<< "Drag based criterion: " << epsilonDrag << endl
            << "Chord based criterion: " << epsilonLiftMax_ << endl
            << "Mesh based criterion " << epsilonMesh << endl
            << "Using interaction sphere of " << projectionRadius
            << " to determine maximum distance for element interaction" << endl;
    }

    // for safety, add 5 % extra on the sphere Radius
    const point  p1 = origin_ + (cylMin_ - 1.05*sphereRadius)*axis_; // lower
    const point  p2 = origin_ + (cylMax_ + 1.05*sphereRadius)*axis_; // upper
    const vector axisVec = p2 - p1;
    const scalar magAxis2 = magSqr(axisVec);
    scalar invMaxAxis2 = 0.0;
    if (magAxis2 > 0)
    {
        invMaxAxis2 = 1.0/magAxis2;
    }
    else
    {
        FatalErrorIn("void turbineALSource::addPointsInCylinder()")
                << "Invalid cylinder size"
                << abort(FatalError);
    }
    const scalar rad2 = sqr(cylRadius_ + 1.05*sphereRadius);

    const vectorField& C = mesh_.C();

    forAll(C, cellI)
    {
        const vector d = C[cellI] - p1;
        const scalar magD = d & axisVec;               // projection onto axis

        // axial bounds
        if (magD > 0 && magD < magAxis2)
        {
            // radial distance squared
            const scalar d2 = magSqr(d) - sqr(magD)*invMaxAxis2;

            if (d2 < rad2)                             // inside the cylinder
            {
                globalToLocal[cellI] = nActive++;
            }
        }
    }
}

// Actuator line with compact fields
void Foam::fv::turbineALSource::updateForceFieldAL
(
    volVectorField &forceField
)
{
    if (activeForceField_.size() > 0)
    {
        forAll(localToGlobal_, forceIndex)
        {
            forceField[localToGlobal_[forceIndex]] +=
                activeForceField_[forceIndex];
        }
    }
}

// Actuator disc specific with possibility for relaxation
void Foam::fv::turbineALSource::updateForceFieldAD
(
    volVectorField &forceField
)
{
    if (relaxForceField_)
    {
        if (filteredForceField_.size() == 0)
        {
            filteredForceField_ = activeForceField_;
            forAll(localToGlobal_, forceIndex)
            {
                forceField[localToGlobal_[forceIndex]] +=
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
                    forceField[localToGlobal_[forceIndex]] +=
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
                    forceField[localToGlobal_[forceIndex]] +=
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
                forceField[localToGlobal_[forceIndex]] +=
                    activeForceField_[forceIndex];
            }
        }
    }
}
// ************************************************************************* //
