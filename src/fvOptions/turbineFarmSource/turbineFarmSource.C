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

#include "turbineFarmSource.H"
#include "axialFlowTurbineALSource.H"
#include "axialFlowTurbineADSource.H"
#include "crossFlowTurbineALSource.H"
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
    defineTypeNameAndDebug(turbineFarmSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        turbineFarmSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::turbineFarmSource::createAxialFlowTurbines()
{
    // Create the turbines
    axialFlowTurbines_.setSize(nAxialFlowTurbines_);
    for (int i = 0; i < nAxialFlowTurbines_; i++)
    {
        word turbineName = axialFlowTurbineNames_[i];
        // Create dictionary items for this blade
        dictionary turbineSubDict = axialFlowTurbinesDict_.subDict(turbineName);
        word modelType
        (
            turbineSubDict.lookupOrDefault<word>("type", "actuatorLine")
        );
        turbineSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));

        // copy farm data if not specified
        word selectionMode = turbineSubDict.lookupOrAddDefault
        (
            "selectionMode",
            coeffs_.getOrDefault<word>("selectionMode", "all")
        );

        if (selectionMode != "all")
        {
            turbineSubDict.lookupOrAddDefault
            (
                "cellSet",
                coeffs_.lookup("cellSet")
            );
        }

        // Do not write force from individual turbines unless specified
        turbineSubDict.lookupOrAddDefault("writeForceField", false);

        // set active as default
        turbineSubDict.lookupOrAddDefault<word>("active", "on");

        if (modelType == "actuatorLine")
        {
            dictionary dict;
            dict.add("axialFlowTurbineALSourceCoeffs", turbineSubDict);
            dict.add("type", "axialFlowTurbineALSource");
            dict.add("active", dict_.lookup("active"));
            axialFlowTurbines_.set
            (
                i,
                new axialFlowTurbineALSource(
                    name_ + "." + turbineName,
                    "axialFlowTurbineALSource",
                    turbineSubDict,
                    mesh_
                )
            );
        }
        else if (modelType == "actuatorDisc")
        {
            dictionary dict;
            dict.add("axialFlowTurbineALSourceCoeffs", turbineSubDict);
            dict.add("type", "axialFlowTurbineALSource");
            dict.add("active", dict_.lookup("active"));
            axialFlowTurbines_.set
            (
                i,
                new axialFlowTurbineADSource(
                    name_ + "." + turbineName,
                    "axialFlowTurbineADSource",
                    turbineSubDict,
                    mesh_
                )
            );
        }
        else
        {
            FatalErrorIn("void turbineFarmSource::createAxialFlowTurbines()")
                << "Unknown model type " << modelType
                << " Valid options are actuatorLine and actuatorDisc"
                << abort(FatalError);
        }
    }
    axialFlowTurbinePowerCoefficients_.setSize(nAxialFlowTurbines_);
    axialFlowTurbineDragCoefficients_.setSize(nAxialFlowTurbines_);
}

void Foam::fv::turbineFarmSource::createCrossFlowTurbines()
{
    // Create the turbines
    crossFlowTurbines_.setSize(nCrossFlowTurbines_);
    for (int i = 0; i < nCrossFlowTurbines_; i++)
    {
        word turbineName = crossFlowTurbineNames_[i];
        // Create dictionary items for this blade
        dictionary turbineSubDict = crossFlowTurbinesDict_.subDict(turbineName);
        word modelType
        (
            turbineSubDict.lookupOrDefault<word>("type", "actuatorLine")
        );
        turbineSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));

        // copy farm data if not specified
        word selectionMode = turbineSubDict.lookupOrAddDefault
        (
            "selectionMode",
            coeffs_.getOrDefault<word>("selectionMode", "all")
        );

        if (selectionMode != "all")
        {
            turbineSubDict.lookupOrAddDefault
            (
                "cellSet",
                coeffs_.lookup("cellSet")
            );
        }

        // Do not write force from individual turbines unless specified
        turbineSubDict.lookupOrAddDefault("writeForceField", false);

        // set active as default
        turbineSubDict.lookupOrAddDefault<word>("active", "on");

        if (modelType == "actuatorLine")
        {
            dictionary dict;
            dict.add("crossFlowTurbineALSourceCoeffs", turbineSubDict);
            dict.add("type", "crossFlowTurbineALSource");
            dict.add("active", dict_.lookup("active"));
            crossFlowTurbines_.set
            (
                i,
                new crossFlowTurbineALSource(
                    name_ + "." + turbineName,
                    "crossFlowTurbineALSource",
                    turbineSubDict,
                    mesh_
                )
            );
        }
        else if (modelType == "actuatorDisc")
        {
            dictionary dict;
            dict.add("crossFlowTurbineALSourceCoeffs", turbineSubDict);
            dict.add("type", "crossFlowTurbineALSource");
            dict.add("active", dict_.lookup("active"));
            crossFlowTurbines_.set
            (
                i,
                new crossFlowTurbineADSource(
                    name_ + "." + turbineName,
                    "crossFlowTurbineADSource",
                    turbineSubDict,
                    mesh_
                )
            );
        }
        else
        {
            FatalErrorIn("void turbineFarmSource::createCrossFlowTurbines()")
                << "Unknown model type " << modelType
                << " Valid options are actuatorLine and actuatorDisc"
                << abort(FatalError);
        }
    }
    crossFlowTurbinePowerCoefficients_.setSize(nCrossFlowTurbines_);
    crossFlowTurbineDragCoefficients_.setSize(nCrossFlowTurbines_);
}

void Foam::fv::turbineFarmSource::getTurbineData()
{
    forAll(axialFlowTurbines_, i)
    {
        axialFlowTurbinePowerCoefficients_[i] =
            axialFlowTurbines_[i].meanPowerCoefficient();
        axialFlowTurbineDragCoefficients_[i] =
            axialFlowTurbines_[i].meanDragCoefficient();
    }
    forAll(crossFlowTurbines_, i)
    {
        crossFlowTurbinePowerCoefficients_[i] =
            crossFlowTurbines_[i].meanPowerCoefficient();
        crossFlowTurbineDragCoefficients_[i] =
            crossFlowTurbines_[i].meanDragCoefficient();
    }
}

void Foam::fv::turbineFarmSource::allocateAL()
{
    forAll(axialFlowTurbines_, i)
    {
        axialFlowTurbines_[i].allocateAL();
    }
    forAll(crossFlowTurbines_, i)
    {
        crossFlowTurbines_[i].allocateAL();
    }
}

void Foam::fv::turbineFarmSource::initializeAL()
{
    // Note that there will be an additional MPI synchronization each iteration
    // However, this is only done once per simulation
    forAll(axialFlowTurbines_, i)
    {
        axialFlowTurbines_[i].initializeAL();
    }
    forAll(crossFlowTurbines_, i)
    {
        crossFlowTurbines_[i].initializeAL();
    }
}

void Foam::fv::turbineFarmSource::setupPositions(bool includeRing)
{
    forAll(axialFlowTurbines_, i)
    {
        axialFlowTurbines_[i].setupPositions(includeRing);
    }
    forAll(crossFlowTurbines_, i)
    {
        crossFlowTurbines_[i].setupPositions(includeRing);
    }
}

void Foam::fv::turbineFarmSource::calculateForces()
{
    forAll(axialFlowTurbines_, i)
    {
        axialFlowTurbines_[i].calculateForces();
    }
    forAll(crossFlowTurbines_, i)
    {
        crossFlowTurbines_[i].calculateForces();
    }
}

void Foam::fv::turbineFarmSource::createForceFieldForChildren
(
    const bool compressible
)
{
    // Create for myself (only if writeForceField_ is true)
    createForceField(false, compressible);

    forAll(axialFlowTurbines_, i)
    {
        axialFlowTurbines_[i].createForceFieldForChildren(compressible);
    }
    forAll(crossFlowTurbines_, i)
    {
        crossFlowTurbines_[i].createForceFieldForChildren(compressible);
    }
}

void Foam::fv::turbineFarmSource::addForce
(
    volVectorField &forceField,
    scalar scale,
    bool compressible
)
{
    forAll(axialFlowTurbines_, i)
    {
        axialFlowTurbines_[i].addForce
        (
            forceField,
            scale,
            compressible
        );
    }
    forAll(crossFlowTurbines_, i)
    {
        crossFlowTurbines_[i].addForce
        (
            forceField,
            scale,
            compressible
        );
    }

    getTurbineData();
    if (Pstream::master())
    {
        writePerf();
    }
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::turbineFarmSource::turbineFarmSource(
    const word &name,
    const word &modelType,
    const dictionary &dict,
    const fvMesh &mesh)
    : actuatorModelBase(name, modelType, dict, mesh),
    nAxialFlowTurbines_(0),
    nCrossFlowTurbines_(0),
    axialFlowTurbinePowerCoefficients_(0),
    axialFlowTurbineDragCoefficients_(0),
    crossFlowTurbinePowerCoefficients_(0),
    crossFlowTurbineDragCoefficients_(0)
{
    read(dict);
    createAxialFlowTurbines();
    createCrossFlowTurbines();

    // set up actuatorLines_ with all lines in simulation
    label nLines = 0;
    forAll(axialFlowTurbines_, i)
    {
        nLines += axialFlowTurbines_[i].actuatorLines().size();
    }
    forAll(crossFlowTurbines_, i)
    {
        nLines += crossFlowTurbines_[i].actuatorLines().size();
    }
    actuatorLines_.setSize(nLines);
    label index = 0;
    forAll(axialFlowTurbines_, i)
    {
        forAll(axialFlowTurbines_[i].actuatorLines(), j)
        {
            actuatorLines_[index++] = axialFlowTurbines_[i].actuatorLines()[j];
        }
    }
    forAll(crossFlowTurbines_, i)
    {
        forAll(crossFlowTurbines_[i].actuatorLines(), j)
        {
            actuatorLines_[index++] = crossFlowTurbines_[i].actuatorLines()[j];
        }
    }

    createOutputFile();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::turbineFarmSource::~turbineFarmSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::turbineFarmSource::addTurbulence
(
    fvMatrix<scalar>& eqn,
    const word fieldName
)
{
    forAll(axialFlowTurbines_, i)
    {
        axialFlowTurbines_[i].addTurbulence(eqn, fieldName);
    }
    forAll(crossFlowTurbines_, i)
    {
       crossFlowTurbines_[i].addTurbulence(eqn, fieldName);
    }
}

void Foam::fv::turbineFarmSource::writePerf()
{
    scalar combinedMean = 0.0;

    if (nAxialFlowTurbines_ + nCrossFlowTurbines_ > 0)
    {
        combinedMean =
            (sum(axialFlowTurbinePowerCoefficients_) +
            sum(crossFlowTurbinePowerCoefficients_))
            / (nAxialFlowTurbines_ + nCrossFlowTurbines_);
    }
    *outputFile_<< mesh_.time().value() << "," << combinedMean;

    // Write power coefficient and drag coefficient for each turbine
    forAll(axialFlowTurbines_, i)
    {
        *outputFile_<< "," << axialFlowTurbinePowerCoefficients_[i]
                    << "," << axialFlowTurbineDragCoefficients_[i];
    }
    forAll(crossFlowTurbines_, i)
    {
        *outputFile_<< "," << crossFlowTurbinePowerCoefficients_[i]
                    << "," << crossFlowTurbineDragCoefficients_[i];
    }

    *outputFile_<< endl;
}

bool Foam::fv::turbineFarmSource::read(const dictionary &dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Get blade information
        axialFlowTurbinesDict_ = coeffs_.subOrEmptyDict("axialFlowTurbines");
        if (axialFlowTurbinesDict_.empty())
        {
            nAxialFlowTurbines_ = 0;
        }
        else
        {
            nAxialFlowTurbines_ = axialFlowTurbinesDict_.keys().size();
            axialFlowTurbineNames_ = axialFlowTurbinesDict_.toc();
        }
        crossFlowTurbinesDict_ = coeffs_.subOrEmptyDict("crossFlowTurbines");
        if (crossFlowTurbinesDict_.empty())
        {
            nCrossFlowTurbines_ = 0;
        }
        else
        {
            nCrossFlowTurbines_ = crossFlowTurbinesDict_.keys().size();
            crossFlowTurbineNames_ = crossFlowTurbinesDict_.toc();
        }

        if (nAxialFlowTurbines_ == 0 && nCrossFlowTurbines_ == 0)
        {
            FatalErrorIn("void turbineFarmSource::read()")
                << "No turbines found in input file. "
                << "Insert axial flow turbines with axialFlowTurbines "
                << "and cross flow turbines with crossFlowTurbines"
                << abort(FatalError);
        }
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::fv::turbineFarmSource::printCoeffs() const
{
    Info<< "Number of axial flow turbines: " << nAxialFlowTurbines_ << endl;
    Info<< "Number of cross flow turbines: " << nCrossFlowTurbines_ << endl;
}

void Foam::fv::turbineFarmSource::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/farms"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/farms"
            / mesh_.time().timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_<< "time,cp";

    forAll(axialFlowTurbines_, i)
    {
        *outputFile_<< ",cp_" << axialFlowTurbineNames_[i];
        *outputFile_<< ",cd_" << axialFlowTurbineNames_[i];
    }
    forAll(crossFlowTurbines_, i)
    {
        *outputFile_<< ",cp_" << crossFlowTurbineNames_[i];
        *outputFile_<< ",cd_" << crossFlowTurbineNames_[i];
    }

    *outputFile_<< endl;
}

void Foam::fv::turbineFarmSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}

// ************************************************************************* //
