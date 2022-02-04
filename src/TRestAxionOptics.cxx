/******************** REST disclaimer ***********************************
 * This file is part of the REST software framework.                     *
 *                                                                       *
 * Copyright (C) 2016 GIFNA/TREX (University of Zaragoza)                *
 * For more information see http://gifna.unizar.es/trex                  *
 *                                                                       *
 * REST is free software: you can redistribute it and/or modify          *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * REST is distributed in the hope that it will be useful,               *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have a copy of the GNU General Public License along with   *
 * REST in $REST_PATH/LICENSE.                                           *
 * If not, see http://www.gnu.org/licenses/.                             *
 * For the list of contributors see $REST_PATH/CREDITS.                  *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
/// TRestAxionOptics is a class that allows to load externally
/// defined optics response files. This metadata class will define few
/// common metadata members to define the optics alignment, position, and
/// load the response data so that it can be used by
/// TRestAxionOpticsResponseProcess to produce the ray tracing of photons
/// moving up to a given plane.
///
///
/// TODO We might include opaque regions in our TRestAxionOptics. So, for
/// example we implement a ring mask, so that only photons going through
/// that mask will be considered.
///
/// TODO We could add an angle (fTiltTheta and fTiltPhi?) that helps to tilt
/// the axis by a given value.
/// We could use for example an auxiliar vector: TVector3 fAxisTilted; //!
///
/// ### RML definition
///
/// We can add any number of magnetic volumes inside the RML definition
/// as shown in the following piece of code,
///
/// Example 1:
/// \code
/// <TRestAxionOptics>
/// 	<parameter name="center" value="(0,0,950)mm" />
/// 	<parameter name="axis" value="(0,0,1)" />
/// <TRestAxionOptics/>
/// \endcode
///
/// Example 2:
/// \code
/// <TRestAxionOptics center="(0,0,950)mm" axis="(0,0,1)" />
/// \endcode
///
/// Example 3:
/// \code
/// <TRestAxionOptics center="(0,0,95)" units="cm" />
/// \endcode
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: First concept and implementation of TRestAxionOptics class.
///            	  Javier Galan
///
/// \class      TRestAxionOptics
/// \author     Javier Galan <javier.galan@unizar.es>
///
/// <hr>
///

#include "TRestAxionOptics.h"

using namespace std;

#include "TRestPhysics.h"
using namespace REST_Physics;

ClassImp(TRestAxionOptics);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionOptics::TRestAxionOptics() : TRestMetadata() { Initialize(); }

///////////////////////////////////////////////
/// \brief Constructor loading data from a config file
///
/// If no configuration path is defined using TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or
/// relative.
///
/// The default behaviour is that the config file must be specified with
/// full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
/// \param name The name of the specific metadata. It will be used to find the
/// corresponding TRestAxionMagneticField section inside the RML.
///
TRestAxionOptics::TRestAxionOptics(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    debug << "Entering TRestAxionOptics constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionOptics::~TRestAxionOptics() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionOptics members
///
void TRestAxionOptics::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fEntrance = fCenter - 0.5 * fLength * fAxis;
    fExit = fCenter + 0.5 * fLength * fAxis;
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionOptics field members through a RML file
///
void TRestAxionOptics::InitFromConfigFile() {
    this->Initialize();

    fCenter = Get3DVectorParameterWithUnits("center", TVector3(0, 0, 0));
    fAxis = Get3DVectorParameterWithUnits("axis", TVector3(0, 0, 2));
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionOptics
///
void TRestAxionOptics::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "Optics entrance: (" << fEntrance.X() << ", " << fEntrance.Y() << ", " << fEntrance.Z() << ")"
             << endl;
    metadata << "Optics center: (" << fCenter.X() << ", " << fCenter.Y() << ", " << fCenter.Z() << ")"
             << endl;
    metadata << "Optics exit: (" << fExit.X() << ", " << fExit.Y() << ", " << fExit.Z() << ")" << endl;
    metadata << "Optics axis: (" << fAxis.X() << ", " << fAxis.Y() << ", " << fAxis.Z() << ")" << endl;
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
