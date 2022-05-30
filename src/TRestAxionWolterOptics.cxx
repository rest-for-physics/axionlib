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
/// TRestAxionWolterOptics is a class that inherits from TRestAxionOptics.
///
/// ToDO: Write what happens here
///
///
/// ### RML definition
///
/// Example 1:
/// \code
/// <TRestAxionWolterOptics name="dummy">
///   	<parameter name="center" value="(0,0,200)mm" />
///		<parameter name="axis" value="(0,0.02,0.98)" />
///		<parameter name="length" value="22cm" />
///
///		<!-- We build mirror shells with 0.1mm thickness -->
///		<parameter name="shellMinRadii" value="5,10,15,20,25" />
///		<parameter name="shellMaxRadii" value="9.9,14.9,19.9,24.9,29.9" />
/// <TRestAxionWolterOptics/>
/// \endcode
///
/// Example 2:
/// \code
/// <TRestAxionWolterOptics center="(0,0,950)mm" axis="(0,0,1)" />
/// \endcode
///
/// Example 3:
/// \code
/// <TRestAxionWolterOptics center="(0,0,95)" units="cm" />
/// \endcode
///
/// \htmlonly <style>div.image img[src="Wolter.png"]{width:800px;}</style> \endhtmlonly
///
/// ![Wolter optics schematic figure](Wolter.png)
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: First concept and implementation of TRestAxionWolterOptics class.
///            	  Johanna von Oy
///
/// 2022-May: Final integration
///            	  Javier Galan
///
/// \class      TRestAxionWolterOptics
/// \author     Johanna von Oy <vonoy@physik.uni-bonn.de>
/// \author     Javier Galan <javier.galan@unizar.es>
///
/// <hr>
///

#include "TRestAxionWolterOptics.h"

using namespace std;
#include <cmath>
#include "TRestPhysics.h"
using namespace REST_Physics;
ClassImp(TRestAxionWolterOptics);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionWolterOptics::TRestAxionWolterOptics() : TRestAxionOptics() { Initialize(); }

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionWolterOptics::~TRestAxionWolterOptics() {}

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
TRestAxionWolterOptics::TRestAxionWolterOptics(const char* cfgFileName, string name)
    : TRestAxionOptics(cfgFileName) {
    RESTDebug << "Entering TRestAxionWolterOptics constructor( cfgFileName, name )" << RESTendl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionWolterOptics members
///
void TRestAxionWolterOptics::Initialize() {
    TRestAxionOptics::Initialize();

    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    // Call any initialization methods here
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionWolterOptics field members through a RML file
///
void TRestAxionWolterOptics::InitFromConfigFile() {
    TRestAxionOptics::InitFromConfigFile();

    /// TODO Initialize metadata members of this class here

    // If we recover the metadata class from ROOT file we will need to call Initialize ourselves
    this->Initialize();
}

///////////////////////////////////////////////
/// \brief This calculates the interaction Point between the X-Ray photons path and a mirror in the
/// predetermined layer
///
TVector3 TRestAxionWolterOptics::GetInteractionPoint(const TVector3& pos, const TVector3& dir) {
    Int_t layer = GetEntranceRing(pos, dir);
    Double_t r1 = fRingsRadii[layer].second;
    Double_t angle = fShellsAngle[layer];
    Double_t xSep = fShellsSep[layer];
    Double_t fLength = GetMirrLength();
    // Double_t distMirrors = 0 for first stack and (fLength + xSep) * cos(angle) for the second // distance
    // of the mirror should be added in for second stack
    Double_t r2 = r1 - fLength * sin(angle);
    // r1 = sqrt((pos[0] + s * dir[0]) * (pos[0] + s * dir[0]) +
    //    (pos[1] + s * dir[1]) * (pos[1] + s * dir[1])) -
    //  ((r2 - r1) * (pos[2] + s * dir[2] - distMirrors) / (cos(angle) * fLength))  //needs to be solved for
    //  s, maybe by incresing s for some value and comparing this to r1
    return TVector3(0, 0, 0);
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionWolterOptics
///
void TRestAxionWolterOptics::PrintMetadata() {
    TRestAxionOptics::PrintMetadata();

    RESTMetadata << "---------" << RESTendl;
    /// Print here metadata members
    /// metadata << "xxx : " < fXXX << endl;
    RESTMetadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << RESTendl;
}
