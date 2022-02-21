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
/// TRestAxionGenericOptics is a class that inherits from TRestAxionOptics.
///
/// ToDO: Write what happens here
///
///
/// ### RML definition
///
/// We can add any number of magnetic volumes inside the RML definition
/// as shown in the following piece of code,
///
/// Example 1:
/// \code
/// <TRestAxionGenericOptics name="XMM" >
///	   <parameter name="center" value="(0,0,200)mm" />
///	   <parameter name="axis" value="(0,0.02,0.98)" />
///	   <parameter name="length" value="22cm" />
///
///	   <!-- We build mirror shells with 0.1mm thickness -->
///	   <parameter name="shellMinRadii" value="5,10,15,20,25" />
///	   <parameter name="shellMaxRadii" value="9.9,14.9,19.9,24.9,29.9" />
/// </TRestAxionGenericOptics>
/// \endcode
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: First concept and implementation of TRestAxionGenericOptics class.
///            	  Johanna von Oy
///
/// \class      TRestAxionGenericOptics
/// \author     Johanna von Oy <vonoy@physik.uni-bonn.de>
///
/// <hr>
///

#include "TRestAxionGenericOptics.h"

using namespace std;

#include "TRestPhysics.h"
using namespace REST_Physics;
ClassImp(TRestAxionGenericOptics);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionGenericOptics::TRestAxionGenericOptics() : TRestAxionOptics() { Initialize(); }

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionGenericOptics::~TRestAxionGenericOptics() {}

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
TRestAxionGenericOptics::TRestAxionGenericOptics(const char* cfgFileName, string name)
    : TRestAxionOptics(cfgFileName) {
    debug << "Entering TRestAxionGenericOptics constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionGenericOptics members
///
void TRestAxionGenericOptics::Initialize() {
    TRestAxionOptics::Initialize();

    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    // Call any initialization methods here
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionGenericOptics field members through a RML file
///
void TRestAxionGenericOptics::InitFromConfigFile() {
    TRestAxionOptics::InitFromConfigFile();

    /// TODO Initialize metadata members of this class here

    // If we recover the metadata class from ROOT file we will need to call Initialize ourselves
    this->Initialize();
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionGenericOptics
///
void TRestAxionGenericOptics::PrintMetadata() {
    TRestAxionOptics::PrintMetadata();

    metadata << "---------" << endl;
    /// Print here metadata members
    /// metadata << "xxx : " < fXXX << endl;
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
