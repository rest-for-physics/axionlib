/*************************************************************************
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
/// TRestAxionTemplate just a class to serve as example of a specific
/// implementation of a metadata class.
///
/// Just copy .cxx and .h files to a new class name,
/// cp TRestAxionTemplate.cxx TRestAxionSpecificMetadata.cxx
/// cp TRestAxionTemplate.h TRestAxionSpecificMetadata.h
///
/// Then, inside each file replace TRestAxionTemplate by
/// TRestAxionSpecificMetadata.
///
/// Re-run "cmake ../" at the build directory to include the new
/// class to the compilation.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionTemplate class.
///             Javier Galan
///
/// \class      TRestAxionTemplate
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionTemplate.h"
using namespace std;

ClassImp(TRestAxionTemplate);
//______________________________________________________________________________
TRestAxionTemplate::TRestAxionTemplate() : TRestMetadata() {
    // TRestAxionTemplate default constructor
  				  Initialize();
}

//______________________________________________________________________________
TRestAxionTemplate::TRestAxionTemplate(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    cout << "Entering TRestAxionTemplate constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    PrintMetadata();
}

//______________________________________________________________________________
TRestAxionTemplate::~TRestAxionTemplate() {
    // TRestAxionTemplate destructor
}

void TRestAxionTemplate::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

//______________________________________________________________________________
void TRestAxionTemplate::InitFromConfigFile() {
    this->Initialize();

    // Initialize the metadata members from a configfile
    fDummyValue = StringToDouble(GetParameter("dummy", "317"));

    if (GetVerboseLevel() >= REST_Debug) PrintMetadata();
}

void TRestAxionTemplate::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << " - Dummy metadata member : " << fDummyValue << endl;
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
