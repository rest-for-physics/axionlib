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

/***************** DOXYGEN DOCUMENTATION ********************************
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// \class      TRestAxionSpectrum
/// \author     Sebastian Hoof
///
/// <hr>
///
 *************************************************************************/

#include "TRestAxionSpectrum.h"

ClassImp(TRestAxionSpectrum);

TRestAxionSpectrum::TRestAxionSpectrum() : TRestMetadata() {
    // TRestAxionSolarModel default constructor
    Initialize();
}

TRestAxionSpectrum::TRestAxionSpectrum(const char* cfgFileName, std::string name)
    : TRestMetadata(cfgFileName) {
    cout << "Entering TRestAxionSpectrum constructor (cfgFileName, name)" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    PrintMetadata();
}

TRestAxionSpectrum::~TRestAxionSpectrum() {
    // TRestAxionSolarModel destructor
}

void TRestAxionSpectrum::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}


void TRestAxionSpectrum::InitFromConfigFile() {
    debug << "Entering TRestAxionSolarModel::InitFromConfigFile" << endl;

    this->Initialize();

    // Initialize the metadata members from a configfile

    // fClassMember = GetParameter( "paramName", "defaultValue" );


    PrintMetadata();
}

void TRestAxionSpectrum::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << " - Metadata from file header : " << fMetaDataFromFileHeader << endl;
}

double TRestAxionSpectrum::GetSolarAxionFlux(double erg_lo, double erg_hi, double er_step_size) { return 0; }
double TRestAxionSpectrum::GetDifferentialSolarAxionFlux(double erg) { return 0; }
