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
/// TRestAxionMCPLOptics is a class that inherits from TRestAxionOptics.
/// It will load the optics response from a MCPL file, and inherit the common
/// optics description.
///
/// This method will implement the pure virtual method PropagatePhoton, that
/// will be responsible to move a particle to the exit plane of the optics.
/// It may use the common methods defined at TRestAxionOptics.
///
/// ### RML definition
///
/// We can add any number of magnetic volumes inside the RML definition
/// as shown in the following piece of code,
///
/// Example 1:
/// \code
/// <TRestAxionMCPLOptics name="mcpl" >
///	   <parameter name="center" value="(0,0,200)mm" />
///	   <parameter name="axis" value="(0,0.02,0.98)" />
///	   <parameter name="length" value="22cm" />
///
///	   <!-- We build mirror shells with 0.1mm thickness -->
///	   <parameter name="shellMinRadii" value="5,10,15,20,25" />
///	   <parameter name="shellMaxRadii" value="9.9,14.9,19.9,24.9,29.9" />
/// </TRestAxionMCPLOptics>
/// \endcode
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: First concept and implementation of TRestAxionMCPLOptics class.
///            	  Javier Galan
///
/// \class      TRestAxionMCPLOptics
/// \author     Javier Galan <javier.galan@unizar.es>
///
/// <hr>
///

#include "TRestAxionMCPLOptics.h"
#include "mcpl.h"

using namespace std;

#include "TRestPhysics.h"
using namespace REST_Physics;

ClassImp(TRestAxionMCPLOptics);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionMCPLOptics::TRestAxionMCPLOptics() : TRestAxionOptics() { Initialize(); }

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
TRestAxionMCPLOptics::TRestAxionMCPLOptics(const char* cfgFileName, string name)
    : TRestAxionOptics(cfgFileName) {
    debug << "Entering TRestAxionMCPLOptics constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionMCPLOptics::~TRestAxionMCPLOptics() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionMCPLOptics members
///
void TRestAxionMCPLOptics::Initialize() {
    TRestAxionOptics::Initialize();

    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    if (fInputMCPLFilename != "" && fOutputMCPLFilename != "") LoadMCPLFiles();
}

///////////////////////////////////////////////
/// \brief Method to load inside the class the MCPL files
///
Int_t TRestAxionMCPLOptics::LoadMCPLFiles() {
    string fname = fInputMCPLFilename;
    if (fname.find("https") == 0) fname = TRestTools::DownloadRemoteFile(fname);
    mcpl_file_t f = mcpl_open_file(fname.c_str());

    string fname2 = fOutputMCPLFilename;
    if (fname2.find("https") == 0) fname2 = TRestTools::DownloadRemoteFile(fname2);
    mcpl_file_t g = mcpl_open_file(fname2.c_str());

    const mcpl_particle_t *p, *q;
    while ((p = mcpl_read(f))) {
        q = mcpl_read(g);

        cout << "Input MCPL file (detector_plane)" << endl;
        cout << "--------------------------------" << endl;
        cout << "X: " << p->position[0] << " Y: " << p->position[1] << " Z: " << p->position[2] << endl;
        cout << "pX: " << p->direction[0] << " pY: " << p->direction[1] << " pZ: " << p->direction[2] << endl;
        cout << "polX: " << p->polarisation[0] << " polY: " << p->polarisation[1]
             << " polZ: " << p->polarisation[2] << endl;
        cout << "Ekin = " << p->ekin << endl;
        cout << "Time = " << p->time << endl;
        cout << "Weight = " << p->weight << endl;
        cout << "User flags = " << p->userflags << endl;
        cout << endl;
        cout << "Output MCPL file (post_optics)" << endl;
        cout << "--------------------------------" << endl;
        cout << "X: " << q->position[0] << " Y: " << q->position[1] << " Z: " << q->position[2] << endl;
        cout << "pX: " << q->direction[0] << " pY: " << q->direction[1] << " pZ: " << q->direction[2] << endl;
        cout << "polX: " << q->polarisation[0] << " polY: " << q->polarisation[1]
             << " polZ: " << q->polarisation[2] << endl;
        cout << "Ekin = " << q->ekin << endl;
        cout << "Time = " << q->time << endl;
        cout << "Weight = " << q->weight << endl;
        cout << "User flags = " << q->userflags << endl;
        cout << endl;

        getchar();
    }
    mcpl_close_file(f);

    return 0;
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionMCPLOptics field members through a RML file
///
void TRestAxionMCPLOptics::InitFromConfigFile() {
    TRestAxionOptics::InitFromConfigFile();

    fInputMCPLFilename = GetParameter("inputMCPLFilename", "none");
    fOutputMCPLFilename = GetParameter("outputMCPLFilename", "none");

    // If we recover the metadata class from ROOT file we will need to call Initialize ourselves
    this->Initialize();
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionMCPLOptics
///
void TRestAxionMCPLOptics::PrintMetadata() {
    TRestAxionOptics::PrintMetadata();

    metadata << "---------" << endl;
    metadata << "Input MCPL file: " << fInputMCPLFilename << endl;
    metadata << "Output MCPL file: " << fOutputMCPLFilename << endl;
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
