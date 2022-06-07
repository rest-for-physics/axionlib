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

#ifndef _TRestAxionMCPLOptics
#define _TRestAxionMCPLOptics

#include <TRestAxionOptics.h>
#include <iostream>

/// A class to load optics response using MCPL files
class TRestAxionMCPLOptics : public TRestAxionOptics {
   private:
    void Initialize() override;

    /// The file containing the input particle list
    std::string fInputMCPLFilename;

    /// The file containing the output particle list
    std::string fOutputMCPLFilename;

   public:
    void PrintMetadata() override;

    void InitFromConfigFile() override;

    /// It returns the entrance Z-position defined by the optical axis.
    Double_t GetEntranceZPosition() override { return 0; }

    /// It returns the exit Z-position defined by the optical axis
    Double_t GetExitZPosition() override { return 0; }

    /// It updates the internal TRestAxionOptics particle positions/directions and returns effieciency
    Double_t PropagatePhoton(const TVector3& pos, const TVector3& dir, Double_t energy) override { return 0; }

    TRestAxionMCPLOptics();
    TRestAxionMCPLOptics(const char* cfgFileName, std::string name = "");
    ~TRestAxionMCPLOptics();

    ClassDefOverride(TRestAxionMCPLOptics, 1);
};
#endif
