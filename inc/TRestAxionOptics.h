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

#ifndef _TRestAxionOptics
#define _TRestAxionOptics

#include <TRestMetadata.h>
#include <iostream>

/// A class to load optics response files (WIP). Perhaps we might inherit later TRestAxionMCPLOptics, ...
class TRestAxionOptics : public TRestMetadata {
   private:
    /// It is the position of the optics system. Defined at the entrance of the optics
    TVector3 fCenter = TVector3(0, 0, 0);

    /// The axis of the optics system
    TVector3 fAxis = TVector3(0, 0, 1);

    void Initialize();

   public:
    void PrintMetadata();

    void InitFromConfigFile();

    TRestAxionOptics();
    TRestAxionOptics(const char* cfgFileName, std::string name = "");
    ~TRestAxionOptics();

    ClassDef(TRestAxionOptics, 1);
};
#endif
