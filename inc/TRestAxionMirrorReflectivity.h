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

#ifndef _TRestAxionMirrorReflectivity
#define _TRestAxionMirrorReflectivity

#include <TRestMetadata.h>
#include <iostream>

/// A metadata class accessing the Henke database to load reflectivity data
class TRestAxionMirrorReflectivity : public TRestMetadata {
   private:
   protected:
   public:
    void Initialize();

    /// Pure abstract method to be implemented at inherited class
    virtual Double_t GetEfficiency(const TVector3& pos, const TVector3& dir) { return 0.0; }

    void PrintMetadata();

    TRestAxionMirrorReflectivity();
    TRestAxionMirrorReflectivity(const char* cfgFileName, std::string name = "");
    ~TRestAxionMirrorReflectivity();

    ClassDef(TRestAxionMirrorReflectivity, 1);
};
#endif
