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

#ifndef _TRestAxionTemplate
#define _TRestAxionTemplate

#include <TRestMetadata.h>

//! A metadata class to serve as example on the implementation of future metadata classes
class TRestAxionTemplate:public TRestMetadata {
    private:
        void Initialize();

        void InitFromConfigFile();

		Double_t fDummyValue; //->

    public:
		
        void PrintMetadata( );

        //Constructors
        TRestAxionTemplate();
        TRestAxionTemplate( const char *cfgFileName, std::string name = "");
        //Destructor
        ~TRestAxionTemplate();

        ClassDef(TRestAxionTemplate, 1); 
};
#endif
