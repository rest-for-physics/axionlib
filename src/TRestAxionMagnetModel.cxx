/*************************************************************************
 * This file is part of the REST software framework.                     *
 *                                                                       *
 * Copyright (C) 2016 GIFNA/TREX (University of Zaragoza)                *
 * For more information see https://gifna.unizar.es/trex                 *
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
 * If not, see https://www.gnu.org/licenses/.                            *
 * For the list of contributors see $REST_PATH/CREDITS.                  *
 *************************************************************************/

/////////////////////////////////////////////////////////////////////////
///                                                                      
///----------------------------------------------------------------------
///                                                                      
/// REST-for-Physics - Software for Rare Event Searches Toolkit 	    
///                                                                      
/// History of developments:                                             
///                                                                      
/// 2024-04: First implementation of TRestAxionMagnetModel
/// Javier Galan
///                                                                      
/// \class TRestAxionMagnetModel                                               
/// \author Javier Galan (javier.galan@unizar.es)
///                                                                      
/// <hr>                                                                 
///                                                                      
#include <iostream>

#include "TRestAxionMagnetModel.h"
ClassImp(TRestAxionMagnetModel);

void TRestAxionMagnetModel::Test(Double_t X, Double_t Y, Double_t dl)
{
    std::vector<Double_t> bX = GetComponentAlongPath(0, TVector3 (X, Y, -10000), TVector3(X, Y, 10000), dl); 
    std::vector<Double_t> bY = GetComponentAlongPath(1, TVector3 (X, Y, -10000), TVector3(X, Y, 10000), dl); 
    std::vector<Double_t> bZ = GetComponentAlongPath(2, TVector3 (X, Y, -10000), TVector3(X, Y, 10000), dl); 

	std::vector<Double_t> z;
	for( int n = 0; n < bX.size(); n++ )
		z.push_back(dl/2 + n * dl);

	TRestAxionMagneticFit bFit;
	bFit.LoadData( z, bX, bY, bZ );
	bFit.Fit();

	fLinearFit.push_back(bFit);

	//std::cout << "Chi2 : " << bFit.GetChi2_X() << std::endl;

}

///////////////////////////////////////////////                          
/// \brief Default constructor                                          
///                                                                      
//TRestAxionMagnetModel::TRestAxionMagnetModel() { }

///////////////////////////////////////////////                          
/// \brief Default destructor                                           
///                                                                      
//TRestAxionMagnetModel::~TRestAxionMagnetModel() { }

