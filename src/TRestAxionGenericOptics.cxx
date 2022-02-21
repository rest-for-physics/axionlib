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



