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
/// Write the class description Here                                     
/// 
/// ### Parameters
/// Describe any parameters this process receives: 
/// * **parameter1**: This parameter ...
/// * **parameter2**: This parameter is ...
/// 
/// 
/// ### Examples
/// Give examples of usage and RML descriptions that can be tested.      
/// \code
///     <WRITE A CODE EXAMPLE HERE>
/// \endcode
/// 
/// ### Running pipeline example
/// Add the examples to a pipeline to guarantee the code will be running 
/// on future framework upgrades.                                        
/// 
/// 
/// Please, add any figure that may help to illustrate the process or metadata.  
/// 
/// \htmlonly <style>div.image img[src="image.png"]{width:500px;}</style> \endhtmlonly
/// ![A figure title description](image.png)             
/// 
/// The png image should be uploaded to the ./images/ directory          
///                                                                      
///----------------------------------------------------------------------
///                                                                      
/// REST-for-Physics - Software for Rare Event Searches Toolkit 	    
///                                                                      
/// History of developments:                                             
///                                                                      
/// YEAR-Month: First implementation of TRestAxionMagneticFit
/// WRITE YOUR FULL NAME 
///                                                                      
/// \class TRestAxionMagneticFit                                               
/// \author: TODO. Write full name and e-mail:        jgalan
///                                                                      
/// <hr>                                                                 
///                                                                      
///

#include <iostream>

#include <TGraph.h>
#include <TAxis.h>
#include <TF1.h>

#include "TRestAxionMagneticFit.h"

///////////////////////////////////////////////                          
/// \brief Default constructor                                          
///                                                                      
TRestAxionMagneticFit::TRestAxionMagneticFit() { }

///////////////////////////////////////////////                          
/// \brief Default destructor                                           
///                                                                      
TRestAxionMagneticFit::~TRestAxionMagneticFit() { }

///////////////////////////////////////////////                          
/// \brief Launches the fitting of all components in all the defined ranges
///                                                                      
void TRestAxionMagneticFit::Fit( )
{
	FitBx( );
}

///////////////////////////////////////////////                          
/// \brief It returns the value of the Bx field component using the fit parameters
///                                                                      
Double_t TRestAxionMagneticFit::Bx( const Double_t &z )
{
	Double_t Z = z;
	if( Z > 10000 )
		Z = 20000 - Z;

	if( Z < 3000 )
	{
		return BxFunction_0( &Z, &fBxParameters[0][0] );
	}
	else if ( Z < 6000 )
	{
		return BxFunction_1( &Z, &fBxParameters[1][0] );
	}
	else if ( Z <= 10000 )
	{
		return BxFunction_2( &Z, &fBxParameters[2][0] );
	}

	return 0;
}

///////////////////////////////////////////////                          
/// \brief Produces a plot with the data and the fitted profile with the present 
/// parameters status
///                                                                      
TCanvas *TRestAxionMagneticFit::DrawComponents( )
{
	if (fCanvas != NULL) {
		delete fCanvas;
		fCanvas = NULL;
	}

	fCanvas = new TCanvas("FitCanvas", "", 800, 600);

	fGraphBx->Draw("AP");
	fCanvas->Update();

	const Double_t *xData = fGraphBx->GetX();
	std::vector<Double_t> bX;
	for( int n = 0; n < fGraphBx->GetN(); n++ )
	{
		bX.push_back( Bx(xData[n]) );
	}

	TGraph *bxG = new TGraph(bX.size(), xData, &bX[0]);
	bxG->SetLineColor(kRed);
	bxG->SetLineWidth(3);
	bxG->Draw("same");

	fCanvas->Update();

	return fCanvas;
}

///////////////////////////////////////////////                          
/// \brief It reinitializes the fit parameter values to the ones defined by default
///                                                                      
void TRestAxionMagneticFit::Initialize( )
{
	fBxParameters.clear();
	Int_t nPars_0 = 5;
	Double_t rangeLow_0 = 0;
	Double_t rangeHigh_0 = 3000;
	fBxParameters.push_back({-9.25005e-06, 0.00201326, (Double_t) nPars_0, -0.000335604, -2.49173e-07, -5.33662e-11, -5.44359e-14, -9.78496e-20});

	Int_t nPars_1 = 0;
	Double_t rangeLow_1 = 3000;
	Double_t rangeHigh_1 = 6000;
	fBxParameters.push_back({5.54085730344933e-08, 15828.6874921521, -0.000340902705115058, 2.76592388034543e-13, -2.01413327616508e-17, -5.04429881940182e-21, -6.17400428116368e-17, 2.87443466840424e-24, -6.37575734100784e-22, 5425.79751717787, -0.000340902705115058, -2.72736208621692e-08, (Double_t) nPars_1 });

	Int_t nPars_2 = 5;
	Double_t rangeLow_2 = 6000;
	Double_t rangeHigh_2 = 10000;
	fBxParameters.push_back({ 0.174666, 9.33125e-6, 5980.9, (Double_t) nPars_2, 21.8151, -0.00998028 ,1.71624e-06 ,-1.30946e-10 , 3.7387e-15 });
}

///////////////////////////////////////////////                          
/// \brief It loads the data into the TGraph objects which store the field profile data
///                                                                      
void TRestAxionMagneticFit::LoadData( const std::vector <Double_t> &z, const std::vector<Double_t> &Bx, const std::vector<Double_t> &By, const std::vector<Double_t> &Bz )
{
	// Creating a TGraph object that contains the Bx profile
	if( fGraphBx != nullptr )
	{
		delete fGraphBx;
		fGraphBx = nullptr;
	}
	fGraphBx = new TGraph(z.size(), &z[0], &Bx[0]);

	// Optionally customize the appearance of the graph
	fGraphBx->SetLineColor(kRed);
	fGraphBx->SetLineWidth(2);
	fGraphBx->SetMarkerStyle(20);
	fGraphBx->SetMarkerSize(0.5);

	// Optionally, add a title and labels
	fGraphBx->SetTitle("");
	fGraphBx->GetXaxis()->SetTitle("Z [mm]");
	fGraphBx->GetYaxis()->SetTitle("B field components");
	//fGraphBx->GetXaxis()->SetLimits(rangeLow,rangeHigh);
	//fGraphBx->GetYaxis()->SetRangeUser(-1, 2 );
	

}

///////////////////////////////////////////////                          
/// \brief It launches the fitting of the magnetic field on the Bx component.
///
void TRestAxionMagneticFit::FitBx( )
{
	Int_t nPars_0 = 5;
	Double_t rangeLow_0 = 0;
	Double_t rangeHigh_0 = 3000;

	Int_t nPars_1 = 0;
	Double_t rangeLow_1 = 3000;
	Double_t rangeHigh_1 = 6000;

	Int_t nPars_2 = 5;
	Double_t rangeLow_2 = 6000;
	Double_t rangeHigh_2 = 10000;

	TF1 *fn = new TF1("fitBx0", &BxFunction_0, rangeLow_0, rangeHigh_0, 3+nPars_0);

	for( int n = 0; n < fBxParameters[0].size(); n++ )
		fn->SetParameter(n, fBxParameters[0][n] );
	fn->FixParameter(2,nPars_0);

	fGraphBx->Fit("fitBx0","R", "", rangeLow_0, rangeHigh_0);

	fBxParameters[0].clear();
	std::cout << "Bx-0 N params : " << fn->GetNpar() << std::endl;
	for( int n = 0; n < fn->GetNpar(); n++ )
	{
		std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
		fBxParameters[0].push_back(fn->GetParameter(n));
	}
	delete fn;

	fn = new TF1("fitBx1", &BxFunction_1, rangeLow_1, rangeHigh_1, 13+nPars_1);

	for( int n = 0; n < fBxParameters[1].size(); n++ )
		fn->SetParameter(n, fBxParameters[1][n] );
	fn->FixParameter(12,nPars_1);

	fGraphBx->Fit("fitBx1","R", "", rangeLow_1, rangeHigh_1);

	fBxParameters[1].clear();
	std::cout << "Bx-1 N params : " << fn->GetNpar() << std::endl;
	for( int n = 0; n < fn->GetNpar(); n++ )
	{
		std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
		fBxParameters[1].push_back(fn->GetParameter(n));
	}
	delete fn;

	/////// 

	fn = new TF1("fitBx2", &BxFunction_2, rangeLow_2, rangeHigh_2, 4+nPars_2);

	for( int n = 0; n < fBxParameters[2].size(); n++ )
		fn->SetParameter(n, fBxParameters[2][n] );
	fn->FixParameter(3,nPars_2);

	fGraphBx->Fit("fitBx2","R", "", rangeLow_2, rangeHigh_2);

	fBxParameters[2].clear();
	std::cout << "Bx-2 N params : " << fn->GetNpar() << std::endl;
	for( int n = 0; n < fn->GetNpar(); n++ )
	{
		std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
		fBxParameters[2].push_back(fn->GetParameter(n));
	}
	delete fn;
}

