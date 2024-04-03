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

#ifndef REST_TRestAxionMagneticFit
#define REST_TRestAxionMagneticFit

#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>

/// A class to fit and store the fit parameters from the magnetic field map
class TRestAxionMagneticFit {
private:

	/// It will be enabled once the fit has been performed
	Bool_t fFitReady = false;  //<

	/// The parameter used to fit the Bx component. Each parameters vector defines different ranges.
	std::vector<std::vector<Double_t>> fBxParameters; //<

	/// A canvas where the components can be plotted
	TCanvas *fCanvas = nullptr;

	/// A graph where to place the data to be fitted
	TGraph *fGraphBx = nullptr;

	static Double_t BxFunction_0(Double_t *x, Double_t *par)
	{

		Double_t X = x[0];

		Double_t func = par[0]*TMath::Exp(par[1]*X);
		Int_t N = (int) par[2];
		for( int n = 0; n < N; n++ )
		{
			Double_t monomio = par[3+n];
			for( int m = 0; m < n; m++ )
				monomio *= x[0];	
			func += monomio;
		}

		return func;
	}

	static Double_t BxFunction_1(Double_t *x, Double_t *par)
	{

		Double_t X = x[0];

		Double_t func = TMath::Exp( -par[0]*(X-par[1])*(X-par[1]) ) * ( par[2] + (par[3] + par[4]*X + par[5]*X*X)/(par[6]+par[7]*X+par[8]*(X-par[9])*(X-par[9])) + par[10] + par[11]*X );

		Int_t N = (int) par[12];
		for( int n = 0; n < N; n++ )
		{
			Double_t monomio = par[13+n];
			for( int m = 0; m < n; m++ )
				monomio *= x[0];	
			func += monomio;
		}

		return func;
	}

	static Double_t BxFunction_2(Double_t *x, Double_t *par)
	{


		Double_t X = x[0];

		Double_t func = par[0] * TMath::Exp( -par[1]*(X-par[2])*(X-par[2]) );

		Int_t N = (int) par[3];
		for( int n = 0; n < N; n++ )
		{
			Double_t monomio = par[4+n];
			for( int m = 0; m < n; m++ )
				monomio *= x[0];	
			func += monomio;
		}

		return func;
	}

public:

	Double_t Bx( const Double_t &z );

	void Initialize( );
	TCanvas *DrawComponents( );

	void FitBx( );
	void Fit( );
	void LoadData( const std::vector <Double_t> &z, const std::vector<Double_t> &Bx, const std::vector<Double_t> &By, const std::vector<Double_t> &Bz );

    TRestAxionMagneticFit();
    ~TRestAxionMagneticFit();
};
#endif
