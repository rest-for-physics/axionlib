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

	/// A canvas where the components can be plotted
	TCanvas *fCanvas = nullptr;

	/// The Z value at which we apply field symmetry. It is updated by the fit.
	Double_t fZx = 10000;
	
	/// Chi-square for Bx component
	Double_t fChiSquareBx = 0;

	/// The parameter used to fit the Bx component. Each parameters vector defines different ranges.
	std::vector<std::vector<Double_t>> fBxParameters; //<
	
	/// A pad to plot the Bx component absolute
	TPad *fPadBxTop = nullptr;
	
	/// A pad to plot the Bx component residuals
	TPad *fPadBxBottom = nullptr;

	/// A graph where to place the data to be fitted
	TGraph *fGraphBx = nullptr;

	static Double_t BxFunction(Double_t *x, Double_t *par);
	static Double_t BxFunction_0(Double_t *x, Double_t *par);
	static Double_t BxFunction_1(Double_t *x, Double_t *par);
	static Double_t BxFunction_2(Double_t *x, Double_t *par);
	
	/// The Z value at which we apply field symmetry for By component. It is updated by the fit.
	Double_t fZy = 10000;

	/// Chi-square for By component
	Double_t fChiSquareBy = 0;

	/// The parameter used to fit the Bx component. Each parameters vector defines different ranges.
	std::vector<std::vector<Double_t>> fByParameters; //<
	
	/// A pad to plot the Bx component absolute
	TPad *fPadByTop = nullptr;
	
	/// A pad to plot the Bx component residuals
	TPad *fPadByBottom = nullptr;

	/// A graph where to place the data to be fitted
	TGraph *fGraphBy = nullptr;

	static Double_t ByFunction(Double_t *x, Double_t *par);
	static Double_t ByFunction_0(Double_t *x, Double_t *par);
	static Double_t ByFunction_1(Double_t *x, Double_t *par);
	static Double_t ByFunction_2(Double_t *x, Double_t *par);
	
	/// The Z value at which we apply field symmetry for Bz component. It is updated by the fit.
	Double_t fZz = 10000;

	/// Chi-square for Bz component
	Double_t fChiSquareBz = 0;

	/// The parameter used to fit the Bx component. Each parameters vector defines different ranges.
	std::vector<std::vector<Double_t>> fBzParameters; //<
	
	/// A pad to plot the Bx component absolute
	TPad *fPadBzTop = nullptr;
	
	/// A pad to plot the Bx component residuals
	TPad *fPadBzBottom = nullptr;

	/// A graph where to place the data to be fitted
	TGraph *fGraphBz = nullptr;

	static Double_t BzFunction(Double_t *x, Double_t *par);
	static Double_t BzFunction_0(Double_t *x, Double_t *par);
	static Double_t BzFunction_1(Double_t *x, Double_t *par);
	static Double_t BzFunction_2(Double_t *x, Double_t *par);

public:

	Double_t Bx( const Double_t &z );
	Double_t By( const Double_t &z );
	Double_t Bz( const Double_t &z );

	void Initialize( );
	TCanvas *DrawComponents( );

	void FitBxFullRange( );
    void FitByFullRange( );
    void FitBzFullRange( );

	void FitBx( );
	void FitBy( );

	void Fit( );
	void LoadData( const std::vector <Double_t> &z, const std::vector<Double_t> &Bx, const std::vector<Double_t> &By, const std::vector<Double_t> &Bz );

    TRestAxionMagneticFit();
    ~TRestAxionMagneticFit();
};
#endif
