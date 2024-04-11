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
/// This class allows to perform a fit of the magnetic field profile along
/// a track parallel to the magnet axis (which is symmetric). A
/// TRestAxionMagnetModel will be built then using interpolation between
/// the different adjacent parallel lines.
///
/// This class is specifically built to reproduce the BabyIAXO magnet. There
/// are plenty of empirical fit parameters and ranges that will probably not
/// work 100% with different field maps. Specially if those maps are
/// of different extension as the one used in BabyIAXO, which is between
/// -10 meters and 10 meters. Thus, in order to work out new field maps
/// following this approach we would need to create new classes with different
/// schemes, becoming this class an abstract class that serves as interface.
/// E.g. TRestAxionMagneticFitBabyIAXO::TRestAxionMagneticFit.
///
/// The present class can be tested as follows:
///
/// \code
///		TRestAxionMagneticField field("fields.rml", "babyIAXO_2024");
///
///		// Extracting the profile for each B-component (must be parallel)
///		Double_t dl = 50;
///		std::vector<Double_t> bX = field.GetComponentAlongPath( 0, TVector3(340, 0, -10000),
/// TVector3(340, 0, 10000), dl ); 		std::vector<Double_t> bY = field.GetComponentAlongPath( 1,
/// TVector3(340, 0, -10000), TVector3(340, 0, 10000), dl ); 		std::vector<Double_t> bZ =
///field.GetComponentAlongPath( 2, TVector3(340, 0, -10000), TVector3(340, 0, 10000), dl );
///
///		std::vector<Double_t> z;
///		for( int n = 0; n < bX.size(); n++ )
///			z.push_back(dl/2 + n * dl);
///
///		TRestAxionMagneticFit bFit;
///		bFit.LoadData( z, bX, bY, bZ );
///
///		bFit.Fit( );
///
///		TCanvas *c = bFit.DrawComponents();
///		c->Print("bFit.png");
///	\endcode
///
///	The above code will produce the following plots with the fitted components and the residuals:
///
/// \htmlonly <style>div.image img[src="magnetFit.png"]{width:1200px;}</style> \endhtmlonly
///
/// ![Fit results](magnetFit.png)
///
///----------------------------------------------------------------------
///
/// REST-for-Physics - Software for Rare Event Searches Toolkit
///
/// History of developments:
///
/// 2024-04: First implementation of TRestAxionMagneticFit
/// Javier Galan
///
/// \class TRestAxionMagneticFit
/// \author Javier Galan (javier.galan@unizar.es)
///
/// <hr>
///

#include "TRestAxionMagneticFit.h"

#include <TAxis.h>
#include <TF1.h>
#include <TGraph.h>

#include <iostream>

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionMagneticFit::TRestAxionMagneticFit() {}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionMagneticFit::~TRestAxionMagneticFit() {}

///////////////////////////////////////////////
/// \brief Launches the fitting of all components in all the defined ranges
///
void TRestAxionMagneticFit::Fit() {
    Initialize();
    FitBxFullRange();
    FitByFullRange();
    FitBzFullRange();
}

////// Starts Bz component methods

///////////////////////////////////////////////
/// \brief It returns the value of the Bz field component using the fit parameters
///
///
Double_t TRestAxionMagneticFit::Bz(const Double_t& z) {
    Double_t signature = 1;
    Double_t Z = z;
    if (Z < 0 || Z > 2 * fZz) return 0;

    if (Z > fZz) {
        signature = -1;
        Z = 2 * fZz - Z;
    }

    if (Z < 4000) {
        return signature * BzFunction_0(&Z, &fBzParameters[0][0]);
    } else if (Z < 6000) {
        return signature * BzFunction_1(&Z, &fBzParameters[1][0]);
    } else if (Z <= fZz) {
        return signature * BzFunction_2(&Z, &fBzParameters[2][0]);
    }

    return 0;
}

Double_t TRestAxionMagneticFit::BzFunction(Double_t* x, Double_t* par) {
    Double_t signature = 1;
    Double_t Z = x[0];
    if (Z < 0 || Z > 2 * par[0]) return 0;

    if (Z > par[0]) {
        signature = -1;
        Z = 2 * par[0] - Z;
    }

    if (Z < 4000) {
        return signature * BzFunction_0(&Z, &par[1]);
    } else if (Z < 6000) {
        return signature * BzFunction_1(&Z, &par[10]);
    } else if (Z <= par[0]) {
        return signature * BzFunction_2(&Z, &par[21]);
    }

    return 0;
}

Double_t TRestAxionMagneticFit::BzFunction_0(Double_t* x, Double_t* par) {
    Double_t func = 0;
    Int_t N = (int)par[0];
    for (int n = 0; n < N; n++) {
        Double_t monomio = par[1 + n];
        for (int m = 0; m < n; m++) monomio *= x[0];
        func += monomio;
    }

    return func;
}

Double_t TRestAxionMagneticFit::BzFunction_1(Double_t* x, Double_t* par) {
    /// Between 0 and 2000
    Double_t X = x[0] - 4000;

    Double_t func = (par[0] + par[1] * X) / (par[2] + par[3] * (X - par[4]) * (X - par[4]));
    Int_t N = (int)par[5];
    for (int n = 0; n < N; n++) {
        Double_t monomio = par[6 + n];
        for (int m = 0; m < n; m++) monomio *= X;
        func += monomio;
    }

    return func;
}

Double_t TRestAxionMagneticFit::BzFunction_2(Double_t* x, Double_t* par) {
    /// Between 2000 and 6000
    Double_t X = x[0] - 4000;

    Double_t func = par[0] * TMath::Exp(-par[1] * (X - par[2])) + par[3] * TMath::Exp(-par[4] * (X - par[5]));

    Int_t N = (int)par[6];
    for (int n = 0; n < N; n++) {
        Double_t monomio = par[7 + n];
        for (int m = 0; m < n; m++) monomio *= X;
        func += monomio;
    }

    return func;
}

////// Starts By component methods

///////////////////////////////////////////////
/// \brief It returns the value of the By field component using the fit parameters
///
Double_t TRestAxionMagneticFit::By(const Double_t& z) {
    Double_t Z = z;
    if (Z < 0 || Z > 2 * fZy) return 0;

    if (Z > fZy) Z = 2 * fZy - Z;

    if (Z < 4000) {
        return ByFunction_0(&Z, &fByParameters[0][0]);
    } else if (Z < 6000) {
        return ByFunction_1(&Z, &fByParameters[1][0]);
    } else if (Z <= fZy) {
        return ByFunction_2(&Z, &fByParameters[2][0]);
    }

    return 0;
}

Double_t TRestAxionMagneticFit::ByFunction(Double_t* x, Double_t* par) {
    Double_t Z = x[0];
    if (Z < 0 || Z > 2 * par[0]) return 0;

    if (Z > par[0]) Z = 2 * par[0] - Z;

    if (Z < 4000) {
        return ByFunction_0(&Z, &par[1]);
    } else if (Z < 6000) {
        return ByFunction_1(&Z, &par[7]);
    } else if (Z <= par[0]) {
        return ByFunction_2(&Z, &par[21]);
    }

    return 0;
}

Double_t TRestAxionMagneticFit::ByFunction_0(Double_t* x, Double_t* par) {
    Double_t X = x[0];

    Double_t func = par[0] * TMath::Exp(par[1] * (X - par[2])) + par[3] * TMath::Exp(par[4] * (X - par[5]));

    return -func;
}

Double_t TRestAxionMagneticFit::ByFunction_1(Double_t* x, Double_t* par) {
    /// Between 0 and 2000
    Double_t X = x[0] - 4000;

    Double_t func = par[0] * TMath::ATan(par[1] * (X - par[2]) + par[3]) + par[4];
    Int_t N = (int)par[5];
    for (int n = 0; n < N; n++) {
        Double_t monomio = par[6 + n] * X;
        for (int m = 0; m < n; m++) monomio *= X;
        func += monomio;
    }

    return -func;
}

Double_t TRestAxionMagneticFit::ByFunction_2(Double_t* x, Double_t* par) {
    /// Between 2000 and 6000
    Double_t X = x[0] - 4000;
    Double_t func = par[0] - par[1] * TMath::Exp(-par[2] * (X - par[3]));
    Int_t N = (int)par[4];
    for (int n = 0; n < N; n++) {
        Double_t monomio = par[5 + n] * X;
        for (int m = 0; m < n; m++) monomio *= X;
        func += monomio;
    }
    return -func;
}

///////////////////////////////////////////////
/// \brief It returns the value of the Bx field component using the fit parameters
///
Double_t TRestAxionMagneticFit::Bx(const Double_t& z) {
    Double_t Z = z;
    if (Z < 0 || Z > 2 * fZx) return 0;

    if (Z > fZx) Z = 2 * fZx - Z;

    if (Z < 3000) {
        return BxFunction_0(&Z, &fBxParameters[0][0]);
    } else if (Z < 6000) {
        return BxFunction_1(&Z, &fBxParameters[1][0]);
    } else if (Z <= fZx) {
        return BxFunction_2(&Z, &fBxParameters[2][0]);
    }

    return 0;
}

Double_t TRestAxionMagneticFit::BxFunction(Double_t* x, Double_t* par) {
    Double_t Z = x[0];
    if (Z > par[0]) Z = 2 * par[0] - Z;

    if (Z < 0 || Z > 2 * par[0]) return 0;

    if (Z < 3000) {
        return BxFunction_0(&Z, &par[1]);
    } else if (Z < 6000) {
        return BxFunction_1(&Z, &par[9]);
    } else if (Z <= par[0]) {
        return BxFunction_2(&Z, &par[23]);
    }

    return 0;
}

Double_t TRestAxionMagneticFit::BxFunction_0(Double_t* x, Double_t* par) {
    Double_t X = x[0];

    Double_t func = par[0] * TMath::Exp(par[1] * X);
    Int_t N = (int)par[2];
    for (int n = 0; n < N; n++) {
        Double_t monomio = par[3 + n];
        for (int m = 0; m < n; m++) monomio *= x[0];
        func += monomio;
    }

    return func;
}

Double_t TRestAxionMagneticFit::BxFunction_1(Double_t* x, Double_t* par) {
    Double_t X = x[0];

    Double_t func = TMath::Exp(-par[0] * (X - par[1]) * (X - par[1]) + par[12] * X) *
                    (par[2] +
                     (par[3] + par[4] * X + par[5] * X * X) /
                         (par[6] + par[7] * X + par[8] * (X - par[9]) * (X - par[9])) +
                     par[10] + par[11] * X);

    Int_t N = (int)par[13];
    for (int n = 0; n < N; n++) {
        Double_t monomio = par[14 + n];
        for (int m = 0; m < n; m++) monomio *= x[0];
        func += monomio;
    }

    return func;
}

Double_t TRestAxionMagneticFit::BxFunction_2(Double_t* x, Double_t* par) {
    Double_t X = x[0];

    Double_t func = par[0] * TMath::Exp(-par[1] * (X - par[2]) * (X - par[2]));

    Int_t N = (int)par[3];
    for (int n = 0; n < N; n++) {
        Double_t monomio = par[4 + n];
        for (int m = 0; m < n; m++) monomio *= x[0];
        func += monomio;
    }

    return func;
}

///////////////////////////////////////////////
/// \brief Produces a plot with the data and the fitted profile with the present
/// parameters status
///
TCanvas* TRestAxionMagneticFit::DrawComponents() {
    ///////////////////// Bx component ///////////////////////////
    /// We translate to meters
	std::cout << "A" << std::endl;
    TGraph* graphToDrawX = new TGraph();
    for (int i = 0; i < fGraphBx->GetN(); ++i) {
        double x, y;
        fGraphBx->GetPoint(i, x, y);
        graphToDrawX->SetPoint(i, (x - 10000) / 1000.0, y);
    }

	std::cout << "B" << std::endl;
    TGraph* graphResidualsX = new TGraph();
    for (int i = 0; i < fGraphBx->GetN(); ++i) {
        double x, y;
        fGraphBx->GetPoint(i, x, y);
        graphResidualsX->SetPoint(i, (x - 10000) / 1000.0, y - Bx(x));
    }

	std::cout << "C" << std::endl;
    if (fCanvas != nullptr) {
        delete fCanvas;
        fCanvas = nullptr;
    }

	std::cout << "D" << std::endl;
    fCanvas = new TCanvas("FitCanvas", "", 3600, 500, 3600, 400);
    fCanvas->Divide(3, 1);
    fCanvas->cd(1);

    if (fPadBxTop != nullptr) {
        delete fPadBxTop;
        fPadBxTop = nullptr;
    }
    fPadBxTop = new TPad("BxTop", "BxTop", 0.0, 0.3, 1.0, 1.0);

    fPadBxTop->SetTopMargin(0.0);
    fPadBxTop->SetBottomMargin(0.0);
    fPadBxTop->SetLeftMargin(0.15);
    fPadBxTop->SetRightMargin(0.02);
    fPadBxTop->SetBorderMode(0);
    fPadBxTop->Draw();

	std::cout << "E" << std::endl;
    if (fPadBxBottom != nullptr) {
        delete fPadBxBottom;
        fPadBxBottom = nullptr;
    }
    fPadBxBottom = new TPad("BxBottom", "BxBottom", 0.0, 0.0, 1.0, 0.3);

    fPadBxBottom->SetTopMargin(0.0);
    fPadBxBottom->SetLeftMargin(0.15);
    fPadBxBottom->SetBottomMargin(0.38);
    fPadBxBottom->SetRightMargin(0.02);
    fPadBxBottom->SetBorderMode(0);
    fPadBxBottom->Draw();

    fPadBxTop->cd();

    /// Drawing data points Bx
    graphToDrawX->SetTitle("");
    graphToDrawX->GetXaxis()->SetLabelSize(20);
    graphToDrawX->GetXaxis()->SetLabelFont(43);
    graphToDrawX->GetXaxis()->SetTitle("Z [m]");
    graphToDrawX->GetXaxis()->SetTitleFont(43);
    graphToDrawX->GetXaxis()->SetTitleSize(20);
    graphToDrawX->GetXaxis()->SetRangeUser(-10, 10);
    graphToDrawX->SetLineWidth(2);
    graphToDrawX->SetMarkerStyle(4);
    graphToDrawX->SetMarkerSize(0.5);
    graphToDrawX->GetYaxis()->SetTitle("B_{x} [T]");
    graphToDrawX->GetYaxis()->SetTitleFont(43);
    graphToDrawX->GetYaxis()->SetTitleSize(20);
    graphToDrawX->GetYaxis()->SetTitleOffset(1.5);
    graphToDrawX->GetYaxis()->SetLabelFont(43);
    graphToDrawX->GetYaxis()->SetLabelSize(20);
    graphToDrawX->GetYaxis()->SetLabelOffset(0.02);
    graphToDrawX->GetYaxis()->SetRangeUser(-0.95, 0.5);

	std::cout << "F" << std::endl;
    //// Drawing fit function Bx
    std::vector<Double_t> xData, bX;
    for (double x = 0; x < 20000; x = x + 10) {
        xData.push_back((x - 10000) / 1000);
        bX.push_back(Bx(x));
    }

    TGraph* bxG = new TGraph(bX.size(), &xData[0], &bX[0]);
    bxG->SetTitle("");
    bxG->GetXaxis()->SetLabelSize(20);
    bxG->GetXaxis()->SetLabelFont(43);
    bxG->GetXaxis()->SetTitle("Z [m]");
    bxG->GetXaxis()->SetTitleFont(43);
    bxG->GetXaxis()->SetTitleSize(20);
    bxG->GetXaxis()->SetRangeUser(-10, 10);
    bxG->SetLineWidth(2);
    bxG->SetMarkerStyle(4);
    bxG->SetMarkerSize(0.5);
    bxG->GetYaxis()->SetTitle("B_{x} [T]");
    bxG->GetYaxis()->SetTitleFont(43);
    bxG->GetYaxis()->SetTitleSize(20);
    bxG->GetYaxis()->SetTitleOffset(1.5);
    bxG->GetYaxis()->SetLabelFont(43);
    bxG->GetYaxis()->SetLabelSize(20);
    bxG->GetYaxis()->SetLabelOffset(0.02);
    bxG->GetYaxis()->SetRangeUser(-0.95, 0.5);
    bxG->SetLineColor(kRed);
    bxG->SetLineWidth(2);
    bxG->Draw("AL");
    graphToDrawX->Draw("Psame");

    fPadBxBottom->cd();

	std::cout << "G" << std::endl;
    /// Drawing residuals
    graphResidualsX->GetXaxis()->SetLabelSize(20);
    graphResidualsX->GetXaxis()->SetLabelFont(43);
    graphResidualsX->GetXaxis()->SetTitle("Z [m]");
    graphResidualsX->GetXaxis()->SetTitleSize(20);
    graphResidualsX->GetXaxis()->SetTitleFont(43);
    graphResidualsX->GetYaxis()->SetTitle("#chi");
    graphResidualsX->GetYaxis()->SetTitleOffset(1.5);
    graphResidualsX->GetYaxis()->CenterTitle();
    graphResidualsX->GetYaxis()->SetTitleSize(20);
    graphResidualsX->GetYaxis()->SetTitleFont(43);
    graphResidualsX->GetYaxis()->SetLabelSize(20);
    graphResidualsX->GetYaxis()->SetLabelFont(43);
    graphResidualsX->GetYaxis()->SetLabelOffset(0.02);
    graphResidualsX->GetXaxis()->SetRangeUser(-10, 10);
    graphResidualsX->GetYaxis()->SetNdivisions(404);
    graphResidualsX->GetYaxis()->SetRangeUser(-0.035, 0.035);
    graphResidualsX->SetLineWidth(2);
    graphResidualsX->SetMarkerStyle(4);
    graphResidualsX->SetMarkerSize(0.5);
    graphResidualsX->Draw("AP");
    fCanvas->Update();

	std::cout << "H" << std::endl;
    ///////////////////// By component ///////////////////////////
    /// We translate to meters
    TGraph* graphToDrawY = new TGraph();
    for (int i = 0; i < fGraphBy->GetN(); ++i) {
        double x, y;
        fGraphBy->GetPoint(i, x, y);
        graphToDrawY->SetPoint(i, (x - 10000) / 1000.0, y);
    }

    TGraph* graphResidualsY = new TGraph();
    for (int i = 0; i < fGraphBy->GetN(); ++i) {
        double x, y;
        fGraphBy->GetPoint(i, x, y);
        graphResidualsY->SetPoint(i, (x - 10000) / 1000.0, y - By(x));
    }
    fCanvas->cd(2);

    if (fPadByTop != nullptr) {
        delete fPadByTop;
        fPadByTop = nullptr;
    }
    fPadByTop = new TPad("ByTop", "ByTop", 0.0, 0.3, 1.0, 1.0);

	std::cout << "I" << std::endl;
    fPadByTop->SetTopMargin(0.0);
    fPadByTop->SetBottomMargin(0.0);
    fPadByTop->SetLeftMargin(0.15);
    fPadByTop->SetRightMargin(0.02);
    fPadByTop->SetBorderMode(0);
    fPadByTop->Draw();

    if (fPadByBottom != nullptr) {
        delete fPadByBottom;
        fPadByBottom = nullptr;
    }
    fPadByBottom = new TPad("ByBottom", "ByBottom", 0.0, 0.0, 1.0, 0.3);

	std::cout << "G" << std::endl;
    fPadByBottom->SetTopMargin(0.0);
    fPadByBottom->SetLeftMargin(0.15);
    fPadByBottom->SetBottomMargin(0.38);
    fPadByBottom->SetRightMargin(0.02);
    fPadByBottom->SetBorderMode(0);
    fPadByBottom->Draw();

    fPadByTop->cd();

	std::cout << "L" << std::endl;
    /// Drawing data points By
    graphToDrawY->SetTitle("");
    graphToDrawY->GetXaxis()->SetLabelSize(20);
    graphToDrawY->GetXaxis()->SetLabelFont(43);
    graphToDrawY->GetXaxis()->SetTitle("Z [m]");
    graphToDrawY->GetXaxis()->SetTitleFont(43);
    graphToDrawY->GetXaxis()->SetTitleSize(20);
    graphToDrawY->GetXaxis()->SetRangeUser(-10, 10);
    graphToDrawY->SetLineWidth(2);
    graphToDrawY->SetMarkerStyle(4);
    graphToDrawY->SetMarkerSize(0.5);
    graphToDrawY->GetYaxis()->SetTitle("B_{y} [T]");
    graphToDrawY->GetYaxis()->SetTitleFont(43);
    graphToDrawY->GetYaxis()->SetTitleSize(20);
    graphToDrawY->GetYaxis()->SetTitleOffset(1.5);
    graphToDrawY->GetYaxis()->SetLabelFont(43);
    graphToDrawY->GetYaxis()->SetLabelSize(20);
    graphToDrawY->GetYaxis()->SetLabelOffset(0.02);

	std::cout << "M" << std::endl;
    //// Drawing fit function By (xData already filled)
    std::vector<Double_t> bY;
    for (double x = 0; x < 20000; x = x + 10) {
        bY.push_back(By(x));
    }

    TGraph* byG = new TGraph(bY.size(), &xData[0], &bY[0]);
    byG->SetLineColor(kRed);
    byG->SetLineWidth(2);
    byG->SetTitle("");
    byG->GetXaxis()->SetLabelSize(20);
    byG->GetXaxis()->SetLabelFont(43);
    byG->GetXaxis()->SetTitle("Z [m]");
    byG->GetXaxis()->SetTitleFont(43);
    byG->GetXaxis()->SetTitleSize(20);
    byG->GetXaxis()->SetRangeUser(-10, 10);
    byG->SetLineWidth(2);
    byG->SetMarkerStyle(4);
    byG->SetMarkerSize(0.5);
    byG->GetYaxis()->SetTitle("B_{y} [T]");
    byG->GetYaxis()->SetTitleFont(43);
    byG->GetYaxis()->SetTitleSize(20);
    byG->GetYaxis()->SetTitleOffset(1.5);
    byG->GetYaxis()->SetLabelFont(43);
    byG->GetYaxis()->SetLabelSize(20);
    byG->GetYaxis()->SetLabelOffset(0.02);
    byG->Draw("AL");
    graphToDrawY->Draw("Psame");

	std::cout << "N" << std::endl;
    fPadByBottom->cd();

    /// Drawing residuals
    graphResidualsY->GetXaxis()->SetLabelSize(20);
    graphResidualsY->GetXaxis()->SetLabelFont(43);
    graphResidualsY->GetXaxis()->SetTitle("Z [m]");
    graphResidualsY->GetXaxis()->SetTitleSize(20);
    graphResidualsY->GetXaxis()->SetTitleFont(43);
    graphResidualsY->GetYaxis()->SetTitle("#chi");
    graphResidualsY->GetYaxis()->SetTitleOffset(1.5);
    graphResidualsY->GetYaxis()->CenterTitle();
    graphResidualsY->GetYaxis()->SetTitleSize(20);
    graphResidualsY->GetYaxis()->SetTitleFont(43);
    graphResidualsY->GetYaxis()->SetLabelSize(20);
    graphResidualsY->GetYaxis()->SetLabelFont(43);
    graphResidualsY->GetYaxis()->SetLabelOffset(0.02);
    graphResidualsY->GetYaxis()->SetNdivisions(404);
    graphResidualsY->GetXaxis()->SetRangeUser(-10, 10);
    graphResidualsY->GetYaxis()->SetRangeUser(-0.035, 0.035);
    graphResidualsY->SetLineWidth(2);
    graphResidualsY->SetMarkerStyle(4);
    graphResidualsY->SetMarkerSize(0.5);
    graphResidualsY->Draw("AP");
    fCanvas->Update();

    ///////////////////// Bz component ///////////////////////////
    /// We translate to meters
    TGraph* graphToDrawZ = new TGraph();
    for (int i = 0; i < fGraphBz->GetN(); ++i) {
        double x, y;
        fGraphBz->GetPoint(i, x, y);
        graphToDrawZ->SetPoint(i, (x - 10000) / 1000.0, y);
    }

    TGraph* graphResidualsZ = new TGraph();
    for (int i = 0; i < fGraphBz->GetN(); ++i) {
        double x, y;
        fGraphBz->GetPoint(i, x, y);
        graphResidualsZ->SetPoint(i, (x - 10000) / 1000.0, y - Bz(x));
    }

    fCanvas->cd(3);

    if (fPadBzTop != nullptr) {
        delete fPadBzTop;
        fPadBzTop = nullptr;
    }
    fPadBzTop = new TPad("BzTop", "BzTop", 0.0, 0.3, 1.0, 1.0);

    fPadBzTop->SetTopMargin(0.0);
    fPadBzTop->SetBottomMargin(0.0);
    fPadBzTop->SetLeftMargin(0.15);
    fPadBzTop->SetRightMargin(0.02);
    fPadBzTop->SetBorderMode(0);
    fPadBzTop->Draw();

    if (fPadBzBottom != nullptr) {
        delete fPadBzBottom;
        fPadBzBottom = nullptr;
    }
    fPadBzBottom = new TPad("BzBottom", "BzBottom", 0.0, 0.0, 1.0, 0.3);

    fPadBzBottom->SetTopMargin(0.0);
    fPadBzBottom->SetLeftMargin(0.15);
    fPadBzBottom->SetBottomMargin(0.38);
    fPadBzBottom->SetRightMargin(0.02);
    fPadBzBottom->SetBorderMode(0);
    fPadBzBottom->Draw();

    fPadBzTop->cd();

    /// Drawing data points Bz
    graphToDrawZ->SetTitle("");
    graphToDrawZ->GetXaxis()->SetLabelSize(20);
    graphToDrawZ->GetXaxis()->SetLabelFont(43);
    graphToDrawZ->GetXaxis()->SetTitle("Z [m]");
    graphToDrawZ->GetXaxis()->SetTitleFont(43);
    graphToDrawZ->GetXaxis()->SetTitleSize(20);
    graphToDrawZ->GetXaxis()->SetRangeUser(-10, 10);
    graphToDrawZ->SetLineWidth(2);
    graphToDrawZ->SetMarkerStyle(4);
    graphToDrawZ->SetMarkerSize(0.5);
    graphToDrawZ->GetYaxis()->SetTitle("B_{z} [T]");
    graphToDrawZ->GetYaxis()->SetTitleFont(43);
    graphToDrawZ->GetYaxis()->SetTitleSize(20);
    graphToDrawZ->GetYaxis()->SetTitleOffset(1.5);
    graphToDrawZ->GetYaxis()->SetLabelFont(43);
    graphToDrawZ->GetYaxis()->SetLabelSize(20);
    graphToDrawZ->GetYaxis()->SetLabelOffset(0.02);

    //// Drawing fit function Bz (xData already filled)
    std::vector<Double_t> bZ;
    for (double x = 0; x < 20000; x = x + 10) {
        bZ.push_back(Bz(x));
    }

    TGraph* bzG = new TGraph(bZ.size(), &xData[0], &bZ[0]);
    bzG->SetLineColor(kRed);
    bzG->SetLineWidth(2);
    bzG->SetTitle("");
    bzG->GetXaxis()->SetLabelSize(20);
    bzG->GetXaxis()->SetLabelFont(43);
    bzG->GetXaxis()->SetTitle("Z [m]");
    bzG->GetXaxis()->SetTitleFont(43);
    bzG->GetXaxis()->SetTitleSize(20);
    bzG->GetXaxis()->SetRangeUser(-10, 10);
    bzG->SetLineWidth(2);
    bzG->SetMarkerStyle(4);
    bzG->SetMarkerSize(0.5);
    bzG->GetYaxis()->SetTitle("B_{z} [T]");
    bzG->GetYaxis()->SetTitleFont(43);
    bzG->GetYaxis()->SetTitleSize(20);
    bzG->GetYaxis()->SetTitleOffset(1.5);
    bzG->GetYaxis()->SetLabelFont(43);
    bzG->GetYaxis()->SetLabelSize(20);
    bzG->GetYaxis()->SetLabelOffset(0.02);
    bzG->Draw("AL");
    graphToDrawZ->Draw("Psame");

    fPadBzBottom->cd();

    /// Drawing residuals
    graphResidualsZ->GetXaxis()->SetLabelSize(20);
    graphResidualsZ->GetXaxis()->SetLabelFont(43);
    graphResidualsZ->GetXaxis()->SetTitle("Z [m]");
    graphResidualsZ->GetXaxis()->SetTitleSize(20);
    graphResidualsZ->GetXaxis()->SetTitleFont(43);
    graphResidualsZ->GetYaxis()->SetTitle("#chi");
    graphResidualsZ->GetYaxis()->SetTitleOffset(1.5);
    graphResidualsZ->GetYaxis()->CenterTitle();
    graphResidualsZ->GetYaxis()->SetTitleSize(20);
    graphResidualsZ->GetYaxis()->SetTitleFont(43);
    graphResidualsZ->GetYaxis()->SetLabelSize(20);
    graphResidualsZ->GetYaxis()->SetLabelFont(43);
    graphResidualsZ->GetYaxis()->SetLabelOffset(0.02);
    graphResidualsZ->GetYaxis()->SetNdivisions(404);
    graphResidualsZ->GetXaxis()->SetRangeUser(-10, 10);
    graphResidualsZ->GetYaxis()->SetRangeUser(-0.035, 0.035);
    graphResidualsZ->SetLineWidth(2);
    graphResidualsZ->SetMarkerStyle(4);
    graphResidualsZ->SetMarkerSize(0.5);
    graphResidualsZ->Draw("AP");
    fCanvas->Update();

    return fCanvas;
}

///////////////////////////////////////////////
/// \brief It reinitializes the fit parameter values to the ones defined by default
///
void TRestAxionMagneticFit::Initialize() {
    fBxParameters.clear();
    Int_t nPars_0 = 5;
    // Range from 0 to 3000
    fBxParameters.push_back({-9.25005e-06, 0.00201326, (Double_t)nPars_0, -0.000335604, -2.49173e-07,
                             -5.33662e-11, -5.44359e-14, -9.78496e-20});

    Int_t nPars_1 = 0;
    // Range from 3000 to 6000
    fBxParameters.push_back({5.54095e-08, 15842.3, -34.7965, 2.76593e-13, -2.00733e-17, -5.03095e-21,
                             -6.239e-17, -1.298e-22, -6.709e-22, 5416.52, -34.7965, 0.0140181, -5.41826e-8,
                             (Double_t)nPars_1});

    Int_t nPars_2 = 5;
    // Range from 6000 to 10000
    fBxParameters.push_back({0.174666, 9.33125e-6, 5980.9, (Double_t)nPars_2, 21.8151, -0.00998028,
                             1.71624e-06, -1.30946e-10, 3.7387e-15});

    fByParameters.clear();
    // Range from 0 to 4000
    fByParameters.push_back({0.00114352, 0.000269059, 3692.87, 0.00235922, 0.00175728, 2700.33});

    Int_t nPars_1Y = 8;
    // Range from 4000 to 6000
    fByParameters.push_back({0.878548, 0.00505303, 1357.39, -0.227327, 1.28017, (Double_t)nPars_1Y,
                             -1.92992e-05, -2.97035e-08, -1.24648e-11, -5.1324e-16, 3.61321e-18, 3.84381e-21,
                             5.15085e-25, 2.13876e-29});

    Int_t nPars_2Y = 4;
    // Range from 6000 to 10000
    fByParameters.push_back({2.67728, 1.33661, 0.00329388, 1364.57, (Double_t)nPars_2Y, 2.41991e-06,
                             1.10521e-09, 8.83841e-14, -3.64002e-17});

    fBzParameters.clear();

    Int_t nPars_0Z = 8;
    // Range from 0 to 2000
    fBzParameters.push_back({(Double_t)nPars_0Z, 0.00231608, 1.25227e-06, 6.68873e-10, -9.78742e-14,
                             1.39871e-16, 9.66594e-22, 4.84975e-25, 2.45398e-28});

    Int_t nPars_1Z = 5;
    // Range from 2000 to 4000
    fBzParameters.push_back({9.82286, -0.00798509, -3.36196, -4.39675e-05, 1404.14, (Double_t)nPars_1Z,
                             0.179615, 0.000249106, -2.90378e-07, 7.39952e-10, -3.20937e-13});

    Int_t nPars_2Z = 4;
    // Range from 4000 to 10000
    fBzParameters.push_back({5.07839, 0.00700772, 1623.44, 0.61256, 0.00155811, 1522.93, (Double_t)nPars_2Z,
                             0.00399185, 1.54701e-06, 4.82946e-16, -6.33159e-14});
}

///////////////////////////////////////////////
/// \brief It loads the data into the TGraph objects which store the field profile data
///
void TRestAxionMagneticFit::LoadData(const std::vector<Double_t>& z, const std::vector<Double_t>& bx,
                                     const std::vector<Double_t>& by, const std::vector<Double_t>& bz) {
    // Creating a TGraph object that contains the Bx profile
    if (fGraphBx != nullptr) {
        delete fGraphBx;
        fGraphBx = nullptr;
    }
    fGraphBx = new TGraph(z.size(), &z[0], &bx[0]);

    // Creating a TGraph object that contains the By profile
    if (fGraphBy != nullptr) {
        delete fGraphBy;
        fGraphBy = nullptr;
    }
    fGraphBy = new TGraph(z.size(), &z[0], &by[0]);

    // Creating a TGraph object that contains the By profile
    if (fGraphBz != nullptr) {
        delete fGraphBz;
        fGraphBz = nullptr;
    }
    fGraphBz = new TGraph(z.size(), &z[0], &bz[0]);
}

///////////////////////////////////////////////
/// \brief It launches the fitting of the magnetic field on the Bx component.
/// Only used for debugging purposes. It fits independently in each range.
///
void TRestAxionMagneticFit::FitBy() {
    Int_t nPars_0 = 0;
    Double_t rangeLow_0 = 0;
    Double_t rangeHigh_0 = 4000;

    Int_t nPars_1 = 8;
    Double_t rangeLow_1 = 4000;
    Double_t rangeHigh_1 = 6000;

    Int_t nPars_2 = 4;
    Double_t rangeLow_2 = 6000;
    Double_t rangeHigh_2 = 10000;

    TF1* fn = new TF1("fitBy0", &ByFunction_0, rangeLow_0, rangeHigh_0, 6);

    for (size_t n = 0; n < fByParameters[0].size(); n++) fn->SetParameter(n, fByParameters[0][n]);

    fGraphBy->Fit("fitBy0", "R", "", rangeLow_0, rangeHigh_0);

    fByParameters[0].clear();
    std::cout << "By-0 N params : " << fn->GetNpar() << std::endl;
    for (int n = 0; n < fn->GetNpar(); n++) {
        std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fByParameters[0].push_back(fn->GetParameter(n));
    }
    delete fn;

    fn = new TF1("fitBy1", &ByFunction_1, rangeLow_1, rangeHigh_1, 6 + nPars_1);

    for (size_t n = 0; n < fByParameters[1].size(); n++) fn->SetParameter(n, fByParameters[1][n]);
    fn->FixParameter(5, nPars_1);

    fGraphBy->Fit("fitBy1", "QR", "", rangeLow_1, rangeHigh_1);

    fByParameters[1].clear();
    std::cout << "By-1 N params : " << fn->GetNpar() << std::endl;
    for (int n = 0; n < fn->GetNpar(); n++) {
        std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fByParameters[1].push_back(fn->GetParameter(n));
    }
    delete fn;

    ///////
    fn = new TF1("fitBy2", &ByFunction_2, rangeLow_2, rangeHigh_2, 5 + nPars_2);

    for (size_t n = 0; n < fByParameters[2].size(); n++) fn->SetParameter(n, fByParameters[2][n]);
    fn->FixParameter(4, nPars_2);

    fGraphBy->Fit("fitBy2", "QR", "", rangeLow_2, rangeHigh_2);

    fByParameters[2].clear();
    std::cout << "By-2 N params : " << fn->GetNpar() << std::endl;
    for (int n = 0; n < fn->GetNpar(); n++) {
        std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fByParameters[2].push_back(fn->GetParameter(n));
    }
    delete fn;
}

///////////////////////////////////////////////
/// \brief It launches the fitting of the magnetic field on the Bx component.
/// Only used for debugging purposes. It fits independently in each range.
///
void TRestAxionMagneticFit::FitBx() {
    Int_t nPars_0 = 5;
    Double_t rangeLow_0 = 0;
    Double_t rangeHigh_0 = 3000;

    Int_t nPars_1 = 0;
    Double_t rangeLow_1 = 3000;
    Double_t rangeHigh_1 = 6000;

    Int_t nPars_2 = 5;
    Double_t rangeLow_2 = 6000;
    Double_t rangeHigh_2 = 10000;

    TF1* fn = new TF1("fitBx0", &BxFunction_0, rangeLow_0, rangeHigh_0, 3 + nPars_0);

    for (size_t n = 0; n < fBxParameters[0].size(); n++) fn->SetParameter(n, fBxParameters[0][n]);
    fn->FixParameter(2, nPars_0);

    fGraphBx->Fit("fitBx0", "QR", "", rangeLow_0, rangeHigh_0);

    fBxParameters[0].clear();
    std::cout << "Bx-0 N params : " << fn->GetNpar() << std::endl;
    for (int n = 0; n < fn->GetNpar(); n++) {
        std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBxParameters[0].push_back(fn->GetParameter(n));
    }
    delete fn;

    fn = new TF1("fitBx1", &BxFunction_1, rangeLow_1, rangeHigh_1, 14 + nPars_1);

    for (size_t n = 0; n < fBxParameters[1].size(); n++) fn->SetParameter(n, fBxParameters[1][n]);
    fn->FixParameter(13, nPars_1);

    fGraphBx->Fit("fitBx1", "QR", "", rangeLow_1, rangeHigh_1);

    fBxParameters[1].clear();
    std::cout << "Bx-1 N params : " << fn->GetNpar() << std::endl;
    for (int n = 0; n < fn->GetNpar(); n++) {
        std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBxParameters[1].push_back(fn->GetParameter(n));
    }
    delete fn;

    ///////

    fn = new TF1("fitBx2", &BxFunction_2, rangeLow_2, rangeHigh_2, 4 + nPars_2);

    for (size_t n = 0; n < fBxParameters[2].size(); n++) fn->SetParameter(n, fBxParameters[2][n]);
    fn->FixParameter(3, nPars_2);

    fGraphBx->Fit("fitBx2", "QR", "", rangeLow_2, rangeHigh_2);

    fBxParameters[2].clear();
    std::cout << "Bx-2 N params : " << fn->GetNpar() << std::endl;
    for (int n = 0; n < fn->GetNpar(); n++) {
        std::cout << "Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBxParameters[2].push_back(fn->GetParameter(n));
    }
    delete fn;
}

///////////////////////////////////////////////
/// \brief It launches the fitting of the magnetic field on the Bx component.
///
void TRestAxionMagneticFit::FitBxFullRange() {
    std::vector<Double_t> parList;
    parList.push_back(10000);

    parList.insert(parList.end(), fBxParameters[0].begin(), fBxParameters[0].end());
    parList.insert(parList.end(), fBxParameters[1].begin(), fBxParameters[1].end());
    parList.insert(parList.end(), fBxParameters[2].begin(), fBxParameters[2].end());

    TF1* fn = new TF1("fitBx", &BxFunction, 0, 20000, 32);

    for (size_t n = 0; n < parList.size(); n++) fn->SetParameter(n, parList[n]);

    fn->SetParLimits(0, 9800, 10200);

    fn->FixParameter(3, 5);
    fn->FixParameter(21, 0);
    fn->FixParameter(26, 5);

    fGraphBx->Fit("fitBx", "R", "", 0, 20000);
    fChiSquareBx = fGraphBx->Chisquare(fn, "R");

    fZx = fn->GetParameter(0);

    fBxParameters[0].clear();
    std::cout << "Bx-0 N params : " << fn->GetNpar() << std::endl;
    for (int n = 1; n < 9; n++) {
        std::cout << "Bx-0 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBxParameters[0].push_back(fn->GetParameter(n));
    }

    fBxParameters[1].clear();
    for (int n = 9; n < 23; n++) {
        std::cout << "Bx-1 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBxParameters[1].push_back(fn->GetParameter(n));
    }

    fBxParameters[2].clear();
    for (int n = 23; n < 32; n++) {
        std::cout << "Bx-2 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBxParameters[2].push_back(fn->GetParameter(n));
    }

    delete fn;
}

///////////////////////////////////////////////
/// \brief It launches the fitting of the magnetic field on the By component.
///
void TRestAxionMagneticFit::FitByFullRange() {
    std::vector<Double_t> parList;
    parList.push_back(fZy);

    parList.insert(parList.end(), fByParameters[0].begin(), fByParameters[0].end());
    parList.insert(parList.end(), fByParameters[1].begin(), fByParameters[1].end());
    parList.insert(parList.end(), fByParameters[2].begin(), fByParameters[2].end());

    TF1* fn = new TF1("fitBy", &ByFunction, 0, 20000, 30);

    for (size_t n = 0; n < parList.size(); n++) {
        std::cout << "Setting parameter p" << n << ": " << parList[n] << std::endl;
        fn->SetParameter(n, parList[n]);
    }

    fn->SetParLimits(0, 9800, 10200);

    fn->FixParameter(12, 8);
    fn->FixParameter(25, 0);
    fn->SetParLimits(24, 1000, 1600);

    fGraphBy->Fit("fitBy", "R", "", 0, 20000);
    fChiSquareBy = fGraphBy->Chisquare(fn, "R");

    fZy = fn->GetParameter(0);
    std::cout << "fZy : " << fZy << std::endl;
    fByParameters[0].clear();
    std::cout << "By-0 N params : " << fn->GetNpar() << std::endl;
    for (int n = 1; n < 7; n++) {
        std::cout << "By-0 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fByParameters[0].push_back(fn->GetParameter(n));
    }

    fByParameters[1].clear();
    for (int n = 7; n < 21; n++) {
        std::cout << "By-1 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fByParameters[1].push_back(fn->GetParameter(n));
    }

    fByParameters[2].clear();
    for (int n = 21; n < 30; n++) {
        std::cout << "By-2 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fByParameters[2].push_back(fn->GetParameter(n));
    }

    delete fn;
}

///////////////////////////////////////////////
/// \brief It launches the fitting of the magnetic field on the By component.
///
void TRestAxionMagneticFit::FitBzFullRange() {
    std::vector<Double_t> parList;
    parList.push_back(fZz);

    parList.insert(parList.end(), fBzParameters[0].begin(), fBzParameters[0].end());
    parList.insert(parList.end(), fBzParameters[1].begin(), fBzParameters[1].end());
    parList.insert(parList.end(), fBzParameters[2].begin(), fBzParameters[2].end());

    TF1* fn = new TF1("fitBz", &BzFunction, 0, 20000, 33);

    for (size_t n = 0; n < parList.size(); n++) {
        std::cout << "Setting parameter p" << n << ": " << parList[n] << std::endl;
        fn->SetParameter(n, parList[n]);
    }

    fn->SetParLimits(0, 9800, 10200);

    fn->FixParameter(1, 8);
    fn->FixParameter(15, 5);
    fn->FixParameter(27, 5);

	fGraphBz->Fit("fitBz","R", "", 0, 20000);
	fChiSquareBz = fGraphBz->Chisquare( fn, "R");

    fZz = fn->GetParameter(0);
    std::cout << "fZz : " << fZz << std::endl;
    fBzParameters[0].clear();
    std::cout << "Bz-0 N params : " << fn->GetNpar() << std::endl;
    for (int n = 1; n < 10; n++) {
        std::cout << "Bz-0 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBzParameters[0].push_back(fn->GetParameter(n));
    }

    fBzParameters[1].clear();
    for (int n = 10; n < 21; n++) {
        std::cout << "Bz-1 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBzParameters[1].push_back(fn->GetParameter(n));
    }

    fBzParameters[2].clear();
    for (int n = 21; n < 33; n++) {
        std::cout << "By-2 Setting parameter " << n << " : " << fn->GetParameter(n) << std::endl;
        fBzParameters[2].push_back(fn->GetParameter(n));
    }

    delete fn;
}

Double_t TRestAxionMagneticFit::GetChi2_X()
{
	Double_t x,y;
	Double_t chi2 = 0;
	for (int i = 0; i < fGraphBx->GetN(); ++i) {
		std::cout << x << " " << y << std::endl;
		fGraphBx->GetPoint(i, x, y);
		Double_t res = (y - Bx(x));
		chi2 += res * res;
	}
	if( fGraphBx->GetN()<=0 ) return 0;
	return chi2/fGraphBx->GetN();
}
