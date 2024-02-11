#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRestAxionBufferGas.h"
#include "TRestAxionField.h"
//*******************************************************************************************************
//*** Description: This script plots the transmission probability as a function of the axion mass for a given
//***
//*** --------------
//*** Usage: restManager PlotResonances [ma_max=0.1] [ma_min=0] [Ea=4.2] [Bmag=2.5] [Lcoh=10000] [gasName=He]
//*** [vacuum=true] [n_ma=10000]
//***
//*** Being all of them optional arguments.
//*** --------------
//*** Author: Javier Galan
//*******************************************************************************************************

int REST_Axion_FieldIntegrationTests(Double_t sX = 10, Double_t sY = 10, Double_t sZ = 50, Double_t dm = 0.01, Double_t tolerance = 0.1, Double_t Ea = 4.2 ){ 

	/// Setting up magnetic field and track to evaluate
	TRestAxionMagneticField magneticField("fields.rml", "babyIAXO_HD");

	for( size_t n = 0; n < magneticField.GetNumberOfVolumes(); n++ )
		magneticField.ReMap( n, TVector3(sX,sY,sZ) );
	magneticField.SetTrack( TVector3(-110,-110,-11000), TVector3(0.01,0.01,1));
	magneticField.PrintMetadata();

	/// Setting up the gas 
	TRestAxionBufferGas* gas = new TRestAxionBufferGas();
	gas->SetGasDensity( "He", 3.267069078540181e-10);
	gas->PrintMetadata();
	Double_t ma = gas->GetPhotonMass(Ea);

	/// Setting up the axion-field and assign gas and magnetic field. 
	TRestAxionField* ax = new TRestAxionField();
	ax->AssignBufferGas(gas);
	ax->AssignMagneticField(&magneticField);
	ax->SetAxionEnergy(Ea);

	/// Enable debugging mode in axion-field
	ax->SetDebug(true);

	std::pair<double,double> prob = ax->GammaTransmissionFieldMapProbability( Ea, ma - dm, tolerance, 100, 25 );

	return 0;
}
