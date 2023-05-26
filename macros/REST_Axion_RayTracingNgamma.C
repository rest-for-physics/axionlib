#include "TRestDataSet.h"
#include "TRestTask.h"

#ifndef RestTask_RayTracingNgamma
#define RestTask_RayTracingNgamma

const std::string BabyIAXOWeights =
    "axionPhoton_probability,axionPhoton_transmission,boreExitGate_transmission,optics_efficiency,window_"
    "transmission";
const std::string BabyIAXOFlux = "SolarFlux*GeneratorArea*(1000*MassInterval)/Nsim";

//*******************************************************************************************************
//*** Description: This macro will add three new columns (Probability/Flux/Ngamma) to the dataset given
//*** in the argument. The probability will be the result of the products of the weights given by argument,
//*** the Flux will be a normalizing factor related to the solar flux and the number of simulated MC events,
//*** and Ngamma will be the result of multiplying the Flux by Probability.
//***
//*** The weights must be existing columns in the given dataset
//*** The normalization factor is a formula that may contain relevant quantites and existint dataset
//*** columns.
//***
//*** --------------
//*** Usage: restManager RayTracingNGamma dataset.root [weight1,weight2] [normalizationFactor]
//***
//*******************************************************************************************************
Int_t REST_Axion_RayTracingNGamma(const std::string& fname, const std::string& weights = BabyIAXOWeights,
                                  const std::string& mcFlux = BabyIAXOFlux) {
    std::vector<std::string> ws = REST_StringHelper::Split(weights, ",");

    TRestDataSet d;
    d.Import(fname);

    std::string probability = "";
    int cont = 0;
    for (const auto& w : ws) {
        if (cont > 0) probability += "*";
        probability += w;
        cont++;
    }
    // std::cout << "Probability: " << probability << std::endl;

    d.Define("Probability", probability);
    d.Define("Flux", mcFlux);                // Per axion mass and second
    d.Define("Ngamma", "Probability*Flux");  // (units: eV x s-1) Must be normalized by (and integrated to) a
                                             // given axion mass interval
                                             // --> TRestComponent

    d.GetDataFrame().Display({"Probability", "Flux", "Ngamma"})->Print();

    std::string tag = fname;
    size_t pos1 = tag.find("_");
    size_t pos2 = tag.find_last_of(".");
    tag = tag.substr(pos1 + 1, pos2 - pos1 - 1);

    std::cout << "Writting to file: "
              << "DataSet_RayTracing_Ngamma_" << tag << ".root" << std::endl;
    d.Export("DataSet_RayTracing_Ngamma_" + tag + ".root");

    return 0;
}
#endif
