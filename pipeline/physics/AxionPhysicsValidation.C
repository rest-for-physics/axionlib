#include <TRestAxionField.h>
#include <TRestAxionMagneticField.h>

Int_t AxionPhysicsValidation() {
    std::cout << endl;

    std::cout << "Evaluating magnetic field for BabyIAXO" << std::endl;
    std::cout << "--------------------------------------" << std::endl;

    TRestAxionMagneticField magneticField("fields.rml", "babyIAXO");

    std::cout << "Field at (0,0,0) : " << magneticField.GetMagneticField(0, 0, 0).Y() << std::endl;
    std::cout << "Field at (0,0,-4000) : " << magneticField.GetMagneticField(0, 0, -4000).Y() << std::endl;
    std::cout << "Field at (0,300,4000) : " << magneticField.GetMagneticField(0, 300, 4000).Y() << std::endl;

    if (std::round(magneticField.GetMagneticField(0, 0, 0).Y() * 1000) != -2007) {
        std::cout << "Field at (0,0,0) is not -2.07 T!" << std::endl;
        return 10;
    }

    if (std::round(magneticField.GetMagneticField(0, 0, -4000).Y() * 1000) != -1700) {
        std::cout << "Field at (0,0,-4000) is not -1.7 T!" << std::endl;
        return 11;
    }

    if (std::round(magneticField.GetMagneticField(0, 300, 4000).Y() * 1000) != -1290) {
        std::cout << "Field at (0,0,4000) is not -1.29 T!" << std::endl;
        return 12;
    }

    std::cout << endl;
    std::cout << "--> Magnetic field evaluation succeeded!" << std::endl;
    std::cout << endl;

    std::cout << "Evaluating dummy gas mixtures" << std::endl;
    std::cout << "-----------------------------" << std::endl;

    TRestAxionBufferGas helium("bufferGases.rml", "helium");

    std::cout << "Gas name : " << helium.GetName() << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Number of gases : " << helium.GetNumberOfGases() << std::endl;
    std::cout << "Photon mass at 3keV : " << helium.GetPhotonMass(3) << " eV" << std::endl;
    std::cout << "Absorption lenght (cm-1) at 3keV : " << helium.GetPhotonAbsorptionLength(3) << std::endl;

    if (std::round(1.e6 * helium.GetPhotonMass(3)) != 1017) {
        std::cout << "Photon mass is not 1.017 meV!" << std::endl;
        return 21;
    }

    if (std::round(1.e12 * helium.GetPhotonAbsorptionLength(3)) != 4720) {
        std::cout << "Absorption length is not 4.72e-9 cm-1!" << std::endl;
        return 22;
    }

    std::cout << std::endl;

    TRestAxionBufferGas mixture("bufferGases.rml", "HeNeMixture");

    std::cout << "Gas name : " << mixture.GetName() << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Number of gases : " << mixture.GetNumberOfGases() << std::endl;
    std::cout << "Photon mass at 3keV : " << mixture.GetPhotonMass(3) << " eV" << std::endl;
    std::cout << "Absorption lenght (cm-1) at 3keV : " << mixture.GetPhotonAbsorptionLength(3) << std::endl;

    if (std::round(1.e3 * mixture.GetPhotonMass(3)) != 651) {
        std::cout << "Photon mass is not 0.651 eV!" << std::endl;
        return 31;
    }

    if (std::round(1.e3 * mixture.GetPhotonAbsorptionLength(3)) != 397) {
        std::cout << "Absorption length is not 0.397 cm-1!" << std::endl;
        return 32;
    }

    std::cout << endl;
    std::cout << "--> Buffer gas evaluation succeeded!" << std::endl;
    std::cout << endl;

    //// Retrieving the magnetic field profile for a dummy trajectory.
    //// We force it to be strongly un-aligned with z

    TVector3 position1(-300, 600, -6000);
    TVector3 position2(600, -600, 6000);

    TVector3 direction = position2 - position1;
    direction.Unit().Print();

    std::vector<TVector3> v = magneticField.GetFieldBoundaries(position1, direction.Unit());
    if (v.size() != 2) {
        std::cout << "The track does not traverse the magnetic volume" << std::endl;
        return 100;
    }

    std::cout << "In position. X: " << v[0].X() << " Y: " << v[0].Y() << " Z: " << v[0].Z() << std::endl;
    std::cout << "Out position. X: " << v[1].X() << " Y: " << v[1].Y() << " Z: " << v[1].Z() << std::endl;

    TVector3 in = v[0];
    TVector3 out = v[1];

    if (std::round(v[0].X()) != -101 || std::round(v[0].Y()) != 335 || std::round(v[0].Z()) != -3350) {
        std::cout << "Wrong boundary determination. The input position is wrong!" << std::endl;
        return 101;
    }

    if (std::round(v[1].X()) != 293 || std::round(v[1].Y()) != -191 || std::round(v[1].Z()) != 1910) {
        std::cout << "Wrong boundary determination. The output position is wrong!" << std::endl;
        return 102;
    }

    Double_t dl = 50;  // 5cm
    std::vector<Double_t> bProfile = magneticField.GetTransversalComponentAlongPath(v[0], v[1], dl);

    //   for (auto& x : bProfile) std::cout << x << " ";
    //   std::cout << std::endl;

    std::cout << "Field elements (5cm step) : " << bProfile.size() << std::endl;
    if (bProfile.size() != 107) {
        std::cout << "Number of field elements is not the right value!" << std::endl;
        return 103;
    }

    std::cout << "Coherence length : " << (bProfile.size() - 1) * 50 << " mm" << std::endl;

    Double_t fieldAverage = std::accumulate(bProfile.begin(), bProfile.end(), 0.0) / bProfile.size();

    std::cout << "Field average (T): " << fieldAverage << std::endl;

    if (std::round(fieldAverage * 1000.) != 1901) {
        std::cout << "Field average is not 1.901 T!" << std::endl;
        return 104;
    }

    std::cout << endl;
    std::cout << "--> Track boundaries and field along path evaluation succeeded!" << std::endl;
    std::cout << endl;

    TRestAxionField axionField;
    axionField.AssignBufferGas(&helium);

    std::cout << "BL(8.8T, 9.26m) : " << axionField.BL(8.8, 9260) << std::endl;
    std::cout << "0.5 x (BL)^2 : " << axionField.BLHalfSquared(8.8, 9260) << std::endl;

    if (std::round(axionField.BL(8.8, 9260) * 1e10) != 81) {
        std::cout << "BL calculation is not 8.1e-9!" << std::endl;
        return 201;
    }

    if (std::round(axionField.BLHalfSquared(8.8, 9260) * 1e20) != 1627) {
        std::cout << "(BL)^2/2 calculation is not 1.627e-17!" << std::endl;
        return 202;
    }

    Double_t B = 8.8;     // T
    Double_t L = 9260.;   // mm
    Double_t ma = 1.e-2;  // eV
    Double_t Ea = 3.0;    // keV

    axionField.SetMagneticField(B);
    axionField.SetCoherenceLength(L);
    std::cout << "Probability (CAST. Helium. ma:10meV - Ea:3keV): "
              << axionField.GammaTransmissionProbability(Ea, ma) << std::endl;

    if (std::round(axionField.GammaTransmissionProbability(Ea, ma) * 1.e20) != 1547) {
        std::cout << "CAST Probability. Wrong axion-photon conversion probability" << std::endl;
        return 301;
    }

    Double_t prob = axionField.GammaTransmissionProbability(bProfile, 50, Ea, ma);
    std::cout << "Probability BabyIAXO : " << prob << std::endl;

    if (std::round(prob * 1.e22) != 2482) {
        std::cout << "BabyIAXO Probability. Wrong axion-photon conversion probability" << std::endl;
        return 401;
    }

    std::cout << endl;
    std::cout << "--> Axion-photon probability calculations succeeded!" << std::endl;
    std::cout << endl;

    return 0;
}
