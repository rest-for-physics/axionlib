
Int_t Validation() {
    TRestSensitivity sens("DummyIAXO.rml", "VacuumPhase");
    sens.GenerateCurve();
    sens.GetExperiment(0)->GetSignal()->PrintMetadata();
    sens.GetExperiment(0)->GetBackground()->PrintMetadata();

    if (sens.GetExperiment(0)->GetSignal()->GetParameterizationNodes().size() != 2) {
        std::cout << "Error! The number of parameterization nodes is not two!" << std::endl;
        return 1;
    }

    if (sens.GetCurve().size() != 2) {
        std::cout << "Error! The generated sensitivity curve should have 2 nodes!" << std::endl;
        return 2;
    }

    if (sens.GetCurve()[0] < 0.00071575 || sens.GetCurve()[0] > 0.00071576) {
        std::cout << "error! the first point on the sensitivity curve should be 0.000715751!" << std::endl;
        std::cout << "present value: " << sens.GetCurve()[0] << std::endl;
        return 3;
    }

    if (sens.GetCurve()[1] < 0.6135 || sens.GetCurve()[1] > 0.6136) {
        std::cout << "error! the second point on the sensitivity curve should be 0.613593!" << std::endl;
        std::cout << "present value: " << sens.GetCurve()[1] << std::endl;
        return 4;
    }
    return 0;
}
