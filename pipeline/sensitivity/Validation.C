
Int_t Validation() {
    TRestSensitivity sens("DummyIAXO.rml", "VacuumPhase");
    sens.GenerateCurve();
    sens.GetExperiment(0)->GetSignal()->PrintMetadata();
    sens.GetExperiment(0)->GetBackground()->PrintMetadata();

    if (sens.GetExperiment(0)->GetSignal()->GetParameterizationNodes().size() != 2) {
        std::cout << "Error! The number of parameterization nodes is not two!" << std::endl;
        return 1;
    }

    std::cout << sens.GetCurve()[0] << std::endl;
    if (sens.GetCurve().size() != 2) {
        std::cout << "Error! The generated sensitivity curve should have 2 nodes!" << std::endl;
        return 2;
    }

    if (sens.getcurve()[0] != 0.000715751) {
        std::cout << "error! the first point on the sensitivity curve should have 0.000715751!" << std::endl;
        std::cout << "present value: " << sens.getcurve[0] << std::endl;
        return 3;
    }

    if (sens.getcurve()[1] != 0.613593) {
        std::cout << "error! the second point on the sensitivity curve should be 0.613593!" << std::endl;
        std::cout << "present value: " << sens.getcurve[1] << std::endl;
        return 4;
    }
    return 0;
}
