
```
TRestSensitivity sens("BabyIAXO.rml", "VacuumPhase");
sens.GetExperiment(0)->GetSignal()->RegenerateParametricNodes(0.001,10,1.02,true);
sens.GenerateCurve()
sens.ExportCurve("vacuum.txt", 0 )
```

