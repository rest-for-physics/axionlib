
Transmission.C will do the calculation using IAXO like signal.

For BabyIAXO:

```
.L TransmissionSensitivity.C
TransmissionSensitivity( "IAXO.rml", "BabyIAXOSignal", 1.5, "limits/BabyIAXO.txt");
```

For IAXO:

```
.L TransmissionSensitivity.C
TransmissionSensitivity( "IAXO.rml", "IAXOSignal", 1.5, "limits/IAXO.txt");
```

For AMELIE-CMS:

```
.L TransmissionSensitivity.C
TransmissionSensitivity( "CMS.rml", "cmsSignal", 1.5, "limits/CMS.txt");
```
