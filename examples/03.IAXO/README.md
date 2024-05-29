
## Available configuration files

- `BabyIAXO.rml`: It will serve to reproduce the sensitivityof BabyIAXO data taking phase. We include vacuum phase and a combination of vacuum and 73 density settings that extend the axion search up to 0.25 eV.
- `IAXO.rml`:

This example is used to reproduce the different scenarios considered inside the "E. Armengaud et al JCAP06(2019)047, Physics Potential of the International Axion Observatory". Those scenarios are summarized in the following table, where some values, concerning detector efficiency, $\epsilon_d$, have been replaced by just window detector efficiency, $\epsilon_w$, since the detailed detector response is considered and constructed using [restG4](https://github.com/rest-for-physics/restG4/tree/master), see example [14.DetectorResponse](https://github.com/rest-for-physics/restG4/tree/master/examples/14.DetectorResponse).

Parameter | Units   | BabyIAXO | IAXO baseline | IAXO upgraded |
  :---:   |  :---:  |  :---:   |     :---:     |      :---:    |
B         |   $T$     |   ~2     |     ~2.5      |      ~3.5     |
L         |   $m$     |   10     |      20       |       22      |
A         |  $m^2$     |  0.77    |      2.3      |       3.9     |
 - |  -  |  -  |  -  | -  |
b            | $keV^{-1}cm^{-2}s^{-1}$ | $1\times10^{-7}$ | $1\times 10^{-8}$ | $1\times10^{-9}$ |
$\epsilon_w$ |         |   0.85      |      0.95     |     0.95      |
$\epsilon_o$ |         |   0.35      |      0.7      |     0.7       |
a            | $cm^2$  |  2 x 0.3    |    8 x 0.15   |   8 x 0.15    |
$\epsilon_t$ |         |   0.5       |      0.5      |     0.5       |
t            | $year$    |   3+3       |      6+6      |    10+10      |


## Vacuum sensitivity curve generation

```
restRoot
[0] TRestSensitivity sens("BabyIAXO.rml", "VacuumPhase");
[1] sens.GenerateCurve()
[2] sens.ExportCurve("output/BabyIAXO_vacuum.txt", 0 )
```
## Combined vacuum and gas phase sensitivity curve generation

We need first to pre-generate the signals for the different density settings.

```
restRoot
[0] .L GenerateSignalComponents.C
[1] GenerateSignalComponents( "BabyIAXO.rml", "GasSignal" );
```

Then, inside `TRestSensitivity` we define an experiment list, with a
column which defines the signal to be used, with common background
and common exposure time.

```
[0] TRestSensitivity sens("BabyIAXO.rml", "CombinedPhase");
[1] sens.GenerateCurve()
[2] sens.ExportCurve("output/BabyIAXO_vacuum.txt", 0 )
```

**Hints**
- A detailed x-ray detector response is included in the signal calculation. The pre-calculated response matrix will be convoluted with the axion energy spectrum resulting in a calculation that is 150 times (number of reponse matrix bins) more expensive. If removed, calculation should achieve better computational timing.
- The number of parameter nodes, mass values at which the signal is calculated is of the order of 500 points, which will lead to a HD curve. Reducing the number of points by increasing the parameter `stepParameterValue` will also reduce the computational cost.