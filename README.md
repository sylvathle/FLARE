# FLARE

Code for Fluence-to-Dose calculation of ICRP145 shielded with aluminium shell.

Contact point: sylvain.blunier@gmail.com


## Running Geant4 simulations

### ICRP145
Data of the ICRP 145 must be set up to run the code, the global environment variable _PHANTOM\_PATH_ is needed.
A simple way is to run the advanced example of the ICRP145\_HumanPhantoms.

Then
```
export PHANTOM_PATH=geant4dir/examples/advanced/ICRP145_HumanPhantoms/build/ICRP145data/
```

### Install
After cloning the project, the usual Geant4 C++ commands must be run:

```
cd FLARE
mkdir build
cd build
cmake ..
make -j8
```

### Generate macros

Macros are generated with a python3 script in the macros folder:

```
cd ../macros
python3 genRunMacros.py scenario nsim thickness
```

- *scenario*: string, takes only the value 'ICRP-naked' for now but it will evolve in a future version of the code.
- *nsim*: integer, corresponds de the number protons to be simulated
- *thickness*: float, corresponds to the thickness of the aluminum shell in millimeters

A file *run_*scenario*_*thickness*.mac* will be created.

example:
```
python3 genRunMacros.py ICRP-naked 100 100
```

### First run
The code can be run using the test.mac macros in the inputs folder

```
cd ../build
./sim ../macros/run_scenario_thickness.mac
```
If everything runs fine a folder called "results/ICRP145/scenario_thickness/" should be created where are created csv files containing the results.
Each run creates 5 files with pattern YYYYMMDD-HHmmss-XXXXXXX\_nt\_Dose.csv, where XXXXXXX is a random string of numbers, the file contains the absorbed dose (proton_AD) and Dose equivalent (proton_DE) for each organ, and each primary energy bin.

### Macros commands

The macros should work with the usual macro commands plus some that have been defined for the simulation of a human phantom in an alumium shell (see test.mac in the inputs folder as example).


```
/SIM/scoring/phantom ICRP145
```
Define the human phantom to be simulated: HP can take the value ICRP145 only for now.

```
/SIM/scoring/putModule 10
```
Indicating the thickness of the alumium shell in millimeters.


```
/SIM/scoring/sampleSize Nsample
```
Nsample stands for the number of particle that should be generated before storing the results, if Nsample if bigger Nsimulated with the /run/beamOn command then no data will be stored.
If Nsample = 10 and 100 particles where simulated, then 10 runs will be saved.

```
/SIM/scoring/resDir ../results/ICRP145/ICRP-naked_10
```
Indicate folder where the data scored will be stored.


```
/SIM/scoring/radbeam radius mm
```
radius stands for the radius of the hemisphere from where primaries are generated. For logical reason the user should make sure it is a bigger number than the radius of the last layer of the dome. (innerRadius plus the 4 thicknesses of the 4 layers).


## Analysis

This sections assumes the geant4 simulations ran following a specific way of naming the folders containing the csv files.
The data has been resumes and shared at https://doi.org/10.5281/zenodo.16036474, they can be downloaded and moved to a data folder inside the analysis folder to run the analysis.

A folder called figures must be created at the same level than the analysis folder (outside the analysis folder).

```
cd FLARE
cd analysis
```

The file https://zenodo.org/records/14622711/files/organsInfo.csv is needed and must be placed in the analysis folder.

```
wget https://zenodo.org/records/14622711/files/organsInfo.csv
```

The analysis related to doses calculations is run with:

```
python3 plotAnalysis.py
```

A file *dose_data.csv* is created with the following columns:

- scenario: as defined above (ICRP-naked)
- thick: Thickness of the spherical spacecraft given in mm, density is taken to 2.710 g.cm-2
- eBin: primary energy bin number
- organId: organId as defined by ICRP: https://www.icrp.org/publication.asp?id=ICRP%20Publication%20145
- group: Organ group following organsInfo.csv
- WT: weithing tissue factor
- mass[g]: mass of organ
- proton_N: number of simulated proton for the corresponding eBin
- proton_AD,proton_AD_b,proton_AD_t: Absorbed dose (AD) (mGy.cm2/proton) multiplied by the surface of the simulated source per simulated proton in the corresponding energy bin and organ, the _b and _t subscripts indicate the bottom and top standard deviations
- proton_DE,proton_DE_b,proton_DE_t: Dose equivalent (DE) (mSv.cm2/proton) multiplied by the surface of the simulated source per simulated proton in the corresponding energy bin and organ, the _b and _t subscripts indicate the bottom and top standard deviations.

