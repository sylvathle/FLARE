# FLARE

Code for Fluence-to-Dose calculation of ICRP145 shielded with an aluminium shell.
Includes ICRP145 re-exported to be simulated with/without the FLARE vest.

Contact point: sylvain.blunier@gmail.com

## Folder structure

1- G4
The most interesting is likely G4/, containing the C++ source for the scene, including the human phantom.
How to use it is explained in the following section.

2- ICRP145_BlenderToG4
ICRP145_BlenderToG4 contains the Python code to be executed inside the Blender environment with the ICRP145 + vest previously loaded. 
This series of codes exports ICRP145 to the poly file format, following the volume structure defined in list_objs_names.py.
This code has not been tested outside its original environment. To be used, you first need the ICRP145 + vest, and some paths need to be adjusted.
To use it without the vest, the CAD model provided by ICRP can be loaded, and little modifications are needed to remove references to the vest in the code.

3- SEP_spectra
This folder contains Python scripts that retrieve the Band parameters of SEPs from different sources.

4- pyflare
This folder contains Python scripts that convert Geant4-produced data into flux-to-dose coefficients and compute dose quantities from particle spectra.
It wasn't tested outside its original working environment.

## Running Geant4 simulations

### ICRP145

With the original ICRP145: 
Data from ICRP 145 must be set up to run the code; the global environment variable _PHANTOM\_PATH_ is required.
A simple way is to run the advanced example in ICRP145\_HumanPhantoms.

Then
```
export PHANTOM_PATH=geant4dir/examples/advanced/ICRP145_HumanPhantoms/build/ICRP145data/
```

With the newly exported ICRP145:
The information of the ICRP145+vest is available on [Zenodo](10.5281/zenodo.16780309) under the compressed folder 	
[ICRP145_vest](https://zenodo.org/api/records/16780309/draft/files/ICRP145_vest.tar.gz/content).
The folder should be uncompressed and copied into the G4/ folder as "scene", or update the route to its new location in the file G4/src/TETModelImport.cc.

### Install
The usual Geant4 C++ commands must be run:

```
cd FLARE/G4
mkdir build
cd build
cmake ..
make -j8
```

### First run
For a first test, you can use the test macro available in G4/: 
```
cd ../build
./sim test.mac
```
If everything runs fine, a folder called "results/ICRP145/scenario_thickness/" should be created, where CSV files containing the results are created.
Each run creates 5 files with pattern YYYYMMDD-HHmmss-XXXXXXX\_nt\_Dose.csv, where XXXXXXX is a random string of numbers. Each file contains the absorbed dose (proton_AD) and Dose equivalent (proton_DE) for each organ, and each primary energy bin.

### Macros commands

The macros should work with the usual macro commands, plus some that have been defined for simulating a human phantom in an aluminium shell.


Define the location of the CSV file that specifies the organs/parts to be included in the scene. 
```
/SIM/scoring/csvBodies ../scene/bodiesAndVest.csv
```
There is an example, bodiesAndVest.csv, in the ICRP145_vest (or the scene folder) that can be targeted directly.
If the original phantom is used, there is no need to call this line (will be ignored).


Define the human phantom to be simulated
```
/SIM/scoring/phantom HP
```
where HP can take the following values:
- *BDRTOG4*: the newly exported phantom
- *ICRP145*: orginal human phantom

In order to include the shell, use
```
/SIM/scoring/putModule 100
```
Indicating the thickness of the aluminum shell in millimeters.

Indicate the directory where the results should be stored:
```
/SIM/scoring/resDir dest_folder
```

Sphere from where the source will be generated.
```
/SIM/scoring/radbeam 1500 mm
```
For logical reasons, the user should ensure it is larger than the radius of the dome's last layer.

Particle to be generated  (H to Ni, e+/e-, pi+/pi-, deuteron, gamma, neutron)
```
/SIM/generate/particle neutron
```

Energy range to be simulated (with log scale distribution), logarithm base 10 of the minimum energy to be used in MeV (e.g. -1 would mean 0.1 MeV)
```
/SIM/generate/logminkE -9
/SIM/generate/logminkE 5
```

```
/run/beamOn 100
```
Number of particles to simulate.

```
/SIM/scoring/resDir ../results/ICRP145/ICRP-naked_10
```
Indicate the folder where the data scored will be stored.

A Python script has been implemented to generate macros:

```
cd ../macros
python3 genRunMacro.py scenario nsim thickness particle logEmin logEmax
```

- *scenario*: string, takes the value 'ICRP-naked' to use the original mesh-type ICRP145 human phantom, "B2G-vest" or "B2G-naked" to use the newly exported human phantom with or without the vest. "B2G" stands for "Blender to Geant4". 
- *nsim*: integer, corresponds to the number of particles to be simulated
- *thickness*: float, corresponds to the thickness of the aluminum shell in millimeters
- *particle*: particle to be used (H to Ni, e+/e-, pi+/pi-, deuteron, gamma, neutron)
- *logEmin*: logarithm base 10 of the minimum energy to be used in MeV (e.g., -1 would mean 0.1 MeV)
- *logEmax*: logarithm base 10 of the maximum energy to be used in MeV (e.g., 3 would mean 1 GeV)

A file *run_*scenario*_*thickness*.mac* will be created.

example:
```
python3 genRunMacro.py BDRTOG4 100 100 neutron -9 5
```


## Analysis

This section assumes the Geant4 simulations were run using a specific folder-naming convention for the CSV files.
The data has been resumed and shared at https://doi.org/10.5281/zenodo.16036474, they can be downloaded and moved to a data folder inside the analysis folder to run the analysis.

A folder named figures must be created at the same level as the analysis folder (outside it).

```
cd FLARE
cd analysis
```

The file https://zenodo.org/records/14622711/files/organsInfo.csv is needed and must be placed in the analysis folder.

```
wget https://zenodo.org/records/14622711/files/organsInfo.csv
```

The analysis related to dose calculations is run with:

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

