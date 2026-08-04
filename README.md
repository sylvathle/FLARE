# FLARE

Code for Fluence-to-Dose calculation of ICRP145 shielded with an aluminum shell.
Includes ICRP145 re-exported to be simulated with/without the FLARE vest.

Contact point: sylvain.blunier@gmail.com

## Folder structure

1- *G4*
The most interesting is likely G4/, containing the C++ source for the scene, including the human phantom.
How to use it is explained in the following section.

2- *ICRP145_BlenderToG4*
ICRP145_BlenderToG4 contains the Python code to be executed inside the Blender environment with the ICRP145 + vest CAD model previously loaded. The model can be downloaded from [Zenodo](10.5281/zenodo.16780309)(https://zenodo.org/api/records/19456639/draft/files/CAD_ICRP145withVest.zip/content)
These scripts export the ICRP145 to the poly file format, following the volume structure defined in list_objs_names.py.
This code has not been tested outside its original environment. To be used, you first need the ICRP145 + vest, and some paths need to be adjusted.
To use it without the vest, the CAD model provided by ICRP can be loaded, and a few modifications are needed to remove vest-related references from the code.

3- *SEP_spectra*
This folder contains Python scripts that retrieve the Band parameters of SEPs from different sources.

4- *pyflare*
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
The .body files (list of tetrahedra that constitute each organ) of the ICRP145+vest are available on [Zenodo](10.5281/zenodo.16780309) under the compressed folder 	
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

### Macros commands

The macros should work with the usual macro commands, plus some that have been defined for simulating a human phantom in an aluminum shell.


Define the location of the CSV file that specifies the organs/parts to be included in the scene. 
```
/SIM/scoring/csvBodies ../scene/organsAndVest.csv
```
There is an example, organsAndVest.csv, in the ICRP145_vest (or the scene folder) that can be targeted directly. 
In the current state of the project,../scene is hard-coded, so the files containing the tetrahedra of each organ (.body) should be in this folder. You must then create a folder G4/scene and copy all the files of ICRP145_vest available in the [dataset](https://zenodo.org/records/19456639) on Zenodo.
If the original phantom from the Geant4 examples is used, there is no need to call this line (it will be ignored); however, the environment variable PHANTOM_PATH needs to be defined and point to the folder examples/advanced/ICRP145_HumanPhantoms/build/ICRP145data/ of your Geant4 installation.


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

Sphere from which the source will be generated.
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

A file *run_*scenario*_*thickness*_*particle*.mac* will be created.


### First run
For a first quick test (assuming all the data about phantoms is correctly setup as described above), you can run:
```
cd ../macros
python3 genRunMacro.py B2G-vest 100 4 proton 1 5
cd ../build
./sim macro.mac
```
This creates a macro placing the phantom with the vest in a 4 mm shell, prepare to shoot 100 proton with an energy that is logaritmically selected between 1 MeV and 100 GeV. 

If everything runs fine, the folders "results/B2G-vest_4_proton should be created in build, where CSV files containing the results are created.


### Resume simulation

This section assumes the Geant4 simulations were run using a specific folder-naming convention for the CSV files as defined in the macro files: (results/scenario_thickness_particle)

```
cd G4/resume_sim
python3 toCSV.py path_to_results
```

A file *dose_data.csv* is created with the following columns:
- scenario: as defined above (ICRP-naked, B2G-vest...)
- thick: Thickness of the spherical spacecraft given in mm, density is taken to 2.710 g.cm-2
- particle: proton, deuteron, Fe56...
- organId: identify the organ that was hit; detailed organId: https://zenodo.org/api/records/16780309/draft/files/organsInfo_IFHP.csv/content
- eBin: primary energy bin number
- DE: Averaged dose equivalent per particle in the corresponding bin (mSv/particle)
- DE_b and DE_t: asymmetric bottom and top standard deviation of DE for the entire simulation
- AD: Averaged absorbed dose per particle in the corresponding bin (mGy/particle)
- AD_b and AD_t: asymmetric bottom and top standard deviation of AD for the entire simulation

