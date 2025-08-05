Band function fits for solar storms from different sources

# Goes data: 

SPE band function fits are prepared in the *FLARE/SEP_spectra* folder:
```
cd FLARE/SEP_spectra
```

Requires to have directories:
```
mkdir -p data & mkdir -p data/Goes
```

and file *SPEs_cycle-25.csv* from https://doi.org/10.5281/zenodo.16036474 must be placed in *data/Goes/*

Then the SPE data, based on *SPEs_cycle-25.csv* information, from Goes satellites can be downloaded and fitted.

```
python3 fit_goes.py
```

A file *data/Goes/band_goes.csv* is created containing the Band function parameters for each SPE and each one of the goes satellite (16,17 and 18) where data was found.
