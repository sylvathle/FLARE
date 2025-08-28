import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import sys

def calc_stats(serie):
  av = np.mean(serie)
  if av==0: return np.nan,np.nan,np.nan


  below_av = [(av-x)**2 for x in serie if x < av]
  above_av = [(av-x)**2 for x in serie if x >= av]
  
  low_std = np.sqrt(sum(below_av)/len(below_av))
  up_std = np.sqrt(sum(above_av)/len(above_av))
  
  return av, low_std, up_std

def getGCRFlux(gcr_path,listE,listdeltaE,E_cutoff=0):
  df_gcr = pd.read_csv(gcr_path)

  #df_gcr = pd.read_csv(gcr_path)
  # Data given in #proton / [MeV m2 s sr]
  #df_gcr = df_gcr[["Energy (MeV/n)"]]

  df_gcr.set_index("Energy (MeV/n)",inplace=True)
  for particle in df_gcr.columns:
    # Convert to 1 day of GCR / m2
    df_gcr[particle] = df_gcr[particle]*60 # Now in #proton / [MeV m2 min sr]


  df_gcr = df_gcr[df_gcr.index>1]

  dict_gcr =  {"E":listE,"deltaE":listdeltaE}

  for particle in df_gcr.columns:
    # Create an interpolation function for the dose values
    f2d_interp = interp1d(df_gcr.index, df_gcr[particle], kind='linear', fill_value="extrapolate")
    flux_ebin = f2d_interp(listE)
    for i in range(len(listE)):
      if flux_ebin[i]<0:
        flux_ebin[i] = np.nan

    dict_gcr[particle] =  flux_ebin

  df_gcr_resampled = pd.DataFrame(dict_gcr)
  df_gcr_resampled.set_index("E",inplace=True)

  for particle in df_gcr.columns:
    df_gcr_resampled[particle] = df_gcr_resampled[particle]*df_gcr_resampled["deltaE"]*4*np.pi/1e4 # #proton/[cm2.min]

  for particle in df_gcr_resampled.columns:
    for e in df_gcr_resampled.index:
      if e<E_cutoff:
        df_gcr_resampled.loc[e,particle] = 0
    #df_gcr_resampled[particle][df_gcr_resampled.index < cutoff] = 0

  df_gcr_resampled.drop("deltaE",axis=1,inplace=True)
  #print (df_gcr_resampled)
  return df_gcr_resampled

def gcrDose(df_gcr,df_coeffs):
  list_particles = df_coeffs["particle"].unique().tolist()

  list_sample = df_coeffs["i_sample"].unique().tolist()
  list_vals_AD = [0 for i in list_sample]
  list_vals_DE = [0 for i in list_sample]
  list_vals_EDE = [0 for i in list_sample]
  for particle in list_particles:
    df_coeff_particle = df_coeffs[df_coeffs["particle"]==particle]
    df_coeff_particle = df_coeff_particle.merge(df_gcr,left_index=True, right_index=True, how="outer")
    for isample in list_sample:
      df_coeff_sample = df_coeff_particle[df_coeff_particle["i_sample"]==isample]
      list_vals_AD[isample] += float((df_gcr[particle]*df_coeff_sample["AD"]).sum())
      list_vals_DE[isample] += float((df_gcr[particle]*df_coeff_sample["DE"]).sum())
      if "EDE" in df_coeff_sample.columns: list_vals_EDE[isample] += float((df_gcr[particle]*df_coeff_sample["EDE"]).sum())
  ad,ad_b,ad_t = calc_stats(list_vals_AD)
  de,de_b,de_t = calc_stats(list_vals_DE)
  ede,ede_b,ede_t = calc_stats(list_vals_EDE)
  return ad,ad_b,ad_t,de,de_b,de_t,ede,ede_b,ede_t


#listE = [10,100,1000,10000]
#listdeltaE = [90,900,9000,90000]
#df_gcr = getGCRFlux("input/GCR-spectra/gcr_1998-01-01.csv",listE,listdeltaE)
