import pandas as pd
import numpy as np
import datetime as dt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

import band
import physics
import traceback


class Sepem:

  def __init__(self,sample=None):
    self.df_fpio = pd.read_csv("input/SEPEM_RDS_v3/SEPEM_RDS_v3-2/RDS3.2_FPDO_FPIO.csv")
    self.df_fpio["epoch"] = pd.to_datetime(self.df_fpio["epoch"])
    self.df_fpio.set_index("epoch",inplace=True)

    self.fpio_cols = ["fpio_"+str(i) for i in range(1,8)]
    self.fpio_min = [5,10,30,50,100,300,500]
    self.log_fpio_min = [np.log(e) for e in self.fpio_min]

    #Flux are given in cm-2.sr-1.s-1 each 5 minutes
    # We convert it in cm-2.min-1
    #if sample=="H": self.df_fpio = self.df.resample("H").sum()
    #self.df_fpio = self.df.copy()
    self.df_fpio = self.df_fpio.loc[~(self.df_fpio==0).all(axis=1)]
    ##### Need to make sure of this 2 factor -> nope it's 4pi
    self.df_fpio = 4.0*np.pi*60*self.df_fpio[self.fpio_cols]

    #print (self.df_fpio)

    #self.df_fpio = self.df_fpio[self.df_fpio.index>dt.datetime(2017,1,10,0,0,0,0)]


  def getSpectrum(self,datemin,datemax,listE):

    dict_spe =  {"E":listE[:-1]}

    listR = [physics.rigidity(e,physics.mass_proton,1) for e in listE]

    last_popt = [600*2*np.pi, 1.0e-01, 5.36086367e+00, 8.83179720e-03]

    for date,row in self.df_fpio.iterrows():

      if date>datemax or date<datemin: continue

      dat_spectrum = row.tolist()
      date_str = date.strftime("%Y%m%d-%H%M")

      log_dat_spectrum = []
      R_fit = []
      for i in range(len(dat_spectrum)):
        if dat_spectrum[i]>0: 
          log_dat_spectrum.append(np.log10(dat_spectrum[i]))
          R_fit.append(physics.rigidity(self.fpio_min[i],physics.mass_proton,1))
      
      #log_E = []
      #for i in range(len(listE)):
      #  log_E.append(np.log(listE[i]))

      try:
        band_popt, pcov = curve_fit(band.logBand, R_fit, log_dat_spectrum, p0=last_popt, bounds=([0.0,-10.0,0.0,0.0], [np.inf,10.0,100.0,10.0]), maxfev=100000000)
        #last_popt = band_popt
      #  print (date,band_popt[0])
        #flux_ebin = band.lBand(listR,*band_popt)
        #listRp1 = [physics.rigidity(e,physics.mass_proton,1) for e in listE]
        flux_ebin = band.getBandSpectrum(listR,band_popt[0],band_popt[1],band_popt[2],band_popt[3])
      #  #print ("Succeed band")

      except:
      #if True:
        #traceback.print_exc()
      #finally:
        f2d_interp = interp1d(self.fpio_min, dat_spectrum, kind='linear', fill_value="extrapolate")
        flux_ebin = f2d_interp(listE[:-1])

      for i in range(len(listE[:-1])):
        if flux_ebin[i]<0: flux_ebin[i] = np.nan
        #else: flux_ebin[i] = 10**flux_ebin[i]
        

      dict_spe[date_str] = flux_ebin  

    df_spe_resampled = pd.DataFrame(dict_spe)
    df_spe_resampled.set_index("E",inplace=True)

    return df_spe_resampled



