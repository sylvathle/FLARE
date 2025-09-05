#from scipy.io import netcdf
#import netCDF4
import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt
#import datetime as dt
#import numpy as np
import pandas as pd
import sys,os

#import warnings
import pandas as pd

import band
import physics

def get_stats(group,vars_list):

    dict_out = {}
    #print (group)

    #for p in list_particles:
    if True:
      
      for q in vars_list:
      
        dict_out[q] = group[q].sum() / len(group)

        low_de_group = group[group[q]<dict_out[q]]
        if len(low_de_group)==0: dict_out[q+"_b"] = np.nan
        else: dict_out[q+"_b"] = np.sqrt(((low_de_group[q] - dict_out[q]) ** 2).sum() / len(low_de_group))

        up_de_group = group[group[q]>=dict_out[q]]
        if len(up_de_group)==0: dict_out[q+"_t"] = np.nan
        else: dict_out[q+"_t"] = np.sqrt(((up_de_group[q] - dict_out[q]) ** 2).sum() / len(up_de_group))

    return pd.Series(dict_out)

class DoseInfo:

  def __init__(self,sampled=True):
    # Binning of primaries (kinetic energies are assumed in MeV/n)
    self.minE = 1
    self.maxE = 5
    self.Nbins = 100

    f_organ_grouped = "input/df_organ_grouped.csv"
    f_body_doses = "input/df_body_doses.csv"

    if (os.path.isfile(f_organ_grouped) & os.path.isfile(f_body_doses)):
      self.df_organ_grouped = pd.read_csv(f_organ_grouped)
      self.df_organ_grouped.set_index("E",inplace=True)

      self.df_body_doses = pd.read_csv(f_body_doses)
      self.listE = sorted(self.df_body_doses["E"].unique().tolist())
      self.listdeltaE = sorted(self.df_body_doses["deltaE"].unique().tolist())
      self.df_body_doses.drop("deltaE",axis=1,inplace=True)
      self.df_body_doses.set_index("E",inplace=True)

      self.list_sample = self.df_body_doses["i_sample"].unique().tolist()
      self.list_scenario = self.df_body_doses["scenario"].unique().tolist()
      self.list_thick = self.df_body_doses["thick"].unique().tolist()
      self.list_particle = self.df_body_doses["particle"].unique().tolist()

      self.list_organs = self.df_organ_grouped["group"].unique().tolist()

      self.listEp1 = self.listE + [10**self.maxE]
      self.listRp1 = [physics.rigidity(e,physics.mass_proton,1) for e in self.listEp1]

    else:
      self.loadDoses(sampled)

  def loadDoses(self,sampled):

    file_organICRP = 'input/organsInfo_ICRP.csv'
    self.df_organs_ICRP = pd.read_csv(file_organICRP)

    file_organB2G = 'input/organsInfo_b2g.csv'
    self.df_organs_b2g= pd.read_csv(file_organB2G)

    if sampled:
      file_path = 'input/doses_samples_v2.csv'
      list_index = ["i_sample","scenario","thick","particle","E","deltaE","group","WT"] # organId removed
    else: 
      file_path = 'input/data_dose_v2.csv'
      list_index = ["scenario","thick","particle","E","deltaE","group","WT"] # organId removed


    self.df_dat = pd.read_csv(file_path)
    self.df_dat["E"] = (self.df_dat["eBin"])/self.Nbins*(self.maxE-self.minE)+self.minE
    self.df_dat["E"] = 10**self.df_dat["E"]
    self.df_dat["deltaE"] = 10**((self.df_dat["eBin"]+1)/self.Nbins*(self.maxE-self.minE)+self.minE)-self.df_dat["E"]

    self.df_dat["organId"] = self.df_dat["organId"].astype(int)
    self.df_dat = self.df_dat[self.df_dat["organId"]<15000]

    self.list_scenario = self.df_dat["scenario"].unique().tolist()

    self.df_dat_merged = pd.DataFrame()

    for scenario in self.list_scenario:
      df_scenario = self.df_dat[(self.df_dat["scenario"]==scenario)]
      if "ICRP" in scenario: df_scenario = df_scenario.merge(self.df_organs_ICRP[['organId', 'group', 'WT','mass[g]']], on='organId', how='left')
      elif "B2G" in scenario: df_scenario = df_scenario.merge(self.df_organs_b2g[['organId', 'group', 'WT','mass[g]']], on='organId', how='left')
      self.df_dat_merged = pd.concat([self.df_dat_merged,df_scenario])

    self.df_dat = self.df_dat_merged

    list_vals = []

    for c in self.df_dat.columns:
      if "AD" in c or "DE" in c:
        list_vals.append(c)
      
        # Important line correcting the m2 unit for the sphere of generation of protons to cm2
        # It also cancel the angular distribution as it is already multiplied in the band fitted fluxes
        self.df_dat[c] *= 1e4/(2*np.pi)/4.0
        # Correct multiplication done in simulations by energy bin size, not needed
        self.df_dat[c] /= self.df_dat["deltaE"] # cm2.mSv/proton
      if "_N" in c:
        self.df_dat.drop(columns=c,inplace=True)

    self.listE = self.df_dat["E"].unique().tolist()
    self.listE.sort()
    self.listEp1 = self.listE + [10**self.maxE]
    self.listRp1 = [physics.rigidity(e,physics.mass_proton,1) for e in self.listEp1]
    self.listdeltaE = self.df_dat["deltaE"].unique().tolist()
    self.listdeltaE.sort()

    self.listR = [physics.rigidity(e,physics.mass_proton,1) for e in self.listE]
    #self.df_dat.drop("deltaE",axis=1,inplace=True)
    #list_index.remove("deltaE")

    #listEp1 = listE + [10**((df_dat["eBin"].iloc[-1]+1)/Nbins*(maxE-minE)+minE)]
    #listRp1 = listR + [rigidity(listEp1[-1],mass_proton,1)]


    #self.list_thick = sorted(self.df_dat["thick"].unique().tolist())

    # group by organId
    self.df_organ_grouped = self.df_dat.copy()
    self.df_organ_grouped.drop(columns='organId',inplace=True)
    for c in list_vals:
      self.df_organ_grouped[c] = self.df_organ_grouped[c]*self.df_organ_grouped['mass[g]']
    self.df_organ_grouped = self.df_organ_grouped.groupby(list_index,as_index=False).sum()
    for c in list_vals:
      self.df_organ_grouped[c] = self.df_organ_grouped[c]/self.df_organ_grouped['mass[g]']

    self.df_organ_grouped.sort_values(by=list_index,inplace=True,ascending=True)

    self.body_mass = self.df_organs_ICRP["mass[g]"].sum()

    pd.options.mode.copy_on_write = True
    self.list_organs = self.df_organ_grouped["group"].unique().tolist()

    self.df_organ_grouped[["i_sample","scenario","thick","particle","group","E","DE","AD"]].to_csv("input/df_organ_grouped.csv",index=False)

    list_stat = ['','_b','_t']
    list_stat = ['']

   

    self.df_body_doses = self.df_organ_grouped.copy()
    for stat in list_stat:
      self.df_body_doses["EDE"+stat] = self.df_body_doses["DE"+stat]*self.df_body_doses["WT"]
      self.df_body_doses["DE"+stat] = self.df_body_doses["DE"+stat]*self.df_body_doses["mass[g]"] / self.body_mass
      self.df_body_doses["AD"+stat] = self.df_body_doses["AD"+stat]*self.df_body_doses["mass[g]"] / self.body_mass

    self.df_body_doses.drop(columns=["mass[g]","WT"],inplace=True)
    list_index.remove("WT")
    list_index.remove("group")
    self.df_body_doses = self.df_body_doses.groupby(list_index,as_index=False).sum()

    #list_index.remove("i_sample")

    #print (self.df_body_doses)
    #self.df_body_doses = self.df_body_doses.groupby(by=list_index,as_index=False).apply(lambda x: pd.concat([get_stats(x,list_vals)], axis=0)).copy()
    #print (self.df_body_doses)

    self.df_body_doses = self.df_body_doses[["i_sample","scenario","thick","particle","deltaE","E","AD","DE","EDE"]]
    #df_organ_grouped = df_organ_grouped[["scenario","thick","group","particle","E","deltaE","AD","AD_b","AD_t","DE","DE_b","DE_t"]]
    self.df_body_doses.sort_values(by=["i_sample","scenario","thick","particle","E"],inplace=True,ascending=True)
    self.df_body_doses.set_index("E",inplace=True)

    self.df_body_doses.to_csv("input/df_body_doses.csv",index=True)


    self.list_sample = self.df_body_doses["i_sample"].unique().tolist()
    self.list_scenario = self.df_body_doses["scenario"].unique().tolist()
    self.list_thick = self.df_body_doses["thick"].unique().tolist()
    self.list_particle = self.df_body_doses["particle"].unique().tolist()

  def getFlux2DoseCoeffs(self,scenario,thick,organ="all"):
    if organ=="all":
      df_spec = self.df_body_doses[(self.df_body_doses["scenario"]==scenario) & (self.df_body_doses["thick"]==thick)]
    else:
      df_spec = self.df_organ_grouped[(self.df_organ_grouped["scenario"]==scenario) & \
                               (self.df_organ_grouped["thick"]==thick) &\
                               (self.df_organ_grouped["group"]==organ)]
    return df_spec.sort_index()#[list_c]

