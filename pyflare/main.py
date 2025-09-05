#from scipy.io import netcdf
#import netCDF4
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import pandas as pd
import sys, os

import warnings
import pandas as pd

import band
import physics
import gcr
import DoseInfo as di
import sepem as sp

import time


def calc_stats(serie):
  av = np.mean(serie)
  if av==0: return np.nan,np.nan,np.nan


  below_av = [(av-x)**2 for x in serie if x < av]
  above_av = [(av-x)**2 for x in serie if x >= av]
  
  low_std = np.sqrt(sum(below_av)/len(below_av))
  up_std = np.sqrt(sum(above_av)/len(above_av))
  
  return av, low_std, up_std

def generateGCRDoses():
  dose = di.DoseInfo(True)

  if not os.path.exists("output"): os.makedirs("output")
  if not os.path.exists("output/GCR"): os.makedirs("output/GCR")
  if not os.path.exists("output/GCR/tserie"): os.makedirs("output/GCR/tserie")

  thick = 185
  d_increment = 183

  gcr_input_path = 'input/GCR-spectra/'
  gcr_output_path = 'output/GCR/tserie/'

  list_qqty = ["AD","DE","EDE"]
  list_stat = ["","_b","_t"]

  for scenario in dose.list_scenario:
    if "ICRP" in scenario: continue

    datemin = dt.datetime(1974,7,1,0,0,0,0)
    while datemin<dt.datetime(2025,8,1,0,0,0,0):
    #while datemin<dt.datetime(1974,7,13,0,0,0,0):
      datemax = datemin + dt.timedelta(days=d_increment)
      csv_filename = scenario+"_"+datemin.strftime("%Y%m%d")+"-"+datemax.strftime("%Y%m%d")+".csv"
      filepath = gcr_output_path+csv_filename

      if os.path.exists(filepath): 
        datemin = datemax
        continue
      with open(filepath, 'a'): os.utime(filepath, None)

      if "B2G" not in scenario: continue
      dict_gcr_doses = {"date":[],"organ":[]}

      for qqty in list_qqty:
        for stat in list_stat:
          dict_gcr_doses[qqty+stat] = []

      for filename in os.listdir(gcr_input_path):
        print(filename)
        # Extract date from filename (assuming format "gcr_%Y-%m-%d.csv")
        gcr_path_file = gcr_input_path+"/"+filename
        try:
          date_str = filename.replace("gcr_", "").replace(".csv", "")
          date = dt.datetime.strptime(date_str, "%Y-%m-%d")
          #print(f"  Extracted date: {date}")

        except ValueError:
          print(f"  Could not parse date from filename: {filename}")
        if date < datemin: continue
        if date > datemax: continue

        df_gcr = gcr.getGCRFlux(gcr_path_file, dose.listE, dose.listdeltaE ,0)

        for organ in dose.list_organs + ["all"]:
          #print (scenario,date,organ)

          dict_gcr_doses["date"].append(date)
          dict_gcr_doses["organ"].append(organ)

          df_f2d = dose.getFlux2DoseCoeffs(scenario,thick,organ)
          ad,ad_b,ad_t,de,de_b,de_t,ede,ede_b,ede_t = gcr.gcrDose(df_gcr,df_f2d)
          dict_gcr_doses["AD"].append(ad)
          dict_gcr_doses["AD_b"].append(ad_b)
          dict_gcr_doses["AD_t"].append(ad_t)
          dict_gcr_doses["DE"].append(de)
          dict_gcr_doses["DE_b"].append(de_b)
          dict_gcr_doses["DE_t"].append(de_t)
          dict_gcr_doses["EDE"].append(ede)
          dict_gcr_doses["EDE_b"].append(ede_b)
          dict_gcr_doses["EDE_t"].append(ede_t)

      df_gcr_doses = pd.DataFrame(dict_gcr_doses)
      df_gcr_doses.set_index("date",inplace=True)
      #for col in df_gcr_doses.columns:
      #for qqty in list_qqty:
      #  for stat in list_stat:
      #    col = qqty+stat
      #    df_gcr_doses[col] = round(df_gcr_doses[col],5) #/(24*60)*1000
      df_gcr_doses.sort_index(inplace=True)

      df_gcr_doses.to_csv(filepath)
      datemin = datemax

def generateSPEDoses(): 

  dose = di.DoseInfo(True)
  ba = band.Band()

  if not os.path.exists("output"): os.makedirs("output")
  if not os.path.exists("output/SPE"): os.makedirs("output/SPE")

  thick = 185

  #gcr_path = 'input/band_parameters.csv'

  list_qqty = ["AD","DE","EDE"]
  list_stat = ["","_b","_t"]

  pd.options.mode.copy_on_write = True


  for scenario in dose.list_scenario:

    #if "B2G" not in scenario: continue
    dict_spe_doses = {"source":[],"date":[],"organ":[]}

    for qqty in list_qqty:
      for stat in list_stat:
        dict_spe_doses[qqty+stat] = []

    for source in ba.list_sources:
      list_storms = ba.getDatesSPE(source)
      for date_storm in list_storms:
        df_storm = ba.getSPEspectrum(source,date_storm,dose.listEp1)
        date = dt.datetime.strptime(date_storm, "%Y%m%d")

        for organ in dose.list_organs + ["all"]:
          #if organ!="all": continue
          print (source,date_storm,date,organ)
          df_f2d = dose.getFlux2DoseCoeffs(scenario,thick,organ)
          df_f2d = df_f2d[df_f2d["particle"]=="H"]
 
          list_sample = df_f2d["i_sample"].unique().tolist()
          list_AD = [0 for i in list_sample]
          list_DE = [0 for i in list_sample]
          list_EDE = [0 for i in list_sample]
          #if organ!="all": list_EDE = [np.nan for i in list_sample]

          for i_sample in list_sample:
            df_f2d_sample = df_f2d[df_f2d["i_sample"]==i_sample]
            #print (df_f2d_sample)
            #print (df_storm[date_storm])
            list_AD[i_sample] += float((df_f2d_sample["AD"] * df_storm[date_storm]).sum())
            list_DE[i_sample] += float((df_f2d_sample["DE"] * df_storm[date_storm]).sum())
            if organ=="all": list_EDE[i_sample] += float((df_f2d_sample["EDE"] * df_storm[date_storm]).sum())
              #df_f2d_sample["EDE"] = df_f2d_sample["EDE"] * df_storm[date_storm]
            #print (df_f2d_sample["DE"].sum())

          dict_spe_doses["source"].append(source)
          dict_spe_doses["date"].append(date)
          dict_spe_doses["organ"].append(organ)
          ad,ad_b,ad_t = calc_stats(list_AD)
          dict_spe_doses["AD"].append(ad)
          dict_spe_doses["AD_b"].append(ad_b)
          dict_spe_doses["AD_t"].append(ad_t)

          de,de_b,de_t = calc_stats(list_DE)
          dict_spe_doses["DE"].append(de)
          dict_spe_doses["DE_b"].append(de_b)
          dict_spe_doses["DE_t"].append(de_t)

          ede,ede_b,ede_t = calc_stats(list_EDE)
          dict_spe_doses["EDE"].append(ede)
          dict_spe_doses["EDE_b"].append(ede_b)
          dict_spe_doses["EDE_t"].append(ede_t)

    df_spe_doses = pd.DataFrame(dict_spe_doses)
    df_spe_doses.set_index("date",inplace=True)
    df_spe_doses.sort_index(inplace=True)

    for qqty in list_qqty:
      for stat in list_stat:
        col = qqty+stat
        df_spe_doses[col] = round(df_spe_doses[col],3) #/(24*60)*1000

    df_spe_doses.to_csv("output/SPE/"+scenario+"_"+str(thick)+".csv")
      
        
          #print (df_storm)
          #print (df_f2d)
     #sys.exit()
    
def generateSPE_time_doses():

  thick = 185
  d_increment  = 183

  dose = di.DoseInfo(True)
  sepem = sp.Sepem()

  for scenario in dose.list_scenario:
    #if scenario != "B2G-naked": continue
    for organ in ["all"]+dose.list_organs:
      datemin = dt.datetime(1974,7,1,0,0,0,0)
      #while datemin<dt.datetime(1974,7,13,0,0,0,0):
      while datemin<dt.datetime(2018,1,1,0,0,0,0):
        datemax = datemin + dt.timedelta(days=d_increment)
        filename = scenario+"_"+organ.replace(" ","-")+"_"+datemin.strftime("%Y%m%d")+"-"+datemax.strftime("%Y%m%d")+".csv"
        filepath = "output/SPE/tserie/"+filename
        if os.path.exists(filepath): 
          datemin = datemax
          continue
        with open(filepath, 'a'): os.utime(filepath, None)
        df_f2d = dose.getFlux2DoseCoeffs(scenario,thick,organ)
        df_f2d = df_f2d[df_f2d["particle"]=="H"]

        list_sample = df_f2d["i_sample"].unique().tolist()
        df_spe = sepem.getSpectrum(datemin,datemax,dose.listEp1)

        dict_time_vals = {"date":[],"AD":[],"AD_b":[],"AD_t":[],"DE":[],"DE_b":[],"DE_t":[]}
        if organ=="all":
          dict_time_vals["EDE"]=[]
          dict_time_vals["EDE_t"]=[]
          dict_time_vals["EDE_b"]=[]

        for storm in df_spe.columns:
          list_AD = [0 for i in list_sample]
          list_DE = [0 for i in list_sample]
          list_EDE = [0 for i in list_sample]
          for i_sample in list_sample:
            df_dose_sample = df_f2d[df_f2d["i_sample"]==i_sample]

            # [mSv cm2 / p ] * [p . cm-2 . min-1] = mSv / min
            list_AD[i_sample] += (df_dose_sample["AD"]*df_spe[storm]).sum() # each data point is for 5 minutes
            list_DE[i_sample] += (df_dose_sample["DE"]*df_spe[storm]).sum() # each data point is for 5 minutes
            if organ=="all": list_EDE[i_sample] += (df_dose_sample["EDE"]*df_spe[storm]).sum() # each data point is for 5 minutes

          date = dt.datetime.strptime(storm,"%Y%m%d-%H%M")
          dict_time_vals["date"].append(date)

          ad,ad_b,ad_t = calc_stats(list_AD)
          dict_time_vals["AD"].append(ad)
          dict_time_vals["AD_b"].append(ad_b)
          dict_time_vals["AD_t"].append(ad_t)

          de,de_b,de_t = calc_stats(list_DE)
          dict_time_vals["DE"].append(de)
          dict_time_vals["DE_b"].append(de_b)
          dict_time_vals["DE_t"].append(de_t)

          if organ=="all":
            ede,ede_b,ede_t = calc_stats(list_EDE)
            dict_time_vals["EDE"].append(ede)
            dict_time_vals["EDE_b"].append(ede_b)
            dict_time_vals["EDE_t"].append(ede_t)
    
        df_time_vals = pd.DataFrame(dict_time_vals)

        if len(df_time_vals)!=0:
          if not os.path.exists("output"): os.makedirs("output")
          if not os.path.exists("output/SPE"): os.makedirs("output/SPE")
          if not os.path.exists("output/SPE/tserie"): os.makedirs("output/SPE/tserie")
          df_time_vals.to_csv("output/SPE/tserie/"+filename,index=False)
        datemin = datemax
  #return
  print ("All file created")

def getDoses(datemin,datemax,qqty="DE"):

  dest_folder = "output/doses_storms/"
  list_organs = ['ureter', 'lymphatic nodes', 'tonsil', 'tongue', 'thyroid', 'thymus', 'teeth', 'stomach', 'muscle', 'trachea', 'spin', 'oral mucosa', 'pancreas', 'pituitary gland', 'oesophagus', 'skin', 'skin-insensitive', 'salivary gland', 'small intestine', 'extrathoracic', 'spleen', 'breast', 'lung', 'kidney', 'heart', 'gonad', 'gallbladder', 'eye-insensitive', 'extra tissue', 'content', 'liver', 'eye', 'adrenal', 'brain', 'bone surface', 'blood', 'bladder-insensitive', 'bladder', 'RBM', 'colon', 'uterus']
  list_stat = ["","_b","_t"]
  filename = qqty+"_"+datemin.strftime("%Y%m%d")+"_"+datemax.strftime("%Y%m%d")+".csv"
  filepath = dest_folder+filename
 # print (filepath)
  if os.path.exists(filepath):
    #print ("next")
    return
  with open(filepath, 'a'): os.utime(filepath, None)
  ndays_add = 10

  df = pd.DataFrame()
  print (datemin)
  for organ in list_organs:
    str_organ = organ.replace(" ","-")
    print ("\t",str_organ)

    file_path = 'output/SPE/final_tserie/B2G-naked_'+str_organ+'.csv'
    df_spe_naked = pd.read_csv(file_path)
    df_spe_naked['date'] = pd.to_datetime(df_spe_naked['date'])
    df_spe_naked = df_spe_naked[(df_spe_naked['date']>=datemin) & (df_spe_naked['date']<=datemax)] 
    df_spe_naked.set_index('date',inplace=True)
    df_spe_naked = df_spe_naked.drop_duplicates()
    for stat in list_stat: df["SPE_naked_"+str_organ+stat] = df_spe_naked[qqty+stat]
  
    file_path = 'output/SPE/final_tserie/B2G-vest_'+str_organ+'.csv'
    df_spe_vest = pd.read_csv(file_path)
    df_spe_vest['date'] = pd.to_datetime(df_spe_vest['date'])
    df_spe_vest = df_spe_vest[(df_spe_vest['date']>=datemin) & (df_spe_vest['date']<=datemax)]
    df_spe_vest.set_index('date',inplace=True)
    df_spe_vest = df_spe_vest.drop_duplicates()
    for stat in list_stat: df["SPE_vest_"+str_organ+stat] = df_spe_vest[qqty+stat]

    file_path = 'output/GCR/final_tserie/B2G-naked.csv'
    df_gcr_naked = pd.read_csv(file_path)
    df_gcr_naked['date'] = pd.to_datetime(df_gcr_naked['date'])
    df_gcr_naked = df_gcr_naked[df_gcr_naked["organ"]==organ] 
    df_gcr_naked = df_gcr_naked[(df_gcr_naked['date']>=datemin-dt.timedelta(days=ndays_add))\
                                & (df_gcr_naked['date']<=datemax+dt.timedelta(days=ndays_add))]
    df_gcr_naked.set_index('date',inplace=True)
    df_gcr_naked.drop("organ",inplace=True,axis=1)
    df_gcr_naked = df_gcr_naked.drop_duplicates()
    df_gcr_naked = df_gcr_naked.resample('5Min').interpolate(method='linear')
    df_gcr_naked = df_gcr_naked[(df_gcr_naked.index>=datemin) & (df_gcr_naked.index<=datemax)]
    for stat in list_stat: df["GCR_naked_"+str_organ+stat] = df_gcr_naked[qqty+stat]

    file_path = 'output/GCR/final_tserie/B2G-vest.csv'
    df_gcr_vest = pd.read_csv(file_path)
    df_gcr_vest['date'] = pd.to_datetime(df_gcr_vest['date'])
    df_gcr_vest = df_gcr_vest[df_gcr_vest["organ"]==organ] 
    df_gcr_vest = df_gcr_vest[(df_gcr_vest['date']>=datemin-dt.timedelta(days=ndays_add))\
                                & (df_gcr_vest['date']<=datemax+dt.timedelta(days=ndays_add))]
    df_gcr_vest.set_index('date',inplace=True)
    df_gcr_vest.drop("organ",inplace=True,axis=1)
    df_gcr_vest = df_gcr_vest.drop_duplicates()
    df_gcr_vest = df_gcr_vest.resample('5Min').interpolate(method='linear')
    df_gcr_vest = df_gcr_vest[(df_gcr_vest.index>=datemin) & (df_gcr_vest.index<=datemax)]
    for stat in list_stat: df["GCR_vest_"+str_organ+stat] = df_gcr_vest[qqty+stat]
  
  df.to_csv(filepath)
     
     
#datemin = dt.datetime(1989,7,1,0,0,0,0)
#datemax = dt.datetime(1989,8,1,0,0,0,0)
#generateSPE_time_doses()
#generateGCRDoses()

#sys.exit()

thick = 185
#organ = sys.argv[1]
#datemin = dt.datetime.strptime(sys.argv[2],"%Y%m%d")
#datemax = dt.datetime.strptime(sys.argv[3],"%Y%m%d")


t1 = time.time()

df_rel = pd.read_csv("input/SEPEM_RDS_v3/SEPEM_REL_v2/RELv2.csv")
df_rel["start"] = pd.to_datetime(df_rel["start date"])
df_rel["end"] = pd.to_datetime(df_rel["end date"])
df_rel = df_rel[["start","end"]]

for index, row in df_rel.iterrows():
  year = row["start"].year
  datemin = row["start"]-dt.timedelta(hours=12)
  datemax = row["end"]
  getDoses(datemin,datemax,"DE")

#df_spe = generateSPE_time_doses()
#generateGCRDoses(datemin,datemax)

print ("time to run:",time.time()-t1)

#scenario = "B2G-vest"
#df_spe_B2G_vest = generateSPE_time_doses(datemin,datemax,scenario,thick,organ,qqty)
#df_spe_B2G_vest.to_csv("output/SPE/B2G_vest_tserie.csv",index=False)


   
#generateSPEDoses()

#datemin = dt.datetime.strptime(sys.argv[1],"%Y%m%d")
#datemax = dt.datetime.strptime(sys.argv[2],"%Y%m%d")

