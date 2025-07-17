import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import numpy as np
import datetime
import math
from scipy.optimize import curve_fit
import sys, os
import re

import netCDF4
import wget

from physics import *
from functions import *


## Download goes nc file 
def getFile(year,month,day,goes="goes16"):

  nofile_list = []
  with open("data/Goes/nodata.txt") as file:
    nofile_list = [line.rstrip() for line in file]

  dateFile = datetime.date(year,month,day)

  folder = "data/Goes/"+goes+"/"
  if not os.path.exists(folder):
    os.makedirs(folder)

  shortgoes = goes.replace("oes","")
  strdate = str(year)
  if month<10: strdate+="0"
  strdate+=str(month)
  if day<10: strdate+="0"
  strdate+=str(day)

  dest_filename = shortgoes+"_"+strdate+".nc"
  if os.path.isfile(folder+dest_filename) or folder+dest_filename in nofile_list: 
    #print (dest_filename)
    return dest_filename


  # Case storm September 2017
  if year < 2018:
    if (goes=="goes16"):
      datemin = datetime.date(2017,9,1)
      datemax = datetime.date(2017,10,1)


      if dateFile>=datemin and dateFile<datemax:
        dest_fileName = shortgoes+"_"+strdate+".nc"
        source_fileName = "se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc"
        url = "https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/goesr/solar_proton_events/sgps_sep2017_event_data/"
        url+=source_fileName
        try: wget.download(url,out=folder+dest_filename)
        except: 
          fnodata = open("data/Goes/nodata.txt","a")
          fnodata.write(folder+dest_filename+"\n")  
          fnodata.close()
          return ""
        return dest_fileName
    else:
      print ("no data")
      return ""

  version = "v1-0-1"
  if dateFile>datetime.date(2023,10,13): version = "v3-0-2"
  elif dateFile>datetime.date(2023,3,30): version = "v3-0-1"
  elif dateFile>datetime.date(2022,4,25): version = "v3-0-0"
  elif dateFile>datetime.date(2021,8,11): version = "v2-0-0"

  source_filename = "sci_sgps-l2-avg5m_"+shortgoes+"_d"+strdate+"_"+version+".nc"
  url = "httpsdata.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/"+goes+"/l2/data/sgps-l2-avg5m/"
  url+=str(year)+"/"
  if month<10: url+="0"
  url+=str(month) + "/" + source_filename
  try: wget.download(url,out=folder+dest_filename)
  except: 
    fnodata = open("data/Goes/nodata.txt","a")
    fnodata.write(folder+dest_filename+"\n")  
    fnodata.close()
    return ""
  
  #!wget -q -O fodler+dest_filename $url
  return dest_filename

def getGoesData(datestart,dateend,goes,tocsv="",nunit=1):

  #folder = "data/goes/"+goes+"/"
  #if not os.path.exists(folder):
  #  os.makedirs(folder)
  date = datetime.datetime(datestart.year,datestart.month,datestart.day,0,0,0,0)

  folder = "data/Goes/"+goes+"/"




  if (datestart>=datetime.datetime(2017,9,1)) and (dateend<datetime.datetime(2017,10,1)):
    if goes == "goes16":
      return get2017StormFlux(nunit)
    else: return pd.DataFrame()

  df_fluxes = pd.DataFrame()
  d_fluxes = {}

  while date<dateend:

    year = date.year
    month = date.month
    day = date.day
    filename = getFile(year,month,day,goes)

    #filename,url = getFileName(year,month,day,"goes16")
    #print (goes,date)

    try: file2read = netCDF4.Dataset(folder+"/"+filename,'r')
    except:
      #print ("File not found")
      date = date + datetime.timedelta(days=1)
      continue

    #file2read = netCDF4.Dataset('sci_sgps-l2-avg5m_g18_d20240510_v3-0-2.nc','r')
    #print (file2read.variables['time'])
    #print (file2read.variables)
    #print (file2read.variables['L2_SciData_TimeStamp'])

    timeKey = "time"
    if "L2_SciData_TimeStamp" in file2read.variables.keys(): timeKey = "L2_SciData_TimeStamp"
    xtime = file2read.variables[timeKey][:]
    xtime = [datetime.datetime(2000,1,1,12,0,0,0)+datetime.timedelta(seconds=x) for x in xtime]

    #break
    #xtime = file2read.variables['time'][:]
    # 5 minutes interval bins
    #xtime = [datetime.datetime(2000,1,1,12,0,0,0)+datetime.timedelta(seconds=x) for x in xtime]
    d_fluxes["time"] = xtime

    lowerE = file2read.variables['DiffProtonLowerEnergy'][:]
    upperE = file2read.variables['DiffProtonUpperEnergy'][:]

    #print (file2read.variables['DiffProtonLowerEnergy'][:])
    #print (upperE)
    #print (file2read.variables["AvgIntProtonFlux"])
    AvgDiffProtonFlux = file2read.variables['AvgDiffProtonFlux'][:] #protons/(cm^2 sr keV s)
    AvgIntProtonFlux = file2read.variables["AvgIntProtonFlux"][:] #protons/(cm^2 sr s)

    #print (AvgDiffProtonFlux[:,1,0])
    it = 0
    for nsensor in range(13):
      d_fluxes[str(lowerE[nunit,nsensor])+"-"+str(upperE[nunit,nsensor])] = AvgDiffProtonFlux[:,nunit,nsensor]

    d_fluxes[">500"] = AvgIntProtonFlux[:,nunit]

    #print (df_fluxes)
    df_fluxes = pd.concat([df_fluxes,pd.DataFrame(d_fluxes)])


    date = date + datetime.timedelta(days=1)
    #break
  #print (file2read.variables.items())
  if len(df_fluxes)==0: return df_fluxes
  df_fluxes.set_index("time",inplace=True)

  df_fluxes = df_fluxes[(df_fluxes.index>=datestart) & (df_fluxes.index<dateend)]
  if tocsv!="": df_fluxes.to_csv(tocsv)
  return df_fluxes

def get2017StormFlux(nunit=0):

  df_fluxes = pd.DataFrame()
  d_fluxes = {}

  year,month,day = 2017,9,15

  getFile(2017,9,15,"goes16")

  goes = "goes16"
  shortgoes = goes.replace("oes","")
  strdate = str(year)
  if month<10: strdate+="0"
  strdate+=str(month)
  if day<10: strdate+="0"
  strdate+=str(day)

  folder = "data/Goes/goes16/"
  dest_fileName = shortgoes+"_"+strdate+".nc"

  file2read = netCDF4.Dataset(folder+'/'+dest_fileName,'r')
  xtime = file2read.variables["L2_SciData_TimeStamp"][:]
  xtime = [datetime.datetime(2000,1,1,12,0,0,0)+datetime.timedelta(seconds=x) for x in xtime]

  d_fluxes["time"] = xtime

  lowerE = file2read.variables['DiffProtonLowerEnergy'][:]
  upperE = file2read.variables['DiffProtonUpperEnergy'][:]

  AvgDiffProtonFlux = file2read.variables['AvgDiffProtonFlux'][:] #protons/(cm^2 sr keV s)
  AvgIntProtonFlux = file2read.variables["AvgIntProtonFlux"][:] #protons/(cm^2 sr s)

  it = 0
  for nsensor in range(13):
    d_fluxes[str(lowerE[nunit,nsensor])+"-"+str(upperE[nunit,nsensor])] = AvgDiffProtonFlux[:,nunit,nsensor]

  d_fluxes[">500"] = AvgIntProtonFlux[:,nunit]

  #print (df_fluxes)
  df_fluxes = pd.DataFrame(d_fluxes)
  df_fluxes.set_index("time",inplace=True)
  return df_fluxes

def convert_gglsht_url(url):
    pattern = r'https://docs\.google\.com/spreadsheets/d/([a-zA-Z0-9-_]+)(/edit#gid=(\d+)|/edit.*)?'
    replacement = lambda m: f'https://docs.google.com/spreadsheets/d/{m.group(1)}/export?' + (f'gid={m.group(3)}&' if m.group(3) else '') + 'format=csv'
    new_url = re.sub(pattern, replacement, url)
    return new_url

def getStormIntFlux(df_fluxes,df_background,datemin,datemax):
  #df_fluxes_short = df_fluxes[(df_fluxes.index>datemin) & (df_fluxes.index<datemax)]
  df_back_av = df_background.mean(axis=0)
  df_storm = df_fluxes.copy()

  for c in df_storm.columns:
    df_storm[c] -= df_back_av[c]
    df_storm[c] = df_storm[c].mask(df_storm[c] < 0).fillna(float('NaN'))

  flux500 = float(df_storm[">500"].sum()*4*np.pi*300)
  df_storm.drop(">500",axis=1,inplace=True)

  namecol = datemin.strftime("%Y-%m-%d")

  df_sum = pd.DataFrame()
  emin = []
  emax = []
  for erange in df_storm.columns:
    if erange==">500": continue
    emin.append(float(erange.split('-')[0]))
    emax.append(float(erange.split('-')[1]))
  df_sum["emin"] = emin
  df_sum["emax"] = emin[1:]+[500000]

  df_sum[namecol] = pd.DataFrame(df_storm.sum(axis=0), columns=[namecol])[namecol].tolist()
  df_sum[namecol] = df_sum[namecol]*(df_sum["emax"]-df_sum["emin"])*4*np.pi*300

  diff = df_sum[namecol].tolist()

  dict_int = {"emin":[]}
  for e in emin:
    dict_int["emin"].append(e/1000)
  dict_int["emin"].append(500)

  dict_int["storm"] = [0 for i in range(len(dict_int["emin"]))]
  dict_int["storm"][-1] = flux500

  #print (diff)

  for i in range(12,-1,-1):
    #print (i,dict_int["storm"])
    dict_int["storm"][i] = dict_int["storm"][i+1] + diff[i]
  #print (i,dict_int["storm"])
  df_int = pd.DataFrame(dict_int)
  return df_int


# Generate a list of 13 colors from blue to red using a colormap
colors = cm.get_cmap('viridis', 13)(np.linspace(0, 1, 13))

# Convert the colors to hexadecimal format if needed
hex_colors = [matplotlib.colors.to_hex(color) for color in colors]

url = 'https://docs.google.com/spreadsheets/d/1TIP0TXqtUwRWR6jGmsS7KNpPAvjlVwkQbvXKq0yvN88/edit?usp=drive_link'
new_url = convert_gglsht_url(url)
print(new_url)

df_storm_dates = pd.read_csv(new_url)

for c in ["init_storm","end_storm","init_background","end_background"]:
  df_storm_dates[c] = pd.to_datetime(df_storm_dates[c],format='%d/%m/%Y %H:%M')
print (df_storm_dates)


list_goes = ["goes16","goes17","goes18"]
color = ["r",'g','b']
last_popt = [1.0e4*300*4*np.pi, 1.0e-01, 5.36086367e+00, 8.83179720e-03]

loglistE = np.linspace(0,3,100)
listE = [10**(e) for e in loglistE]
listR = [rigidity(e,mass_proton,1) for e in listE]

#stormName = "Sep-17"

list_df_storms = {}
#row = df_storm_dates[df_storm_dates["name"]==stormName].iloc[0]

list_band_parameters = {"Date":[],"end date":[],"goes":[],"J0":[],"gamma1":[],"gamma2":[],"R0":[]}

for i,row in df_storm_dates.iterrows():
#if True:
  stormName = row["name"]

  print (stormName)
  ngoes=2
  if stormName in ["Mar-22","Feb-23"]: ngoes = 3

  #stormRow = df_storm_dates[df_storm_dates["name"]==stormName].iloc[0]
  datemin = min(row["init_storm"],row["init_background"])-datetime.timedelta(days=1)
  datemax = max(row["end_storm"],row["end_background"])+datetime.timedelta(days=1)

  date_storm_init = row["init_storm"]
  date_storm_end = row["end_storm"]

  date_background_init = row["init_background"]
  date_background_end = row["end_background"]

  nunit = 0

  list_df_storm = {}

  for ig, goes in enumerate(list_goes):

    print ("-",goes)

    df_storm = getGoesData(datemin,datemax,goes,nunit=nunit)

    #continue
    if len(df_storm)==0: continue
    #print (df_storm[(df_storm.index<date_background_end)])
    df_background = df_storm[(df_storm.index>=date_background_init) & (df_storm.index<date_background_end)]
    #print (len(df_background))
    df_storm_cote = df_storm[(df_storm.index>=date_storm_init) & (df_storm.index<date_storm_end)]

    if ((len(df_storm_cote)==0) or (len(df_background)==0)): continue
    df_int = getStormIntFlux(df_storm_cote,df_background,date_storm_init,date_storm_end)
    df_int = df_int[df_int["storm"]>0]
    if (len(df_int)==0): continue

    R_fit = [rigidity(e,mass_proton,1) for e in df_int["emin"].tolist()]
    total_fit = df_int["storm"].tolist()
    #print (total_fit)
    log_total_fit = np.log10(total_fit)

    band_popt, pcov = curve_fit(logBand, R_fit, log_total_fit, p0=last_popt, bounds=([0.0,-10.0,0.0,0.0], [np.inf,10.0,100.0,10.0]), maxfev=100000000)
    list_band_parameters["Date"].append(date_storm_init)
    list_band_parameters["end date"].append(date_storm_end)
    list_band_parameters["goes"].append(goes)
    list_band_parameters["J0"].append(band_popt[0])
    list_band_parameters["gamma1"].append(band_popt[1])
    list_band_parameters["gamma2"].append(band_popt[2])
    list_band_parameters["R0"].append(band_popt[3])
    #print (goes,band_popt)
    fitband = lBand(listR,*band_popt)

    poly_popt, pcov = curve_fit(logpolynomial, R_fit, log_total_fit, maxfev=100000000)
    #print (goes,poly_popt)
    fitpoly = lpolynomial(listR,*poly_popt)

    list_df_storm[goes] = {}
    list_df_storm[goes]["df_storm"] = df_storm
    list_df_storm[goes]["df_int"] = df_int
    list_df_storm[goes]["fitband"] = fitband
    list_df_storm[goes]["fitpoly"] = fitpoly
  list_df_storms[stormName] = list_df_storm

df_list_goes_band = pd.DataFrame(list_band_parameters)
df_list_goes_band.to_csv("data/Goes/band_goes.csv",index=False)

def plotStorms():
  row = df_storm_dates[df_storm_dates["name"]==stormName].iloc[0]

  for i,row in df_storm_dates.iterrows():

  #if True:

    stormName = row["name"]

    datemin = min(row["init_storm"],row["init_background"])-datetime.timedelta(days=1)
    datemax = max(row["end_storm"],row["end_background"])+datetime.timedelta(days=1)

    date_storm_init = row["init_storm"]
    date_storm_end = row["end_storm"]

    date_background_init = row["init_background"]
    date_background_end = row["end_background"]


    #if stormName not in list_df_storms.keys(): continue
 
    list_df_storm = list_df_storms[stormName]
    ngoes = len(list_df_storm.keys())

    #fig,ax = plt.subplots(ngoes*2,figsize=(30,10*ngoes),sharex=True)
    fig = plt.figure(figsize=(45,20))
    gs = fig.add_gridspec(ngoes*2, 2, width_ratios=[1, 1])
    plt.suptitle(date_storm_init.strftime("%Y-%m-%d")+" - "+date_storm_end.strftime("%Y-%m-%d"),fontsize=30,y=0.91)
    ax3 = fig.add_subplot(gs[:, 1])
    iplot = 0
  
  
    mind = datetime.datetime(2300,1,1,0,0,0,0)
    maxd = datetime.datetime(1900,1,1,0,0,0,0)
  
    for goes,dict_goes in list_df_storm.items():
  
      df_storm = dict_goes["df_storm"]
      df_int = dict_goes["df_int"]
      fitband = dict_goes["fitband"]
      fitpoly = dict_goes["fitpoly"]
  
      if df_storm.index[0]<mind: mind = df_storm.index[0]
      if df_storm.index[-1]>maxd: maxd = df_storm.index[-1]
  
      ax1 = fig.add_subplot(gs[2*iplot, 0])
      ax2 = fig.add_subplot(gs[2*iplot+1, 0])
      #if iplot%2==0:
  
      #else:
        #ax3 = fig.add_subplot(gs[2:, 1])
  
      #if iplot==0:
      #  ax1.set_title(date_storm_init.strftime("%Y-%m-%d")+" - "+date_storm_end.strftime("%Y-%m-%d"))
  
      for i,c in enumerate(df_storm.columns):
        if ">500"==c: continue
        lab = str(float(c.split('-')[0])/1000)+"-"+str(float(c.split('-')[1])/1000)
        ax1.plot(df_storm.index,df_storm[c],label=lab,color=hex_colors[i])
      ax2.plot(df_storm.index,df_storm[">500"])
      ax1.set_yscale('log')
      ax1.text(0.5, 0.95, goes, transform=ax1.transAxes, ha="center", va="top",fontsize=17)
      ax1.grid(True)
      ax2.grid(True)
  
      ax1.set_ylabel("diff flux #/(cm2.s.sr.keV)",fontsize=18)
      ax2.set_ylabel("integrated flux #/(cm2.sr.s)",fontsize=18)
  
      ax1.axvline(x=date_storm_init,color='darkred',lw=1)
      ax1.axvline(x=date_storm_end,color='darkred',lw=1)
      ax2.axvline(x=date_storm_init,color='darkred',lw=1)
      ax2.axvline(x=date_storm_end,color='darkred',lw=1)
  
      ax1.axvspan(date_storm_init, date_storm_end, color='darkred', alpha=0.1)
      ax2.axvspan(date_storm_init, date_storm_end, color='darkred', alpha=0.1)
  
      ax1.axvline(x=date_background_init,color='darkgreen',lw=1)
      ax1.axvline(x=date_background_end,color='darkgreen',lw=1)
      ax2.axvline(x=date_background_init,color='darkgreen',lw=1)
      ax2.axvline(x=date_background_end,color='darkgreen',lw=1)
  
      ax1.axvspan(date_background_init, date_background_end, color='darkgreen', alpha=0.1)
      ax2.axvspan(date_background_init, date_background_end, color='darkgreen', alpha=0.1)
  
      ax1.tick_params(axis='both', which='major', labelsize=17)
      ax1.tick_params(axis='both', which='major', labelsize=17)
  
    #ax[0].set_title(datemin.strftime("%Y-%m-%d")+" - "+datemax.strftime("%Y-%m-%d"))
      if iplot==0: ax1.legend(ncol=3,fontsize=12)
      iplot = iplot+1
  
  
    #for goes in list_goes:
      #df_storm_cote = df_storm[(df_storm.index>=date_storm_init) & (df_storm.index<date_storm_end)]#getGoesData(date_storm_init,date_storm_end,goes,"sgps_"+goes+".csv",nunit=nunit)
      #df_background = df_storm[(df_storm.index>=date_background_init) & (df_storm.index<date_background_end)]
      #df_storm = df_period[(df_period.index>date_storm_init) & (df_period.index<datemax)]
      #if len(df_storm)==0: continue
      #df_int = getStormIntFlux(df_storm_cote,df_background,date_storm_init,date_storm_end)
  
      #logBand = np.log10(df_int["emin"].tolist())
      #R_fit = [rigidity(e,mass_proton,1) for e in df_int["emin"].tolist()]
      #total_fit = df_int["storm"].tolist()
      #log_total_fit = np.log10(total_fit)
  
      ax3.plot(df_int["emin"],df_int["storm"],'X',markersize=14,label=goes+", unit="+str(nunit),color=color[list_goes.index(goes)])
  
      #band_popt, pcov = curve_fit(logBand, R_fit, log_total_fit, p0=last_popt, bounds=([0.0,-10.0,0.0,0.0], [np.inf,10.0,100.0,10.0]), maxfev=100000000)
      #print (goes,band_popt)
      #fitband = lBand(listR,*band_popt)
      ax3.plot(listE,fitband,color=color[list_goes.index(goes)],label="band fit")
  
  
      #poly_popt, pcov = curve_fit(logpolynomial, R_fit, log_total_fit, maxfev=100000000)
      #print (goes,poly_popt)
      #fitpoly = lpolynomial(listR,*poly_popt)
      ax3.plot(listE,fitpoly,color=color[list_goes.index(goes)],linestyle=":",label="logPoly3 fit")
  
    #
    #print (mind,maxd)
  
    for ig in range(1,2*len(list_df_storm)+1):
      ax1 = fig.axes[ig]
      ax1.set_xlim(mind,maxd)
  
  
    ax3.tick_params(axis='both', which='major', labelsize=20)
    ax3.set_ylabel("Event integrated flux #/cm2",fontsize=20)
    ax3.set_xlabel("Energy (MeV)",fontsize=20)
    ax3.legend(fontsize=20)
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.grid(True)
  
    plt.subplots_adjust(wspace=0.11)
    plt.subplots_adjust(hspace=0)
    plt.savefig("Storms/"+stormName+".png",bbox_inches='tight')
    plt.close()

list_df_integrated_fluxes = {}
list_goes = ["goes16","goes17","goes18"]
nunit = 0

for goes in list_goes:
  list_df_integrated_fluxes[goes] = pd.DataFrame()

for i,row in df_storm_dates.iterrows():
  nameStorm = row["name"]
  #print (rowStorm)

  date_storm_init = row["init_storm"]
  date_storm_end = row["end_storm"]

  date_background_init = row["init_background"]
  date_background_end = row["end_background"]
  for goes in list_goes:

    df_storm = getGoesData(date_storm_init,date_storm_end,goes,nunit=nunit)
    df_background = getGoesData(date_background_init,date_background_end,goes,nunit=nunit)
    if len(df_storm)==0: continue
    df_int = getStormIntFlux(df_storm,df_background,date_storm_init,date_storm_end)
    if "emin" not in list_df_integrated_fluxes[goes].columns: list_df_integrated_fluxes[goes]["emin"] =  df_int["emin"]

    list_df_integrated_fluxes[goes][nameStorm] = df_int["storm"].tolist()


for goes in list_df_integrated_fluxes.keys():
  list_df_integrated_fluxes[goes].to_csv("data/Goes/"+goes+"/integrated_fluxes.csv",index=False)
