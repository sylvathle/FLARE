import pandas as pd
import math
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from functions import *
from physics import *

# Generate a list of 13 colors from blue to red using a colormap
colors = cm.get_cmap('viridis', 13)(np.linspace(0, 1, 13))

# Convert the colors to hexadecimal format if needed
hex_colors = [matplotlib.colors.to_hex(color) for color in colors]



df_band = pd.read_csv("data/band_tylka-hu.csv")
df_band["Date"] = pd.to_datetime(df_band["Date"])
enddate = []
for i,row in df_band.iterrows():
  if math.isnan(row["Duration"]):
    enddate.append(None)
  else:
    enddate.append(row["Date"]+datetime.timedelta(hours=row["Duration"]))
df_band["end date"] = enddate
df_band = df_band[["Date","end date","Source","GLE","J0","gamma1","gamma2","R0"]]

df_hu = df_band[df_band["Source"]=="Hu"]

listGLE = df_hu["GLE"].unique()

sepem_folder = "data/SEPEM_RDS_v3/"

df_rel = pd.read_csv(sepem_folder+"SEPEM_REL_v2/"+"RELv2.csv")
df_rel["start date"] = pd.to_datetime(df_rel["start date"])
df_rel.rename(columns={"start date":"Date"},inplace=True)
df_rel["end date"] = pd.to_datetime(df_rel["end date"])

df = pd.read_csv(sepem_folder+"SEPEM_RDS_v3-2/"+"RDS3.2_FPDO_FPIO.csv")
df["epoch"] = pd.to_datetime(df["epoch"])
df.set_index("epoch",inplace=True)

loglistE = np.linspace(0,3,100)
listE = [10**(e) for e in loglistE]
listR = [rigidity(e,mass_proton,1) for e in listE]

#dateminS = datetime.datetime(2003,10,26,0,0,0,0)
#datemaxS = datetime.datetime(2003,11,1,0,0,0,0)
fpio_cols = ["fpio_"+str(i) for i in range(1,8)]
fpio_min = [5,10,30,50,100,300,500]


#Flux are given in cm-2.sr-1.s-1 each 5 minutes
# We convert it in cm-2.5min-1
df_fpio = df[fpio_cols]*4.0*np.pi*300
print (df_fpio)

sys.exit()


fpdo_cols = ["fpdo_"+str(i) for i in range(1,15)]
fpdo_ebin = [[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],[21.87,31.62]]
fpdo_ebin += [[31.62,45.73],[45.73,66.13],[66.13,95.64],[95.64,138.3],[138.3,200.0]]
fpdo_ebin += [[200.0,289.2],[289.2,418.3],[418.3,604.9],[604.9,874.7]]

fpdo_meanE = [6.01,8.70,12.57,18.18,26.30,38.03,54.99,79.53,115.01,166.31,240.51,347.81,502.97,727.36]

df_fpdo = df[fpdo_cols]

df_storms = df_fpio[fpio_cols]
#df_storms = df_fpio[fpdo_cols]
#df_storms = df_storms[(df_storms.index>dateminS) & (df_storms.index<datemaxS)]


band_J0 = []
band_gamma1 = []
band_gamma2 = []
band_R0 = []
band_lin_mse = []
band_log_mse = []

rpow_a = []
rpow_b = []
rpow_lin_mse = []
rpow_log_mse = []

epow_a = []
epow_b = []
epow_lin_mse = []
epow_log_mse = []

#df_rel = df_rel.sort_values("Fluence",ascending=False)

df_rel_short = df_rel.copy() #.iloc[:10]
#df_rel_short = df_rel_short[(df_rel_short["start date"]>dateminS) & (df_rel_short["start date"]<datemaxS)]


print (df_rel_short)

plot = True

#list_GLE = df_tylka["Official GLE No"].unique()

  #fig,ax = plt.subplots(len(df_tylka),figsize=(10,3*len(df_tylka)),sharex=True)
  #lE = [10**e for e in np.linspace(0,3)]
  #lR = [rigidity(e,mass_proton,1) for e in lE]

nplot =0

datesmin = []
datesmax = []



#for gle in listGLE:
for irow, row_rel in df_rel.iterrows():
 # df_gle = df_band[df_band["GLE"]==gle]
  #row_rel = df_gle.iloc[0]
  #print (row_rel)
  #datemin = row_rel["start date"]
  datemin = row_rel["Date"]
  datemax = row_rel["end date"]

  datesmin.append(datemin)
  datesmax.append(datemax)

  strdate = row_rel["Date"].strftime("%Y-%m-%d")
  date = row_rel["Date"]
  strdate = date.strftime("%Y-%m-%d")


  #for i, row in df_tylka.iterrows():
  #  temp_hu = lBand(lR,row["J0"],row["gamma1"],row["gamma2"],row["R0"])
  #  for i in range(len(lR)): huband[i] += temp_hu[i]

#for id,row_rel in df_rel_short.iterrows():
#for id,row_rel in df_tylka.iterrows():
  #print (datemin,datemax)

  df_storm = df_storms [(df_storms.index>datemin) & (df_storms.index<datemax)]
  #print (df_storm)
  #break
  # Calculate the sum of each column
  totals = df_storm.sum().tolist()

  R = [rigidity(e,mass_proton,1) for e in fpio_min]

  #last_popt = [row_rel["Fluence"], 1.0e-01, 5.36086367e+00, 8.83179720e-03]
  last_popt = [1.0e4*300*4*np.pi, 1.0e-01, 5.36086367e+00, 8.83179720e-03]

  E_fit = []
  R_fit = []
  total_fit = []
  for i in range(len(totals)):
    if totals[i]>0:
      E_fit.append(fpio_min[i])
      R_fit.append(R[i])
      total_fit.append(totals[i])

#
  band_popt, pcov = curve_fit(logBand, R_fit, np.log10(total_fit), p0=last_popt, bounds=([0.0,-10.0,0.0,0.0], [np.inf,10.0,100.0,10.0]), maxfev=100000000)
  if plot:
    fig,ax = plt.subplots(1,2,figsize=(15,6),sharex=False)
    fitband = lBand(listR,*band_popt)
    #huband = lBand(lR,row_rel["J0"],row_rel["gamma1"],row_rel["gamma2"],row_rel["R0"])

  band_J0.append(float(band_popt[0]))
  band_gamma1.append(float(band_popt[1]))
  band_gamma2.append(float(band_popt[2]))
  band_R0.append(float(band_popt[3]))
  band_lin_mse.append(mse_band(R_fit,total_fit,lBand,band_popt))
  band_log_mse.append(mse_band(R_fit,np.log10(total_fit),logBand,band_popt))

  fit_hu = [0 for r in R_fit]
  temp_band = lBand(R_fit,band_popt[0],band_popt[1],band_popt[2],band_popt[3])
  for i in range(len(temp_band)):
    fit_hu[i] += temp_band[i]

  log_fit_hu = np.log10(fit_hu)

  err_hu = mse(fit_hu,total_fit)
  log_err_hu = mse(log_fit_hu,np.log10(total_fit))

  rpow_popt,pcov = curve_fit(logpowerLaw,R_fit, np.log10(total_fit))
  if plot: fit_rpow = powerLaw(listR,*rpow_popt)
  rpow_a.append(rpow_popt[0])
  rpow_b.append(rpow_popt[1])
  rpow_lin_mse.append(mse_band(R_fit,total_fit,powerLaw,rpow_popt))
  rpow_log_mse.append(mse_band(R_fit,np.log10(total_fit),logpowerLaw,rpow_popt))

  epow_popt,pcov = curve_fit(logpowerLaw,E_fit, np.log10(total_fit))
  if plot: fit_epow = powerLaw(listE,*epow_popt)
  epow_a.append(epow_popt[0])
  epow_b.append(epow_popt[1])
  epow_lin_mse.append(mse_band(E_fit,total_fit,powerLaw,epow_popt))
  epow_log_mse.append(mse_band(E_fit,np.log10(total_fit),logpowerLaw,epow_popt))

  if plot:
    for ic,col in enumerate(fpio_cols):
      lab = ">" + str(fpio_min[ic])+" MeV"
      ax[0].plot(df_storm.index,df_storm[col],label=lab,color=hex_colors[ic])
    ax[1].plot(fpio_min,totals,'x','b')
    ax[1].plot(listE,fitband,'r',label = "band linErr="+str(round(band_lin_mse[-1]*100,1))+"% logErr="+str(round(band_log_mse[-1]*100,1))+"%")
    #ax[nplot].plot(listE,huband,'g',label = "hufit band linErr="+str(round(err_hu*100,1))+"% logErr="+str(round(log_err_hu*100,1))+"%")
    ax[1].plot(listE,fit_rpow,'g',label = "rpow linErr="+str(round(rpow_lin_mse[-1]*100,1))+"% logErr="+str(round(rpow_log_mse[-1]*100,1))+"%")
    ax[1].plot(listE,fit_epow,'b',label = "epow linErr="+str(round(epow_lin_mse[-1]*100,1))+"% logErr="+str(round(epow_log_mse[-1]*100,1))+"%")
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[0].set_yscale('log')
    ax[0].tick_params(axis='x', labelrotation=30)
    ax[0].set_ylabel(r'Integrated proton flux $1/(cm^2.sr.s)$')
    ax[1].set_ylabel(r'Event integrated proton flux $1/(cm^2.sr)$')
    ax[0].legend(loc=1)
    ax[1].legend(loc=3)
    nplot += 1

    plt.suptitle(datemin.strftime("%Y-%m-%d %H:%M")+" - "+datemax.strftime("%Y-%m-%d %H:%M"),y=0.92)
#if plot:
    plt.subplots_adjust(wspace=0.15)
    plt.savefig("Storm_sepem/"+strdate+".png",bbox_inches="tight")
    plt.close()


dict_hu = {}
dict_hu["Date"] = datesmin
dict_hu["end date"] = datesmax
dict_hu["J0"] = band_J0
dict_hu["gamma1"] = band_gamma1
dict_hu["gamma2"] = band_gamma2
dict_hu["R0"] = band_R0

df_hu = pd.DataFrame(dict_hu)
df_hu.to_csv("band_sepem_fit.csv")
  #break

  #print ("linear mse = ",str(round(mse(R_fit,total_fit,lBand,popt)*100,1))+"%")
  #print ("logar mse = ",str(round(mse(R_fit,np.log10(total_fit),logBand,popt)*100,1))+"%")
# Create a new row with the totals
#total_row = pd.DataFrame(totals).transpose()

# Append the total row to the DataFrame
#df_storm_with_total = pd.concat([df_storm, total_row], ignore_index=True)

#Rename the last row
#df_storm_with_total.rename(index={len(df_storm_with_total)-1:"Total"},inplace=True)
#print (df_storm_with_total)

