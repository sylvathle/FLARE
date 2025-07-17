import pandas as pd
import math
import datetime

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


print (df_band)

#df_band["Date"] = pd.to_datetime(df_band["Date"])
df_band["end date"] = pd.to_datetime(df_band["end date"])

df_tylka = df_band[df_band["Source"]=="Tylka"]

print (df_tylka)
df_hu = df_band[df_band["Source"]=="Hu"]

listGLE = df_hu["GLE"].unique()

#listLogE = np.linspace(0,4,100)
#listE = [float(10**e) for e in listLogE]
#listR = [rigidity(e,mass_proton,1) for e in listE]

#dict_hu = {"E":listE,"R":listR}

print (listGLE)

for gle in listGLE:
  listrow = df_hu[df_hu["GLE"]==gle]
  row = listrow.iloc[0]
  strdate = row["Date"].strftime("%Y-%m-%d")
  dict_hu[strdate] = [0 for r in listR]
  for ir,row in listrow.iterrows():
    temp_band = lBand(listR,row["J0"],row["gamma1"],row["gamma2"],row["R0"])
    for i in range(len(temp_band)):
      dict_hu[strdate][i] += temp_band[i]

print (dict_hu.keys())
df_hu = pd.DataFrame(dict_hu)
df_hu
#

print (df_hu)
