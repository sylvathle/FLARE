import pandas as pd
import math
import datetime

list_cols = ["Source","GLE","Date","end date","J0","gamma1","gamma2","R0"]

df_sepem_hulim = pd.read_csv("band_sepem_withHUlims_fit.csv")
df_sepem_hulim["Source"] = "sepem_hulim"
df_sepem_hulim["GLE"] = -1
df_sepem_hulim["Date"] = pd.to_datetime(df_sepem_hulim["Date"])
df_sepem_hulim["end date"] = pd.to_datetime(df_sepem_hulim["end date"])
df_sepem_hulim = df_sepem_hulim[list_cols]

df_sepem = pd.read_csv("band_sepem_fit.csv")
df_sepem["Source"] = "sepem"
df_sepem["GLE"] = -1
df_sepem["Date"] = pd.to_datetime(df_sepem["Date"])
df_sepem["end date"] = pd.to_datetime(df_sepem["end date"])
df_sepem = df_sepem[list_cols]

# Get Storms from sepem that are before the Hu windows
df_sepem1 = df_sepem[df_sepem["Date"]<min(df_sepem_hulim["Date"])]

# Get Storms from sepem that are after the Hu windows
df_sepem2 = df_sepem[df_sepem["Date"]>max(df_sepem_hulim["end date"])]

df_tylka_hu = pd.read_csv("data/band_tylka-hu.csv")
df_tylka_hu["Date"] = pd.to_datetime(df_tylka_hu["Date"])
enddate = []
for i,row in df_tylka_hu.iterrows():
  if math.isnan(row["Duration"]):
    enddate.append(None)
  else:
    enddate.append(row["Date"]+datetime.timedelta(hours=row["Duration"]))
df_tylka_hu["end date"] = enddate
df_tylka_hu = df_tylka_hu[["Date","end date","Source","GLE","J0","gamma1","gamma2","R0"]]
print (df_tylka_hu[list_cols])

df_goes = pd.read_csv("data/Goes/band_goes.csv")
df_goes["Date"] = pd.to_datetime(df_goes["Date"])
df_goes["end date"] = pd.to_datetime(df_goes["end date"])
df_goes["GLE"] = -1
df_goes.rename(columns={"goes":"Source"},inplace=True)

print (df_goes[list_cols])

df_sep = pd.concat([df_tylka_hu,df_sepem,df_sepem_hulim,df_goes])
df_sep = df_sep.sort_values(by=["Date"])
df_sep.to_csv("band_parameters.csv",index=False)
