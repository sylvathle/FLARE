import pandas as pd
import datetime
from os import listdir
from os.path import isfile, join


directory = "output/SPE/tserie/"
#list_files = [f for f in listdir(directory) if isfile(join(directory, f))]

#print (list_files)

dict_df = {}
for f in listdir(directory):
  fsplit = f.split("_")
  scenario = fsplit[0]
  organ = fsplit[1]
  k  = scenario+"_"+organ
  if k not in dict_df.keys(): 
    dict_df[k]=pd.DataFrame()
  str_dates = fsplit[2].split("-")
  str_datemin = str_dates[0]
  datemin = datetime.datetime.strptime(str_datemin,"%Y%m%d")
  try: df_k = pd.read_csv(directory+f)
  except: continue
  dict_df[k] = pd.concat([dict_df[k],df_k])

for k,df in dict_df.items():
  print (k)
  if len(df)==0: continue
  df["date"] = pd.to_datetime(df["date"])
  df.sort_values(by=["date"],inplace=True)
  df.to_csv("output/SPE/final_tserie/"+k+".csv",index=False)
