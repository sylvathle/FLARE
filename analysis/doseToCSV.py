import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
import sys,os

def extract_prefixes(lst):
    prefixes = set()
    for string in lst:
        prefix = string.split('_')[0]
        prefixes.add(prefix)
    return list(prefixes)

def extract_column_names(csv_file):
    column_names = []
    with open(csv_file, 'r') as f:
        line_num = 1
        for line in f:
            if len(line) > 0 and line[0]=='#':
    	        # Assuming the column names are in the 5th line after the symbol '#'
    	        if line_num >= 5:
    	            column_names.append(line.split()[-1])  # Last word of the line contains the column names
            line_num+=1                
    return column_names
    
def get_stats(group,list_particles):

    dict_out = {}
    #print (group)

    for p in list_particles:
      
      dict_out[p+"_N"] = group[p+'_N'].sum() 
      
      for q in ["DE","AD"]:
      
        dict_out[p+"_"+q] = (group[p+"_"+q] * group[p+'_N']).sum() / group[p+'_N'].sum()
        low_de_group = group[group[p+"_"+q]<dict_out[p+"_"+q]]

        if len(low_de_group)==0: dict_out[p+"_"+q+"_b"] = np.nan
        else: dict_out[p+"_"+q+"_b"] = np.sqrt((low_de_group[p+'_N'] * (low_de_group[p+"_"+q] - dict_out[p+"_"+q]) ** 2).sum() / low_de_group[p+'_N'].sum())

        up_de_group = group[group[p+"_"+q]>=dict_out[p+"_"+q]]
        if len(up_de_group)==0: dict_out[p+"_"+q+"_t"] = np.nan
        else: dict_out[p+"_"+q+"_t"] = np.sqrt((up_de_group[p+'_N'] * (up_de_group[p+"_"+q] - dict_out[p+"_"+q]) ** 2).sum() / up_de_group[p+'_N'].sum())

    return pd.Series(dict_out)

res_dir = "../results/ICRP145/"
res_dir = "../results/"
list_exp = [x[0] for x in os.walk(res_dir) if x[0].count("/")==3]


N_av = 5
list_av_df_doses = [pd.DataFrame() for i in range(N_av)]

for dir_path in list_exp:

  dir_path_split = dir_path.split("/")
  dir_exp = dir_path_split[-1]

  print (dir_exp)

  list_prefixes =  sorted(set([f.split('_')[0] for f in listdir(dir_path) if isfile(join(dir_path, f))]))

  for iprefix, prefix in enumerate(list_prefixes):

    dose_csv = dir_path + "/" + prefix + "_nt_Doses.csv"
    if not os.path.exists(dose_csv): continue

    cols_dose = extract_column_names(dose_csv)
    df_dose = pd.read_csv(dose_csv,names=cols_dose,skiprows=len(cols_dose)+4)
    if len(df_dose)==0: continue
    df_dose["scenario"] = dir_exp
    
    list_av_df_doses[iprefix%N_av] = pd.concat([list_av_df_doses[iprefix%N_av],df_dose])

    #break
  #break
    

#sys.exit()

df_dose = pd.DataFrame()
for df in list_av_df_doses:
  #df = df[df["proton_N"]>1000]
  #print (df[df["eBin"]==60])
  print (df)
  if len(df)==0: continue
  df_grouped = df.groupby(["scenario","organId","eBin"],as_index=False).sum()
  #print (df_grouped[df_grouped["eBin"]==60])
  for c in df_grouped.columns: 
    if "AD" in c:
      particle = c.split("_")[0]
      df_grouped[c] = df_grouped[c]/df_grouped[particle+"_N"]
      df_grouped[particle+"_DE"] = df_grouped[particle+"_DE"]/df_grouped[particle+"_N"]
  #print (df_grouped[df_grouped["eBin"]==60])
  #print ("\n")
  df_dose = pd.concat([df_dose,df_grouped])
#print (df_dose[df_dose["eBin"]==60])
#sys.exit()

list_particles = [p.split("_")[0] for p in df_dose.columns if "AD" in p]

df_dose = df_dose.groupby(by=["scenario","organId","eBin"],as_index=False).apply(lambda x: pd.concat([get_stats(x,list_particles)], axis=0))
  
df_organs = pd.read_csv("organsInfo.csv")
df_dose = df_dose.merge(df_organs[['organId', 'group', 'WT','mass[g]']], on='organId', how='left')

list_q = ["_AD","_AD_b","_AD_t","_DE","_DE_b","_DE_t"]

new = df_dose["scenario"].str.split("_",n=1,expand=True)
df_dose["scenario"] = new[0]
df_dose["thick"] = new[1]

cols_in_order = ["scenario","thick","eBin",'organId', 'group', 'WT','mass[g]']
for p in list_particles:
  cols_in_order.append(p+"_N")
for p in list_particles:
  for q in list_q:
    cols_in_order.append(p+q)

#MeVToJoules = 1.60218e-13
#for p in list_particles:
#  for q in list_q:
#    df_dose[p+q] = df_dose[p+q]/df_dose['mass[g]']*1e3*MeVToJoules

print (df_dose)


df_dose[cols_in_order].to_csv("data_dose.csv",index=False)




  
