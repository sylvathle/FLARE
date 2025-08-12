import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
import sys,os

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
    
def get_stats(group,vars_list):

    dict_out = {}
    #print (group)

    #for p in list_particles:
    if True:
      
      dict_out["N"] = group['N'].sum() 
      
      for q in vars_list:
      
        dict_out[q] = (group[q] * group['N']).sum() / group['N'].sum()
        low_de_group = group[group[q]<dict_out[q]]

        if len(low_de_group)==0: dict_out[q+"_b"] = np.nan
        else: dict_out[q+"_b"] = np.sqrt((low_de_group['N'] * (low_de_group[q] - dict_out[q]) ** 2).sum() / low_de_group['N'].sum())

        up_de_group = group[group[q]>=dict_out[q]]
        if len(up_de_group)==0: dict_out[q+"_t"] = np.nan
        else: dict_out[q+"_t"] = np.sqrt((up_de_group['N'] * (up_de_group[q] - dict_out[q]) ** 2).sum() / up_de_group['N'].sum())

    return pd.Series(dict_out)

#res_dir = "../results/IcruSphere/rego30/"
res_dir = "../results/"

df_all_dose = pd.DataFrame()
#df_all_dose_sample = pd.DataFrame()
df_all_flux = pd.DataFrame()

#i_sample = -1

for scenario in listdir(res_dir):
  print (scenario)
  #continue
  #i_sample = i_sample+1
  scenario_dir = res_dir + scenario + "/"

  list_prefix = [scenario_dir+f.split("_")[0]+"_nt_" for f in listdir(scenario_dir) if isfile(join(scenario_dir, f)) and "Dose" in f]

  N_av = 5
  list_av_df_doses = [pd.DataFrame() for i in range(N_av)]
  #list_av_df_flux = [pd.DataFrame() for i in range(N_av)]
  list_sum_df_N = [pd.DataFrame() for i in range(N_av)]

  for iprefix, f in enumerate(list_prefix):
    cols_dose = extract_column_names(f+"Doses.csv")
    df_dose = pd.read_csv(f+"Doses.csv",names=cols_dose,skiprows=len(cols_dose)+4)
    #cols_flux = extract_column_names(f+"InnerFlux.csv") 
    #df_flux = pd.read_csv(f+"InnerFlux.csv",names=cols_flux,skiprows=len(cols_flux)+4)
    cols_N = extract_column_names(f+"N.csv") 
    df_N = pd.read_csv(f+"N.csv",names=cols_N,skiprows=len(cols_N)+4)

  #df_dose = df_dose.merge(df_N,left_on="eBin",right_on="ikE")
  #df_dose.drop("ikE",axis=1,inplace=True)


  #df_flux = df_flux.merge(df_N,left_on="ikE",right_on="ikE")

    if len(df_dose)==0 or len(df_N)==0: continue
  
    list_av_df_doses[iprefix%N_av] = pd.concat([list_av_df_doses[iprefix%N_av],df_dose])
    #list_av_df_flux[iprefix%N_av] = pd.concat([list_av_df_flux[iprefix%N_av],df_flux])
    list_sum_df_N[iprefix%N_av] = pd.concat([list_sum_df_N[iprefix%N_av],df_N])

  # Sum N by group of average
  df_N = pd.DataFrame()
  for i in range(len(list_sum_df_N)):
    df = list_sum_df_N[i]
    if len(df)==0: continue
    list_sum_df_N[i] = df.groupby(["ikE"],as_index=False).sum()
    #print (list_sum_df_N[i])

  ## Divide values by N and get mean, std_up and std_down for doses
  df_dose = pd.DataFrame()
  for i,df in enumerate(list_av_df_doses):
    #print (df)
    if len(df)==0: continue
    df_grouped = df.groupby(["eBin","organId"],as_index=False).sum()
    df_grouped = df_grouped.merge(list_sum_df_N[i],left_on="eBin",right_on="ikE")
    df_grouped.drop("ikE",axis=1,inplace=True)
    #print (df_grouped)

    df_grouped["DE"] = df_grouped["DE"] / df_grouped["N"]
    df_grouped["AD"] = df_grouped["AD"] / df_grouped["N"]
    df_grouped["i_sample"] = i

    df_dose = pd.concat([df_dose,df_grouped])

  var_list = ["DE","AD"]
  df_dose["scenario"] = scenario
  print (df_dose)
  #df_all_dose_sample = pd.concat([df_all_dose_sample,df_dose])
  df_dose_av = df_dose.groupby(by=["eBin","scenario","organId"],as_index=False).apply(lambda x: pd.concat([get_stats(x,var_list)], axis=0))
  df_all_dose = pd.concat([df_all_dose,df_dose_av])

  ## Divide values by N and get mean, std_up and std_down for fluxes
  #df_flux = pd.DataFrame()
  #for i,df in enumerate(list_av_df_flux):
  #  if len(df)==0: continue
  #  df_grouped = df.groupby(["ikE","okE"],as_index=False).sum()
  #  df_grouped = df_grouped.merge(list_sum_df_N[i],left_on="ikE",right_on="ikE")
#
#    df_grouped["down"] = df_grouped["down"] / df_grouped["N"]
#    df_grouped["up"] = df_grouped["up"] / df_grouped["N"]
#
#    df_flux = pd.concat([df_flux,df_grouped])

  #var_list = ["up","down"]
  #df_flux = df_flux.groupby(by=["ikE","okE"],as_index=False).apply(lambda x: pd.concat([get_stats(x,var_list)], axis=0))
  #df_flux.drop("N",axis=1,inplace=True)

  #df_flux.sort_values(by=["ikE","okE"],ascending=True,inplace=True)
  #df_flux["scenario"] = scenario
  #df_all_flux = pd.concat([df_all_flux,df_flux])
  #df_flux.to_csv("flux.csv",index=False)

#df_all_dose_sample = df_all_dose_sample[["scenario","i_sample","eBin","N","DE","AD"]]
#df_all_dose_sample[['scenario', 'thick', 'particle']] = df_all_dose_sample['scenario'].str.split('_', expand=True)
#df_all_dose_sample.to_csv("doses_samples.csv",index=False)

df_all_dose[['scenario', 'thick', 'particle']] = df_all_dose['scenario'].str.split('_', expand=True)
df_all_dose = df_all_dose[["scenario",'thick','particle',"organId","eBin","DE","DE_b","DE_t","AD","AD_b","AD_t"]]
df_all_dose.to_csv("doses.csv",index=False)

#df_all_flux = df_all_flux[["scenario","ikE","okE","up","up_b","up_t","down","down_b","down_t"]]
#df_all_flux.to_csv("flux.csv",index=False)

sys.exit()
