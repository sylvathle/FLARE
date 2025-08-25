import urllib.request
import os.path
import pandas as pd


def getSunSpots(forceDownload=False):

  ssfile = "input/SunSpots.csv"
  if os.path.isfile(ssfile):
    urllib.request.urlretrieve("https://www.sidc.be/SILSO/INFO/sndtotcsv.php",ssfile)

  df_ss_tot = pd.read_csv(ssfile,sep=";")
  cols = ['year','month','day','Dec_date','SS','std','n_obs','indicator']
  df_ss_tot.columns = cols

  df_ss_tot['date'] = pd.to_datetime(df_ss_tot[['year', 'month', 'day']])

  # Set the new datetime column as the index
  df_ss_tot.set_index('date', inplace=True)
  return df_ss_tot

print (getSunSpots())
