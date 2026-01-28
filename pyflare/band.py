import numpy as np
import pandas as pd
import physics

def BandFunction(R,J0,gamma1,gamma2,R0):
  if R<=(gamma2-gamma1)*R0:
    return float(J0*R**(-gamma1)*np.exp(-R/R0))
  else:
    diffgam = gamma2-gamma1
    return float(J0*R**(-gamma2)*(diffgam*R0)**diffgam*np.exp(-diffgam))

def logBand(R,J0,gamma1,gamma2,R0):
  lreturn = []
  for r in R:
    if r<=(gamma2-gamma1)*R0:
      lreturn.append(np.log10(float(J0*r**(-gamma1)*np.exp(-r/R0))))
    else:
      diffgam = gamma2-gamma1
      #print ("pars",J0,gamma1,gamma2,R0,diffgam)
      #print ("exps",diffgam,r**(-gamma2),(diffgam*R0)**diffgam,np.exp(-diffgam))
      lreturn.append(np.log10(float(J0*r**(-gamma2)*(diffgam*R0)**diffgam*np.exp(-diffgam))))

  return lreturn

def lBand(R,J0,gamma1,gamma2,R0):
  lreturn = []
  for r in R:
    if r<=(gamma2-gamma1)*R0:
      lreturn.append(float(J0*r**(-gamma1)*np.exp(-r/R0)))
    else:
      diffgam = gamma2-gamma1
      lreturn.append(float(J0*r**(-gamma2)*(diffgam*R0)**diffgam*np.exp(-diffgam)))
  return lreturn

def diffBand(R,J0,gamma1,gamma2,R0):
  lreturn = []
  for r in R:
    if r<=(gamma2-gamma1)*R0:
      dJdr = -J0*r**(-gamma1)*np.exp(-r/R0)* (gamma1/r + 1.0/R0)
      lreturn.append(dJdr)
    else:
      diffgam = gamma2-gamma1
      dJdr = -gamma2*J0*r**(-gamma2-1)*(diffgam*R0)**diffgam*np.exp(-diffgam)
      lreturn.append(dJdr)

  return lreturn

def quadratic_interpolate_series(series, value):                                                                                                                       # Check if the value is within the range of the index
    if value < series.index.min() or value > series.index.max():
        return np.nan

    # Find the indices of the three closest points
    index = series.index
    lower_index = index[index <= value].max()
    upper_index = index[index >= value].min()

    if lower_index == upper_index:
        return series[lower_index]

    # Find the third point
    if lower_index == index.min():
        third_index = index[index > upper_index].min()
        points = [lower_index, upper_index, third_index]
    elif upper_index == index.max():
        third_index = index[index < lower_index].max()
        points = [third_index, lower_index, upper_index]
    else:
        third_index_lower = index[index < lower_index].max()
        third_index_upper = index[index > upper_index].min()
        if abs(value - third_index_lower) < abs(value - third_index_upper):
            points = [third_index_lower, lower_index, upper_index]
        else:
            points = [lower_index, upper_index, third_index_upper]

    # Get the corresponding values
    x = np.array(points)
    y = np.array([series[p] for p in points])

    # Perform quadratic interpolation
    coefficients = np.polyfit(x, y, 2)
    polynomial = np.poly1d(coefficients)
    interpolated_value = polynomial(value)

    return interpolated_value

def getCumul(flist):
  cumlist = [flist[-1]]
  for i in range(len(flist)-2,-1,-1):
    cumlist.append(cumlist[-1]+flist[i])
  return cumlist[::-1]

# Convert integral to differential spectrum
def getBandSpectrum(listR,J0,gamma1,gamma2,R0):
  bandFunction = []
  eplot = []
  rplot = []
  for ie in range(len(listR)-1):
    Rup = listR[ie+1]#*1.0001
    Rdown = listR[ie]#*0.9999
    #Eup = rigidity2Energy(Rup,mass_proton,1)
    #Edown = rigidity2Energy(Rdown,mass_proton,1)

    dbandup = 0
    dbanddown = 0

    dbandup += BandFunction(Rup,J0,gamma1,gamma2,R0)
    dbanddown += BandFunction(Rdown,J0,gamma1,gamma2,R0)
    dband = -dbandup+dbanddown
    #print (dband)
    # Not allowing negative values, unphysical
    if dband<0:
      bandFunction.append(float('NaN'))
    else:
      #dband /= (listE[ie+1]-listE[ie])
      #dband /= Eup-Edown
      #eplot.append(listE[ie])
      #rplot.append(listR[ie])
      bandFunction.append(dband)

  return bandFunction#,eplot

# Convert integral to differential spectrum
def getBandDiffSpetrum(listR,J0,gamma1,gamma2,R0):
  bandFunction = []
  eplot = []
  rplot = []
  for ie in range(len(listR)):
    Rup = listR[ie]*1.0001
    Rdown = listR[ie]*0.9999
    Eup = rigidity2Energy(Rup,mass_proton,1)
    Edown = rigidity2Energy(Rdown,mass_proton,1)

    dbandup = 0
    dbanddown = 0

    dbandup += Band(Rup,J0,gamma1,gamma2,R0)
    dbanddown += Band(Rdown,J0,gamma1,gamma2,R0)
    dband = -dbandup+dbanddown
    #print (dband)
    # Not allowing negative values, unphysical
    if dband<0:
      bandFunction.append(float('NaN'))
    else:
      #dband /= (listE[ie+1]-listE[ie])
      dband /= Eup-Edown
      #eplot.append(listE[ie])
      #rplot.append(listR[ie])
      bandFunction.append(dband)

  return bandFunction#,eplot


class Band:
  def __init__(self):
    self.band_path = 'input/band_parameters.csv'

    self.df_band = pd.read_csv(self.band_path)

    self.df_band["Date"] = pd.to_datetime(self.df_band["Date"])
    self.df_band["Date"] = self.df_band["Date"].dt.floor('D')
    self.df_band["Date"] = self.df_band["Date"].dt.strftime("%Y%m%d")
    self.df_band["end date"] = pd.to_datetime(self.df_band["end date"])

    self.df_band = self.df_band[self.df_band["GLE"]!="24alt"]

    self.df_band["GLE"] = self.df_band["GLE"].astype(int)
    self.list_sources = self.df_band["Source"].unique().tolist()
    self.list_spe = self.df_band["Date"].unique().tolist()

  def getDatesSPE(self,source):
    #print (self.df_band[self.df_band["Source"]==source])
    return self.df_band[self.df_band["Source"]==source]["Date"].unique().tolist()

  def getSPEspectrum(self,source,date,listEp1):
    df_spe = self.df_band[self.df_band["Source"]==source]
    df_spe = df_spe[df_spe["Date"]==date]

    if len(df_spe)==0: return None
  
    df_spectrum = pd.DataFrame()
    df_spectrum["E"] = listEp1[:-1]
    
    datestr = date
    df_spectrum[datestr] = 0

    for irow,row in df_spe.iterrows():
      listRp1 = [physics.rigidity(e,physics.mass_proton,1) for e in listEp1]
      temp_list = getBandSpectrum(listRp1,row["J0"],row["gamma1"],row["gamma2"],row["R0"])
      df_spectrum[datestr] += temp_list

    df_spectrum.set_index("E",inplace=True)
    return df_spectrum

  def getListSources(self):
    return self.list_sources
    
  def getListSPE(self):
    return self.list_spe
