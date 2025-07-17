#import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
#import datetime
import math
from scipy.optimize import curve_fit
import sys, os


def lpolynomial(R,a,b,c,d):
  lreturn = []
  for r in R:
    lr = np.log10(r)
    lreturn.append(10**(float(a*lr**3+b*lr**2+c*lr+d)))
  return lreturn

def lBand(R,J0,gamma1,gamma2,R0):
  lreturn = []
  for r in R:
    if r<=(gamma2-gamma1)*R0:
      lreturn.append(float(J0*r**(-gamma1)*np.exp(-r/R0)))
    else:
      diffgam = gamma2-gamma1
      #print ("pars",J0,gamma1,gamma2,R0,diffgam)
      #print ("exps",diffgam,r**(-gamma2),(diffgam*R0)**diffgam,np.exp(-diffgam))
      lreturn.append(float(J0*r**(-gamma2)*(diffgam*R0)**diffgam*np.exp(-diffgam)))

  return lreturn

def logpolynomial(R,a,b,c,d):
  lreturn = []
  for r in R:
    lr = np.log10(r)
    lreturn.append(float(a*lr**3+b*lr**2+c*lr+d))
  return lreturn

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

def dBand(R,J0,gamma1,gamma2,R0):
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

def powerLaw(E,a,b): return a*E**b

def logpowerLaw(E,a,b): return np.log10(a*E**b)


# Evaluation of the fit

def mse_band(R,vals,function,popt):
  fit = function(R,*popt)
  err = 0
  for i in range(len(vals)):
    err += 0.5*abs(vals[i]-fit[i])/abs(vals[i]+fit[i])
  err = err/len(vals)
  return err

def mse(Y1,Y2):
  err = 0
  for i in range(len(Y1)):
    err += 0.5*abs(Y1[i]-Y2[i])/abs(Y1[i]+Y2[i])
  err = err/len(Y1)
  return err
