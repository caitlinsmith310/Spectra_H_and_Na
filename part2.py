# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:49:39 2020

@author: caitl
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg
from uncertainties import ufloat
from scipy import stats
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
from scipy import signal
from scipy.interpolate import interp1d
import math
import cmath as cm
import pandas as pd
import numpy as np
from numpy import array
from scipy.optimize import curve_fit
from numpy import where
from scipy.optimize import fmin
import matplotlib.pyplot as plt
import os


#from pandas import readcsv
#%%

R_y=2.18*10**-18
c=3*10**8
h=6.626*10**-34
Z=11
const=h*c/(R_y*Z**2)
d=0
for n1 in range(1,20):
    for dn in range (1,20):
        n2=n1+dn
        wavelength=const*(n1-d)**2*(n2-d)**2/((n2-d)**2-(n1-d)**2)*10**9
        if 380 <= wavelength <=740:
            print("N1=",n1,"  N2=",n2, "Lambda,nm=",wavelength)
        
        
#%%
            
directory =r'C:\Users\caitl\Documents\390\Hydrogen'
os.chdir(directory)


df = pd.read_csv('sodium_data.csv')
df = df.values

wavelength_exp =np.array( df[:,0])
intensity_exp = df[:,1]
wavelength=np.arange(min(wavelength_exp)-1, max(wavelength_exp)+1)
#plt.plot(wavelength_exp,intensity_exp)
mean=np.array([0])
stdev=np.array([0])

for n in range(322,747):
    index=[i for i in range(0,len(wavelength_exp)) if wavelength_exp[i] == n]
    intensity_temp=intensity_exp[index]
    stdev_temp=np.std(intensity_temp)
    mean_temp=np.mean(intensity_temp)
    mean=np.append(mean, mean_temp)
    stdev=np.append(stdev,stdev_temp)
    #n=n+1
wavelength = wavelength[1:]
stdev=stdev[1:]
mean=mean [1:]
plt.figure("raw")
plt.scatter(wavelength, mean, marker=".", s=50)
plt.errorbar(wavelength, mean, stdev,ls='dotted')
plt.ylabel("Intensity")
plt.xlabel("Wavelength, nm")
plt.title("Raw data for sodium spectum")



#%%
    
def Gaussian(x, A,  x0, sigma, b, c):
    y = A*np.exp(-(x-x0)**2/sigma**2) + b*x + c  #Defining a function for one Gaussian peak
    return y



center=595
spread=20
lowerlim=center-spread-322
upperlim=center+spread-322
print("Bounds:",lowerlim+322,":",upperlim+322)
midwavelength=wavelength[lowerlim:upperlim]     #Selecting only the middle points of data
midwavelength =np.concatenate((midwavelength[0:13],midwavelength[29:39]), axis=None)
midmean=mean[lowerlim:upperlim]
midmean =np.concatenate((midmean[0:13],midmean[29:39]), axis=None)

sy = np.sqrt(midmean)

sy_data=stdev[lowerlim:upperlim]
#sy=sy_data
sy_data=np.concatenate((sy_data[0:13],sy_data[29:39]), axis=None)

a= 171.878313776
x0=596.070080468
sigma=3.48530345807
b=0
c=0.5

pguess = [a,x0,sigma,b,c]   ###  A,xo,sigma ,b ,c

testing_array = np.arange(min(midwavelength), max(midwavelength), 0.01)

#Our guess is made by varying parameters and observing the effect on the theoretical curve 
p,cov = curve_fit(Gaussian, midwavelength, midmean, sigma=sy, p0=pguess, )
yestimate = Gaussian(testing_array , *p)
#print(yestimate)
#print(p)

#plt.figure("Gaussian Fit on 515 peak")
#plt.plot(midwavelength, midmean, label = 'data')
#plt.plot(testing_array, yestimate, label = 'Gaussian Fit')
#plt.legend()


plt.figure("Gaussian Fit")
plt.scatter(midwavelength, midmean, marker=".", s=50, label="Experimental")
plt.errorbar(midwavelength, midmean, sy_data, ls='none')
#plt.plot(midwavelength, midmean, label = 'data')
plt.plot(testing_array, yestimate, label = 'Gaussian Fit')
plt.ylabel("Intensity")
plt.xlabel("Wavelength, nm")
plt.title("Gaussian Fit for peak around "+str(center)+"nm")
plt.legend()

print("A:", p[0],"Error in A:", np.sqrt(cov[0][0]))
print("X0:", p[1],"Error in X0:", np.sqrt(cov[1][1]))
print("Sigma:", p[2],"Error in sigma:", np.sqrt(cov[2][2]))
print("B:", p[3],"Error in B:", np.sqrt(cov[3][3]))
print("C:", p[4],"Error in C:", np.sqrt(cov[4][4]))






#%%
c=3*10**8
h=6.626*10**-34
R_y=2.18*10**-18

L_ob=ufloat(596.846,0.168)
L_pre=0


hi=(h*c)*((1/L_ob)-(1/L_pre))*10**9/(2*R_y*11**2)
print(hi)
