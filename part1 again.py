# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 20:58:17 2020

@author: caitl
"""


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

R_y=9.12*10**-8
Z=11
c=R_y/Z**2
da=5
db=4
for n1 in range(3,4):
    for dn in range (1,7):
        n2=n1+dn
        wavelength=c*(n1-da)**2*(n2-db)**2/((n2-db)**2-(n1-da)**2)*10**9
        if 0 <= wavelength <=1000000:
            print("N1=",n1,"  N2=",n2, "Lambda,nm=",wavelength)
   
        #%% 
 

n1=3
n2=5

#%%
c=3*10**8
h=6.262*10**-34    
R_y=2.18*10**-18
    

for n1 in range(1,5):
    for dn in range (1,5): 
        n2=n1+dn
        wavelength= (h*c/R_y)*((n1)**2*(n2)**2/((n2)**2-(n1)**2))
        print("N1=",n1,"  N2=",n2, "Lambda,nm=",wavelength*10**9)


#print(wavelength*10**-9)

#%%%
R_y=9.12*10**-8
Z=11
c=R_y/Z**2
dp=0
ds=10


wavelengtha=c*(3-dp)**2*(5-ds)**2/((5-ds)**2-(3-dp)**2)*10**9
wavelength0=c*(3-dp)**2*(6-ds)**2/((6-ds)**2-(3-dp)**2)*10**9


def emission(wavelength, c, n1, n2):
    wavelength =c*(n1-d1)**2*(n2-d2)**2/((n2-d2)**2-(n1-d1)**2)*10**9
    return wavelength


wavelength=np.array

100*10**-9=c(1/(5-a)**2-1/(6-a)**2)
100*10**-9/c=(1/(5-a)**2-1/(6-a)**2)

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
midwavelength=wavelength[lowerlim:upperlim]  
ml_r=midwavelength[14:29]   #Selecting only the middle points of data
midwavelength =np.concatenate((midwavelength[0:14],midwavelength[29:39]), axis=None)
midmean=mean[lowerlim:upperlim]
mm_r=midmean[14:29]
midmean =np.concatenate((midmean[0:14],midmean[29:39]), axis=None)



sy = np.sqrt(midmean)

sy_data=stdev[lowerlim:upperlim]
sy_r=sy_data[14:29]
sy_data=np.concatenate((sy_data[0:14],sy_data[29:39]), axis=None)
sy=sy_data

a= 148.5339773
x0=595.84587
sigmab=3.5983822
b=0.00018921319
c=-0.09363492952

pguess = [a,x0,sigmab,b,c]   ###  A,xo,sigma ,b ,c

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
plt.scatter(ml_r, mm_r, marker=".", s=50, label="Experimental, removed", color="red")
plt.errorbar(ml_r, mm_r, sy_r, ls='none',color="red")
plt.scatter(midwavelength, midmean, marker=".", s=50, label="Experimental")
plt.errorbar(midwavelength, midmean, sy_data, ls='none')
#plt.plot(midwavelength, midmean, label = 'data')
plt.plot(testing_array, yestimate, label = 'Gaussian Fit')
plt.ylabel("Intensity")
plt.yscale("symlog")
plt.xlabel("Wavelength, nm")
plt.title("Gaussian Fit for peak around "+str(center)+"nm")
plt.legend()

print("A:", p[0],"Error in A:", np.sqrt(cov[0][0]))
print("X0:", p[1],"Error in X0:", np.sqrt(cov[1][1]))
print("Sigma:", p[2],"Error in sigma:", np.sqrt(cov[2][2]))
print("B:", p[3],"Error in B:", np.sqrt(cov[3][3]))
print("C:", p[4],"Error in C:", np.sqrt(cov[4][4]))






#%%
h=6.626*10**-34
#h=4.135667696*10**-15
c=299792458	
#n=np.array([35/6,16/3,100/21,9/2,196/45])
G_wave_n=np.array([661.99,487.46,434.14,410.20])
G_wave_s=np.array([0.0066,0.22,0.25,0.32])

n2=np.array([3,4,5,6])
n=(4*n2**2)/(n2**2-4)
stuff=(h*c*10**9)

c,cov=np.polyfit(n,G_wave_n,1, w=G_wave_s, cov=True)
a=c[0]
b=c[1]
sa=np.sqrt(cov[0][0])
sb=np.sqrt(cov[1][1])

grad=ufloat(a,sa)
inter=ufloat(b,sb)
gradient=ufloat(92.72,0.34)
print("intercept:",inter,"Gradient:",grad)
Ry=stuff/gradient
print(Ry)
plt.scatter(n, G_wave_n, marker=".", s=5, label="Experimental")
plt.errorbar(n, G_wave_n, G_wave_s,ls='none',capthick=1)
plt.plot(n,n*grad.n+inter.n, linewidth=0.5,label="Line of best fit $\pm1 \sigma$")
plt.fill_between(n,n*(grad.n+grad.s)+(inter.n+inter.s),n*(grad.n-grad.s)+(inter.n-inter.s), alpha=0.2, )
plt.xlabel(r"$(n_2^2n_1^2)/(n_2^2-n_1^2)$")
plt.ylabel("Wavelength, nm")
plt.title("Gaussian Fit")
plt.grid()
plt.legend()
