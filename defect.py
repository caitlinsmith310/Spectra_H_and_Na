# -*- coding: utf-8 -*-
"""
Created on Mon May 11 12:58:12 2020

@author: caitl
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt
import os
from uncertainties import ufloat
from uncertainties import unumpy
import scipy.linalg as linalg
from scipy.optimize import curve_fit
from matplotlib.patches import Polygon
from scipy.integrate import quad
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
import numpy as np
from numpy import linalg
import os
from uncertainties import ufloat
from uncertainties import unumpy
from scipy.optimize import curve_fit
#%%

c=299792458
h=6.626*10**-34
R_y=2.18*10**-18
Z=11
R_eff=R_y

L_obs_n=np.array([498.029,515.095,568.498,615.976])*10**-9
L_obs_error=np.array([0.013,0.019,0.027,  0.312])*10**-9
L_obs=unumpy.uarray(L_obs_n, L_obs_error)
n=np.array([5,6,4,5])


defect=n-((1/(2.117**2)-(h*c/(R_eff*L_obs))))**-0.5
print(defect)
ave=(defect[1]+defect[3])/2
print(ave)
print(R_y*(1/((3-ave)**2))/1.6e-19)