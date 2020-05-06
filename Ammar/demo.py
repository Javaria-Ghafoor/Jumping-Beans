# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:09:09 2020

@author: Home

references:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
    http://kitchingroup.cheme.cmu.edu/blog/2011/08/31/Solving-an-ode-for-a-specific-solution-value/
"""

#############################
from math import *
from mpmath import *
import numpy as np 
from numpy import *
#from scipy.optimize import *
from sympy.interactive import printing
printing.init_printing(use_latex = True)

from matplotlib import*
from sympy import *
import sympy as sp
from scipy import *

from scipy.integrate import odeint 
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

#############################


############## data ######################



L = 0.0126       # length of capsule
#R = 0.0024       # radius of ball bearing
#Lb = 0.9         # length of board
#m = 0.0004541    # mass of ball
#M = 0.0000624    # mass of capsule
#meu_1 = 0.032    # put meu_1 and meu_2 = 0 for no friction / ideal case solution
#meu_2 = 0.34
alpha = np.pi/4 
g = 9.8 
#
#A_cyl = 2 * pi * R * (L - 2 * R)      # surface area of clyinder 
#A_ss = 2 * pi * R ** 2                # surface area hemisphere
#sigma = M / (A_cyl + 2 * A_ss)        # mass density 
#
#I_cyl = sigma * A_cyl * (1/4 * R ** 2 + 1/3 * (L - 2 * R)**2)  
#I_lss = sigma * A_ss * (1/12 * R ** 2 + (1/2 * R + (L - 2 * R))**2)
#I_rss = sigma * A_ss * (1/3 * R ** 2)
#I_ball = 2/5 * m * R ** 2
#
#I_cap = I_cyl + I_lss + I_rss
#I_sys = I_ball + I_cap
#
#
#d = (M * L/2 + m * R)/(M + m)  # center of mass

############################################

def pend(y, t):
     g = 9.8
     x, v = y
     dydt = [v, g*np.sin(alpha)/3]
     return dydt

y0 = [0.0, 0.0]

#We will generate a solution at 100 evenly spaced samples in the interval 0 <= t <= 1. So our array of times is:
t = np.linspace(0, 1, 100)

sol = odeint(pend, y0, t)


plt.plot(t, sol[:, 0], 'b', label='x(t)')
#plt.plot(t, sol[:, 1], 'g', label='v(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()

y_sol = sol[:, 0]
v_sol = sol[:, 1]
t_func = interp1d(y_sol, t, 'cubic')
v_func = interp1d(t, v_sol , 'cubic')


print("time taken to cover lenght L = ",t_func(L))




