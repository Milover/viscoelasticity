#!/usr/bin/env python3

#------------------------------------------------------------------------------

import csv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#------------------------------------------------------------------------------
# Constants

diameter    = 2.0           # pipe diameter
density     = 1060.0        # density

#------------------------------------------------------------------------------
# Functions

def usage():
    print('\nUsage: effective_viscosity.py [FILE]\n')

def decomment(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield row

# linear fit
def func_l(x, *p):
    return (p[0] + p[1] * x)

# quadratic fit
def func_q(x, *p):
    return p[0] + p[1] * x + p[2] * pow(x, 2)

#------------------------------------------------------------------------------
# Work

if len(sys.argv) != 2:
    usage()
    exit()

# process data
filename = os.path.abspath(sys.argv[1])

re = np.array([], dtype=float)
v  = np.array([], dtype=float)
mu = np.array([], dtype=float)

with open(filename, 'r', newline='') as csvfile:
    reader = csv.reader(decomment(csvfile), delimiter='\t')

    for row in reader:
        re = np.append(re, float(row[0]))
        v = np.append(v, float(row[1]))
        mu = np.append(mu, density*diameter*v[-1]/re[-1])

xp = np.linspace(re[0], re[-1], 100)
popt, pcov = curve_fit(func_q, re, mu, p0=np.ones(3))

for p in tuple(popt):
    print("{:.9e}".format(p))

plt.plot(re, mu, 'k.', label='data', markersize=8)
plt.plot(xp, func_q(xp, *popt), 'r-',
         label=r'$p_0 + p_1 x + p_2 x^2$'
                '\n'
               r'$p_0=%5.9f, p_1=%5.9f, p_2=%5.9f$' % tuple(popt))
plt.grid(True, 'both')
plt.xlabel('Re')
plt.ylabel('$\mu$, Pa s')
#plt.legend(loc='upper right', bbox_to_anchor=(0.5, 1.0))
plt.legend(loc='upper right')
plt.show()

#------------------------------------------------------------------------------
