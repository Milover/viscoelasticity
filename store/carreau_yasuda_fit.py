#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


#xp = np.linspace(10, 100, 100)                                                              # domain
#x = np.array([50.25, 67.0, 83.75, 100.5, 115.625, 140.0, 164.375, 188.75, 237.5, 286.25, 335.0])         # actual
#x = np.divide(x, 335.0 / 100.0)                                                                     # scale
#y = np.array([354.2439, 149.0270, 72.1838, 39.9015, 26.3725, 13.8310, 8.4640, 5.7024, 2.8352, 1.6831, 1.1015])

xp = np.logspace(-3, 4, 1000)
x = np.array([5.704246e-02,
              2.715733e-01,
              1.268560e+00,
              6.039422e+00,
              2.821135e+01,
              7.032593e+01,
              9.534940e+01,
              1.268550e+02,
              4.992838e+02,
              6.038903e+02])
y = np.array([1.972667e+02,
              1.120654e+02,
              5.110685e+01,
              2.407842e+01,
              1.276425e+01,
              9.143931e+00,
              8.293010e+00,
              7.645302e+00,
              6.603757e+00,
              6.549989e+00])
y = y / 1000.0      # y is in cP = mPa.s, we need Pa.s

# p[0] - mu_inf
# p[1] - mu_0
# p[2] - k
# p[3] - a
# p[4] - n

def func(x, *p):
    return p[0] + (p[1] - p[0])*np.float_power(1.0 + np.float_power(p[2]*x, p[3]), ((p[4] - 1)/p[3]))

popt, pcov = curve_fit(func, x, y, p0=np.array([0.0065,0.2,1.9,1.25,0.22]))


plt.grid(True, 'both')
plt.plot(x, y, 'k.', label='data', markersize=8)

plt.plot(xp, func(xp, *popt), 'r-', linewidth=1.5,
         label=r'$\eta_\infty + (\eta_0 - \eta_\infty) \left[1 + (k \dot{\gamma})^a\right]^\frac{n - 1}{a}$'
                '\n'
               r'$\eta_\infty=%5.3f, \eta_0=%5.3f, k=%5.3f, a=%5.3f, n=%5.3f$' % tuple(popt))


# residuals
print(np.linalg.norm(y - func(x, *popt)))

# coeffs
print(tuple(popt))

plt.xscale('log')
plt.xlabel('$\dot{\gamma}$, [s$^{-1}$]')
plt.xlim(0.01, 1000)
#plt.locator_params(nbins=15, axis='x')

plt.yscale('log')
plt.ylabel('$\eta$, [Pa s]')
plt.ylim(0.001, 1)
#plt.locator_params(nbins=25, axis='y')

plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0))
plt.show()

