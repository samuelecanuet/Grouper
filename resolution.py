from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def linear(x, m):
    return m*x

def sqrtt(x, m, n):
    return m*np.sqrt(x) + n

def affine(x, m, n):
    return m*x + n
res_Sr = [224, 477, 1337, 2370]
res_Ar = [280, 556, 950]
calib = [4.8]
endpoint_exp_Ar = [6000]
endpoint_sim_Ar = [5625, 6000, 6300]
endpoint_exp_Sr = [3000]
endpoint_sim_Sr = [2550, 2900, 3750, 4300]

fig, ax1 = plt.subplots()
truc=4
# First plot
ax1.scatter(res_Ar, endpoint_sim_Ar, color='red', label='Ar', alpha=0.5)
ax1.scatter(res_Sr, endpoint_sim_Sr, color='blue', label='Sr', alpha=0.5)
popt, err = curve_fit(sqrtt, res_Sr, endpoint_sim_Sr)
ax1.plot(np.linspace(0, 3000, 1000), sqrtt(np.linspace(0, 3000, 1000), popt[0], popt[1]), color='blue', label='Sr', alpha=0.5)
ax1.set_ylim(0, 7000)
ax1.plot(np.linspace(0, 3000, 1000), sqrtt(np.linspace(0, 3000, 1000), popt[0], 2.65*popt[1]), color='red', label='Sr', alpha=0.5)
print(2.65*popt[1])
# Second plot with a shared y-axis and a separate x-axis
ax2 = ax1.twiny()

for i in np.linspace(1, 8, 16):
    ax2.plot(np.linspace(0, 2000, 1000), linear(np.linspace(0, 2000, 1000), 4.9), color='red', label='Ar', alpha=(i-1)/8)
    ax2.plot(np.linspace(0, 2000, 1000), linear(np.linspace(0, 2000, 1000), i), color='blue', label='Sr', alpha=(i-1)/8)
ax2.set_ylim(0, 7000)

plt.show()