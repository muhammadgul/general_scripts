#!/usr/bin/python
#import ROOT
import numpy as np
import matplotlib.pyplot as plt
x, y = np.loadtxt("scan_run_[04-15].txt",skiprows = 1, usecols=(1,2), unpack=True)
plt.plot(x,y,label='$M_{\chi} vs \sigma$')
plt.xlabel('$M_\chi$')
plt.ylabel('$\sigma$')
plt.show()
plt.savefig('plots/Mxd_xsec.png')
print("x:  ", x)
print("y:  ", y)
