import matplotlib.pyplot as plt
#import plotly.plotly as py
import numpy as np
from scipy.signal import argrelextrema
import scipy.fftpack
from spf import *
# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

Lc_factor = 1.0

v_c = 1.
l_c = 1.e-1 #*(S + 1)

t_c = l_c / v_c # 6.5443e-2

t, y = np.loadtxt('data/fft_probe.dat', skiprows=100, unpack=True)

#t = t * t_c

t = t - t[0]

N = len(t)

T = t[N-1]/N

yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T),N/2)


fig, ax = plt.subplots(2, 1)
#ax[0].plot(t,y,'ko-', linewidth = 1.,markersize=1.5)
ax[0].plot(t,y,'k-', linewidth = 1.,markersize=1.5)
ax[0].set_xlabel('$t$')
ax[0].set_ylabel('$v$')
ax[1].set_xticks(np.arange(0,50,10))
ax[1].plot(xf[1:],2.0 / N * np.abs(yf[:N//2])[1:],'k-', linewidth = 1.) # plotting the spectrum
ax[1].set_xlabel('Strouhal Number St')
ax[1].set_ylabel('Power Spectrum')
ax[1].set_xlim(0.0,2.0)
ax[1].set_xticks(np.arange(0,2.5,0.5))
SPF()
plt.savefig('output/fft.eps')
plt.show()

a = 2.0 / N * np.abs(yf[:N//2])[1:]

for i in range (0,len(xf[1:]),1):
    if a[i] == np.amax(a):
        print 'St=',xf[i+1]
        print 'f=',xf[i+1]*t_c

#fig, ax = plt.subplots(2, 1)
#ax[0].plot(t,y,'k-', linewidth = 1.)
#ax[0].set_xlabel('Time')
#ax[0].set_ylabel('v')
#ax[1].plot(xf[1:],2.0 / N * np.abs(yf[:N//2])[1:],'k-', linewidth = 1.) # plotting the spectrum
#ax[1].set_xlabel('$St$')
#ax[1].set_ylabel('$|Y(St)|$')
#ax[1].set_xlim(0.0,10.0)
#plt.savefig('output/fft.png', bbox_inches='tight')
#plt.show()
