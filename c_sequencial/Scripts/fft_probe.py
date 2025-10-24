import matplotlib.pyplot as plt
#import plotly.plotly as py
import numpy as np
from scipy.signal import argrelextrema
import scipy.fftpack
from spf import *
# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

Lc_factor = 1.0

Nx ,Ny , y0 = np.loadtxt('data/mesh.dat', skiprows=0, unpack=True)

Nx = int(Nx)
Ny = int(Ny)
y0 = int(y0)

xc,yc = np.loadtxt('transient/data/grid.dat', skiprows=0, unpack=True)

xc = np.array(xc[0: Nx])* Lc_factor
yc = np.array(yc[0::Nx])* Lc_factor


N = int(input('enter the final iteration:  '))
O = int(input('enter the step:  '))

M=102
N=N+102
NM=int(int(N-M)/O)

u =     ['']*NM
v =     ['']*NM
speed = ['']*NM
C =     ['']*NM
P =     ['']*NM
T =     ['']*NM
Xf =    ['']*NM
Yf =    ['']*NM
print ('data size:', NM)
t = np.zeros(NM)# ['']*N
y = np.zeros(NM)# ['']*N

height = float(input('enter the y location of the probe:  '))

for k in range (0,NM*O,O):
	i = int(k/O)
	print (i)
	u[i] = np.loadtxt('transient/data/u'+str(k+M)+'.dat', skiprows=0, unpack=True)
	v[i] = np.loadtxt('transient/data/v'+str(k+M)+'.dat', skiprows=0, unpack=True)
	P[i] = np.loadtxt('transient/data/P'+str(k+M)+'.dat', skiprows=0, unpack=True)
	T[i] = np.loadtxt('transient/data/T'+str(k+M)+'.dat', skiprows=0, unpack=True)
	C[i] = np.loadtxt('transient/data/C'+str(k+M)+'.dat', skiprows=0, unpack=True)
	t[i] = np.loadtxt('transient/data/time'+str(k+M)+'.dat', skiprows=0, unpack=True)
#	u[i],v[i],Z[i],P[i],T[i] = np.loadtxt('transient/data/plot_field'+str(k+M)+'.dat', skiprows=0, unpack=True)
	t[i] = np.loadtxt('transient/data/time'+str(k+M)+'.dat', skiprows=0, unpack=True)
	v[i] = [iter(v[i])]*Nx
	T[i] = zip(*[iter(T[i])]*Nx)
	P[i] = zip(*[iter(P[i])]*Nx)
	v[i] = np.array(v[i])
	P[i] = np.array(P[i])
	T[i] = np.array(P[i])
	t[i] = np.array(t[i])
        #for j in range (1,len(xc),1):
        #    if xc[j] >= 0.0 and xc[j-1] < 0.0:
        #        for l in range (1,len(yc),1):
        #            if yc[l] > 0.0:
        #                for m in range (l,len(yc),1):
        #                    if v[i][m,j] >= 0.0 and v[i][m-1,j] < 0.0:
        #                        y[i] = yc[m]
	for j in range (1,len(xc)-1,1):
		if xc[j] >= -0.5 and xc[j-1] < -0.5:
			for l in range (1,len(yc),1):
				if yc[l-1] <= height:
					if yc[l] > height:
						y[i] = v[i][l,j]


						col_format = "{:<16}" * 2 + "\n"   # 2 left-justfied columns with 16 character width

with open("data/fft_probe.dat", 'w') as of:
    for x in zip(t, y):
        of.write(col_format.format(*x))

v_c = 0.6
l_c = 5.2905e-2 #*(S + 1)
t_c = l_c / v_c # 6.5443e-2
t, y = np.loadtxt('data/fft_probe.dat', skiprows=0, unpack=True)
t = t * t_c
t = t - t[0]

plt.figure()
plt.plot(t,y,'k^-',fillstyle='none')
plt.xlabel('$t$')
plt.ylabel('Propertie')
SPF()
plt.show()
