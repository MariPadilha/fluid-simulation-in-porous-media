import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import scipy.fftpack


H_probe = 2.5 #altura da sonda no dominio em y, em x está no centro (imax/2)
tc = 4.44e-2 # tempo caracteristico [s]

def SPF():
	F = plt.gcf()
	A = plt.gca()
	
	F.set_size_inches(6.26,1.26)
	#F.set_size_inches(2.55906,2.28012246)
	
	L = plt.legend(loc='best')
	L.get_frame().set_alpha(0.0)
	L.get_frame().set_edgecolor('w')
	
	A.tick_params(axis='y', which='minor',direction='in')
	A.tick_params(axis='x', which='minor',direction='in')
	A.tick_params(axis='both',direction='in')
	plt.tight_layout(pad=0.15)

Mi =  input('enter initial iteration (this is the initial physical time step):')

M=int(Mi)+100
N = input('enter final iteration (number of physical time steps for probing):')
N = int(N)
N= N + 101
O=1
#x = ['']*N
#y = ['']*N
u = ['']*N
v = ['']*N
speed = ['']*N
C = ['']*N
P = ['']*N
T = ['']*N
rho = ['']*N
Xf = ['']*N
Yf = ['']*N
t = ['']*N


x = np.loadtxt('data/x.dat', skiprows=0, unpack=True)
y = np.loadtxt('data/y.dat', skiprows=0, unpack=True)

probe = []
probe_time = []


count = 100
count = int(count)
print('H_probe=',H_probe)
for i in range (M,N,O):
        count = count + 1
        print(i-100)
        u[i] = np.loadtxt('transient/data/u'+str(i)+'.dat', skiprows=0, unpack=True)
        v[i] = np.loadtxt('transient/data/v'+str(i)+'.dat', skiprows=0, unpack=True)
        P[i] = np.loadtxt('transient/data/P'+str(i)+'.dat', skiprows=0, unpack=True)
        T[i] = np.loadtxt('transient/data/T'+str(i)+'.dat', skiprows=0, unpack=True)
        C[i] = np.loadtxt('transient/data/C'+str(i)+'.dat', skiprows=0, unpack=True)
        t[i] = np.loadtxt('transient/data/time'+str(i)+'.dat', skiprows=0, unpack=True)
        u[i] = np.array(u[i])
        v[i] = np.array(v[i])
        T[i] = np.array(T[i])
        C[i] = np.array(C[i])
        for k in range (1,len(y)-1,1):
                if y[k] > H_probe:
                    T2 = v[i][k,int(len(x)/2)]   #MUDAR A VARIAVEL DE INTERESSE AQUI
                    T1 = v[i][k-1,int(len(x)/2)] #MUDAR A VARIAVEL DE INTERESSE AQUI
                    y2 = y[k]
                    y1 = y[k-1]
                    Tint = T2 - ((y2-H_probe)/(y2-y1))*(T1-T2)
                    probe.append(Tint)
                    break
        probe_time.append(t[i])
print('Probe location:')
print('y=',y[k])
print('x=',x[int(len(x)/2)])

probe = np.array(probe)
probe_time = np.array(probe_time)

t = probe_time#*tc
y = probe

t = t - t[0]

N = len(t)

T = t[N-1]/N

yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T),int(N/2))


fig, ax = plt.subplots(1, 2)
ax[0].plot(t*tc,y,'k-', linewidth = 1.,markersize=3.5,fillstyle='none')
ax[0].set_xlabel('$t \ [s]$')
ax[0].set_ylabel('$v$')
#ax[0].set_yticks(np.arange(-0.20,-0.12,0.01))
#ax[0].set_xticks(np.arange(-10.0,40,5))
ax[0].set_xlim(np.min(t*tc),np.max(t*tc))
#ax[0].set_ylim(-0.17,-0.135)
ax[0].text(27., 0.09,'(a)',size=11,bbox=dict(facecolor='w',edgecolor='none', alpha=0.5,pad=0.5))    #axs[k].text(3.5, -9, letras[i],size=10,bbox=dict(facecolor='white',edgecolor='none', alpha=1.,pad=0.5))
ax[1].plot(xf[1:]/tc,2.0 / N * np.abs(yf[:N//2])[1:],'k-', linewidth = 1.) # plotting the spectrum
ax[1].set_xlabel(r'$f \ [Hz]$')
ax[1].set_ylabel('Power \n Spectrum')
ax[1].set_xticks(np.arange(0,5,1))
ax[1].set_yticks([])
ax[1].text(0.85, 0.004,'(b)',size=11,bbox=dict(facecolor='w',edgecolor='none', alpha=0.5,pad=0.5))    #axs[k].text(3.5, -9, letras[i],size=10,bbox=dict(facecolor='white',edgecolor='none', alpha=1.,pad=0.5))
ax[1].set_xlim(0.0,3.0)
SPF()
F = plt.gcf()
F.set_size_inches(6.26,1.26)
plt.savefig('output/fft.png',dpi=300)
plt.show()


a = 2.0 / N * np.abs(yf[:N//2])[1:]

for i in range (0,len(xf[1:]),1):
    if a[i] == np.amax(a):
        print ('St=',xf[i+1])
        print ('f=',xf[i+1]/tc,'Hz')

import csv
with open("data/probe.dat", 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(t,y))

