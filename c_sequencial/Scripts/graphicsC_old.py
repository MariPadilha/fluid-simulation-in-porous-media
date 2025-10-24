import matplotlib.pyplot as plt
import numpy as np
from numpy import arange
from scipy.interpolate import RectBivariateSpline
from strp2 import *

S = 4.65

rad = 1.0




Nx ,Ny , y0 = np.loadtxt('data/mesh.dat', skiprows=0, unpack=True)
x = np.loadtxt('data/x.dat', skiprows=0, unpack=True)
y = np.loadtxt('data/y.dat', skiprows=0, unpack=True)
#	u = np.loadtxt('transient/data/u'+str(i)+'.dat', skiprows=0, unpack=True)
#	v = np.loadtxt('transient/data/v'+str(i)+'.dat', skiprows=0, unpack=True)
	
P   = np.loadtxt('data/P.dat', skiprows=0, unpack=True)
T   = np.loadtxt('data/T.dat', skiprows=0, unpack=True)
C   = np.loadtxt('data/C.dat', skiprows=0, unpack=True)

cont = 100

for i in range (100,1000,1):
#	print (i)
	cont = cont+1
 	#print (i)

	C   = np.loadtxt('transient/data/C'+str(i)+'.dat', skiprows=0, unpack=True)
#itc, error_u, error_v, error_p, mass, residual_p = np.loadtxt('data/error.dat', skiprows=0, unpack=True)

	x_grid,y_grid= np.loadtxt('data/grid.dat', skiprows=0, unpack=True)
#x_grid_u,y_grid_u= np.loadtxt('data/grid_u.dat', skiprows=0, unpack=True)
#x_grid_v,y_grid_v= np.loadtxt('data/grid_v.dat', skiprows=0, unpack=True)
#x_grid2,y_grid2= np.loadtxt('data/grid_boundary.dat', skiprows=0, unpack=True)
#x_grid3,y_grid3= np.loadtxt('data/grid_boundary_side.dat', skiprows=0, unpack=True)
#x_grid4,y_grid4= np.loadtxt('data/grid_droplet.dat', skiprows=0, unpack=True)

	Nx = int(Nx)
	Ny = int(Ny)
	y0 = int(y0)

	x = np.array(x)
	y = np.array(y)
#	u = np.array(u)
#	v = np.array(v)
	C = np.array(C)

#transform to uniform mesh to be able to plot using pythons streamplot function
#	x_int = x
#	y_int = y
#	Nx = int(len(x))
#	Ny = int(len(y))

#	x_int = np.linspace(min(x), max(x), num=Nx,endpoint='True')
#	y_int = np.linspace(min(y), max(y), num=Ny,endpoint='True')

#	X,Y = np.meshgrid(x,y)
#	X2,Y2 = np.meshgrid(x_int,y_int)
#	interp_spline = RectBivariateSpline(y,x,v)
#	v_int = interp_spline(y_int,x_int)

#	X,Y = np.meshgrid(x,y)
#	X2,Y2 = np.meshgrid(x_int,y_int)
#	interp_spline = RectBivariateSpline(y,x,u)
#	u_int = interp_spline(y_int,x_int)

#	speed = np.sqrt(u*u + v*v)

########### TEMPERATURE
	plt.figure(figsize=(8, 8))
	#CS = plt.contourf(x,y,speed,50,cmap='RdBu_r') 
#	CS.set_clim(0,10)
#	plt.contour(x,y,C,colors=('k',),linestyles=('--',),linewidths=(2,), levels=[1/(S+1)])
	plt.contourf(x,y,C,100, cmap=plt.cm.get_cmap('RdBu_r'))
#plt.plot( x_grid  ,y_grid  ,'k-', linewidth=0.2, label='Grid')
	plt.colorbar()
	plt.axis('scaled')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Concentration')
	circle = plt.Circle((0, 0), rad, color='k', linewidth=2.0, fill=False)
#	circle2 = plt.Circle((0, 3), rad, color='k', linewidth=2.0, fill=False)
#	circle3 = plt.Circle((0, 8), rad, color='k', linewidth=2.0, fill=False)
	plt.gca().add_patch(circle)
#plt.gca().add_patch(circle2)
#plt.gca().add_patch(circle3)
	plt.savefig('output/concentration'+str(cont)+'.png', bbox_inches='tight')


