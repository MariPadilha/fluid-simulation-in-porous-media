import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import numpy as np
from scipy.interpolate import RectBivariateSpline
from strp2 import *

S = 4.65
rad = 1.0

print("COMEÇO GRAPHICS.PY")

Nx ,Ny , y0 = np.loadtxt('data/mesh.dat', skiprows=0, unpack=True)
x = np.loadtxt('data/x.dat', skiprows=0, unpack=True)
y = np.loadtxt('data/y.dat', skiprows=0, unpack=True)

P   = np.loadtxt('data/P.dat', skiprows=0, unpack=True)
T   = np.loadtxt('data/T.dat', skiprows=0, unpack=True)


speedmin = 1000
speedmax = 0
for i in range (100,110,1):
	u = np.loadtxt('transient/data/u'+str(i)+'.dat', skiprows=0, unpack=True)
	v = np.loadtxt('transient/data/v'+str(i)+'.dat', skiprows=0, unpack=True)
	u = np.array(u)
	v = np.array(v)
	speed = np.sqrt(u*u + v*v)
	for speedaux in speed:
		auxmin = min(speedaux)
		auxmax = max(speedaux)
		if auxmin < speedmin:
			speedmin = auxmin
		if auxmax > speedmax:
			speedmax = auxmax
speedmin = round(speedmin-0.5)
speedmax = round(speedmax+0.5)
for i in range (100,110,1):
	print(i)
	u = np.loadtxt('transient/data/u'+str(i)+'.dat', skiprows=0, unpack=True)
	v = np.loadtxt('transient/data/v'+str(i)+'.dat', skiprows=0, unpack=True)

	x_grid,y_grid= np.loadtxt('data/grid.dat', skiprows=0, unpack=True)
	x_grid2,y_grid2= np.loadtxt('data/grid_boundary.dat', skiprows=0, unpack=True)
	x_grid3,y_grid3= np.loadtxt('data/grid_boundary_side.dat', skiprows=0, unpack=True)
	x_grid4,y_grid4= np.loadtxt('data/grid_droplet.dat', skiprows=0, unpack=True)

	Nx = int(Nx)
	Ny = int(Ny)
	y0 = int(y0)

	x = np.array(x)
	y = np.array(y)
	u = np.array(u)
	v = np.array(v)
	
	x_int = x
	y_int = y
	Nx = int(len(x))
	Ny = int(len(y))

	x_int = np.linspace(min(x), max(x), num=Nx,endpoint='True')
	y_int = np.linspace(min(y), max(y), num=Ny,endpoint='True')

	X,Y = np.meshgrid(x,y)
	X2,Y2 = np.meshgrid(x_int,y_int)
	interp_spline = RectBivariateSpline(y,x,v)
	v_int = interp_spline(y_int,x_int)

	X,Y = np.meshgrid(x,y)
	X2,Y2 = np.meshgrid(x_int,y_int)
	interp_spline = RectBivariateSpline(y,x,u)
	u_int = interp_spline(y_int,x_int)

	speed = np.sqrt(u*u + v*v)


############# GRID
	plt.figure(figsize=(8, 8))
	plt.plot( x_grid  ,y_grid  ,'k-', linewidth=1.0, label='Grid')
	plt.plot( x_grid2  ,y_grid2  ,'bo', linewidth=0.50, label='Boundary')
	plt.plot( x_grid3  ,y_grid3  ,'ro', linewidth=0.50, label='Boundary_side')
	plt.plot( x_grid4  ,y_grid4  ,'go', linewidth=0.50, label='Boundary_side')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.axis('scaled')
	plt.legend( loc= 'upper right')
	plt.savefig('output/mesh.png', bbox_inches='tight')
	plt.close()
########### PRESSURE
	plt.figure(figsize=(8, 8))
	plt.contourf(x,y,P,50, cmap=plt.cm.get_cmap('RdBu_r'))
	plt.colorbar()
	plt.axis('scaled')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Pressure - P')
	plt.savefig('output/pressure.png', bbox_inches='tight')

########### TEMPERATURE
	plt.figure(figsize=(8, 8))
#	plt.contour(x,y,Z,colors=('k',),linestyles=('--',),linewidths=(2,), levels=[1/(S+1)])
	plt.contourf(x,y,T,50, cmap=plt.cm.get_cmap('RdBu_r'))
	plt.colorbar()
	plt.axis('scaled')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Temperature - T')
	plt.savefig('output/temperature.png', bbox_inches='tight')
#	plt.show()

############# STREAMLINES
	print(speedmin)
	print(speedmax)
	X_stream = np.linspace(0,40,100)
	Y_stream = X_stream*+9 -10 #y[len(y)/2]

	plt.figure(figsize=(8, 8))
	plt.streamplot(x_int, y_int, u_int, v_int,density=3, linewidth=0.8, color='k',arrowstyle='-')
	plt.contourf( x,y,speed,50, cmap=plt.cm.get_cmap('RdBu_r'))
	plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=speedmin, vmax=speedmax), cmap=plt.cm.get_cmap('RdBu_r')))	
	plt.title('Velocity Magnitude')
	plt.axis('scaled')
	plt.ylim(min(y),max(y))
	plt.xlim(min(x),max(x))
	plt.savefig('output/streamlines'+str(i)+'.png', bbox_inches='tight')
