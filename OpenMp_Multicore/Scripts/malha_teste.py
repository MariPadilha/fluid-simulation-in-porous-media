import numpy as np
import matplotlib.pyplot as plt

# Define os limites e o número de pontos na malha
x_start, x_end = 0, 10
y_start, y_end = 0, 10
n_points = 100

# Cria um espaço linear para x e y
Nx, Ny, y0 = np.loadtxt('data/mesh.dat', skiprows=0, unpack=True)
x = np.loadtxt('data/x.dat', skiprows=0, unpack=True)
y = np.loadtxt('data/y.dat', skiprows=0, unpack=True)

Nx = int(Nx)
Ny = int(Ny)
Nx = int(len(x))
Ny = int(len(y))
y0 = int(y0)
x = np.array(x)
y = np.array(y)	
# Cria uma malha 2D
X, Y = np.meshgrid(x, y)

# Plotando a malha
plt.figure(figsize=(8, 8))
plt.scatter(X, Y, s=1,linestyle='-')  # s é o tamanho dos pontos
plt.title('Malha Computacional 2D')
plt.xlabel('Eixo X')
plt.ylabel('Eixo Y')
plt.grid(True)
plt.savefig('output/malha.png', bbox_inches='tight')

