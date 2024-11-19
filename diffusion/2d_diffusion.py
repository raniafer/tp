# DrZGan tutorial
import numpy as np
from matplotlib import pyplot as plt, cm
from mpl_toolkits.mplot3d import Axes3D

nx=10; ny=10
nt=50; tmax=100
lx=1 ; ly=2

nu=1.e-3

x=np.linspace(0., lx, nx)
y=np.linspace(0., ly, ny)
dx=x[1]-x[0]
dy=y[1]-y[0]
dt=tmax/nt

# conditions initiales
# u  en n+1
# un en n
u=np.zeros((nx, ny))
un=np.zeros((nx, ny))

#u[int(0.25*lx):int(0.75*lx) , : ]=1
u[int(0.5/dx):int(1/dx+1),int(0.5/dy):int(1/dy+1)] = 2
print ("initialisation", u)

def diffusion(nt, dt, dx, dy, u) :

	for n in range(nt) :
		un=u.copy()
		u[1:-1,1:-1]=un[1:-1,1:-1]+ nu*dt/dx**2*( un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1] ) + nu*dt/dy**2*( un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2] )
	print ("resultat", u)
	return u

uk=diffusion(nt, dt, dx, dy, u)

fig = plt.figure(figsize = (11,7), dpi=100)
ax = fig.add_subplot(111, projection='3d')
# The '111' means a grid of 1 row and 1 column and this subplot is the first one.
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,uk,cmap=cm.viridis)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$');
plt.show()
