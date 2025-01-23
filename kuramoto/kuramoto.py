# equation Kuramoto Sivashinsky
# Rania Ferhat

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.pyplot import cm

# uk1 u en n+1
# uk u en n
# vk v en n
# vkm v en n-1

def desaliasage(uk,k):
    uk[abs(k)>kmax] = 0;
    return

def stabilite_FG(k, V, dt):
        Lh=-1j*np.diag(k)*V-np.diag(k**2)
        lamb=np.linalg.eigvals(Lh)
        X=dt*lamb.real
        Y=dt*lamb.imag
        plt.plot(X, Y)
        plt.show()

#paramteres numeriques et physiques
nx=100
nt=10000
V=1.

# pas d'espace et de temps
x=np.linspace(-2*np.pi, 2*np.pi, nx)
dx=x[1]-x[0]
dt=0.01

# nombre d'onde
k=np.fft.fftfreq(nx, dx)
kmax=2./3.*max(abs(k))

# zone de stabilite
stabilite_FG(k, V, dt)

# condition initiale et sa transformee de fourier
u0=np.exp(-x**2/0.5**2)
uk=np.fft.fft(u0)
desaliasage(uk, k)

# initialisation
vk=np.fft.fft(u0**2)
vkm=np.zeros(nx)
uk1=np.zeros(nx)

#print ('uk initiale', uk)
#print ('vk initiale', vk)
u_values=np.zeros((nt, nx))

# boucle en temps
for i in range(nt) :
	# pour la premiere iteration
	if (i==0) :
		uk=u0
		vkm=u0**2
	# sinon, mettre a jour u et v
	else :
		uk=uk1
		vkm=vk

	u=np.fft.ifft(uk)
	v=u**2
	vk=np.fft.fft(v)
	desaliasage(vk, k)
	# AB2CN
	uk1=(uk*(1.+dt/2*(k**2-k**4))-vk*3/4*1j*k*dt+vkm/4*1j*k*dt)/(1.-k**2/2*dt+k**4/2*dt)
	#print ('uk', uk)
	#print ('Solving for uk1', uk1)
	#print ('Solving for vk', vk)
	#print ('vkm', vkm)

	# stock des valeur dans une matrice
	u_values[i, :] = np.fft.ifft(uk).real

u=np.real(np.fft.ifft(uk))


plt.plot(x, u0, label='initial')
plt.plot(x, u, label='result')
plt.legend()
plt.show()
