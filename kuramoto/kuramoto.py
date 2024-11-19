import numpy as np
import matplotlib.pyplot as plt 
import integrale_temporelle as it
from matplotlib.pyplot import cm

n=64
nt=100
dt=0.001
V=1.
nu=0.01


x=np.linspace(0, 2*np.pi, n)
t = np.linspace(start=0, stop=5, num=nt)
#x = np.linspace(start=0, stop=L, num=n)
dx=x[1]-x[0]
k=np.fft.fftfreq(n)*n

#the stability area
it.stabilite_FG(k, V, nu, dt)

#u0=1+np.exp(-x**2/0.5**2)	# condition initiale
#u0=np.zeros(n)+5
#u0=np.cos(x/16)*(1+np.sin(x/16))    # https://people.maths.ox.ac.uk/trefethen/pdectb/kuramoto2.pdf
#u0 = np.cos((2 * np.pi * x) / L) + 0.1 * np.cos((4 * np.pi * x) / L)
u0=-np.sin(x/10) # https://www.researchgate.net/publication/350160417_Evolutional_Deep_Neural_Network
uk=np.fft.fft(u0)

# uk1 uk en n+1
# uk  uk en n
# uk0 uk en n-1

# vk1 vk en n+1
# vk  vk en n
# vk0 vk en n-1

vk1=np.zeros(n)

def CFL(dx, u):
	dt=0.8*dx/u
	return(dt)
	
u_values = np.zeros((nt, n))	

# Crank Nicholson bleme: comment avoir vk1
for i in range(nt):
	for j in range(n):
		u=np.fft.ifft(uk)
		v=u**2
		vk=np.fft.fft(v)
		uk1=(uk*(1+k**2*dt-k**4*dt)+vk1*(-1j*k*dt/4)+vk*(-1j*k*dt/4))/(1-k**2/2*dt+k**4*dt)
		#dt=CFL(dx, max(np.fft.ifft(uk).real))
	#print ('CFL = ', dt*max(np.fft.ifft(uk).real)/dx)
	#print('dt = ', dt)
	u_values[i, :] = np.fft.ifft(uk).real
	uk0=uk.copy()
	uk=uk1.copy()
	vk=vk1.copy()

u=np.fft.ifft(uk).real


plt.plot(x, u0, label='initial')
plt.plot(x, u, label='result')
for i in range(nt):
	plt.plot(x, u_values[i, :])
plt.show()

# plot the result
fig, ax = plt.subplots(figsize=(10,8))

xx, tt = np.meshgrid(x, t)
levels = np.arange(-3, 3, 0.01)
cs=ax.contourf(xx, tt, u_values, cmap=cm.jet)
fig.colorbar(cs)
'''
plt.plot(x,u0,label='initial');
plt.plot(x,u,label='result');
plt.show()
'''
ax.set_xlabel("x")
ax.set_ylabel("t")
plt.show()
