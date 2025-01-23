import numpy as np
import matplotlib.pyplot as plt
'''
n=128 # nombre de points 
nt=75 # pas en temps
nu=0.01
k=-n/2
V=1.

dt=0.01

k=np.fft.fftfreq(n)*n
'''
def stabilite_FG(k, V, dt):
# ajoute le pas de temsp
	Lh=-1j*np.diag(k)*V-np.diag(k**2)
	lamb=np.linalg.eigvals(Lh)
	X=dt*lamb.real
	Y=dt*lamb.imag
	plt.plot(X, Y)
	plt.show()

def AB_CN(k):
	x=np.linspace(-np.pi,np.pi,n)
	#uk=np.random.randn(n)
	u0 = np.exp(-(x/0.5)**2)
	uc = np.fft.fft(u0)
	ucm= uc.copy()
	ucp= 0*uc
	plt.plot(x, u0)
	for i in range(nt):
		for j in range(n):
			ucp[j]=((1.-dt/2*nu*k[j]**2)*uc[j]-dt*1j*k[j]*V/2*(3*uc[j]-ucm[j]))/(1+dt/2*nu*k[j]**2)
		ucm = uc.copy()
		uc  = ucp.copy()
		
	u=np.fft.ifft(ucp).real
	return(u)
'''
x=np.linspace(-np.pi,np.pi,n)
U=AB_CN(k)
plt.plot(x, U)
plt.show()

N=128
Nt=10000
dt=0.00001
'''
def RK_CN():
	x=np.cos(np.pi*np.arange(N)/(N-1))
	D=np.zeros((N,N))
	for i in range(N):
		for j in range(N):
			ci = 1.; cj = 1.;
			if(i==0 or i==N-1): ci = 2.
			if(j==0 or j==N-1): cj = 2.
			if(i!=j):
				D[i,j] = (ci/cj)*(-1)**(i+j)/(x[i]-x[j]);
		for i in range(1,N-1):
			D[i,i] = -x[i]/(2.*(1.-x[i]**2))
		D[0,0] = (2.*(N-1)**2+1.)/6.
		D[-1,-1] = -(2.*(N-1)**2+1.)/6.
	D2=np.dot(D,D)
	A  = np.eye(N) - 0.5*dt*nu*D2	
	A[0,:] = 0;
	A[0,0] = 1.;
	A[N-1,:] = 0;
	A[N-1,N-1] = 1.;
	A2 = np.linalg.inv(A)
	u0 = np.exp(-(x/0.1)**2)#np.cos(0.5*np.pi*x);  # condition initiale
	ui = u0.copy();
	plt.plot(x,u0)
	for n in range(Nt):
		f1     = np.dot(np.eye(N) - dt*V*D + nu*0.5*dt*D2 , ui ); # premier second membre
		f1[0] = 0.; f1[N-1] = 0.; # conditions aux limites homogènes
		utilde = np.dot(A2,f1);                                   # predication RK2
		f2     = np.dot(np.eye(N) + nu*0.5*dt*D2 , ui) - 0.5*dt*V*np.dot(D,ui+utilde); # deuxieme second membre
		f2[0] = 0.; f2[N-1] = 0.; # conditions aux limites homogènes
		ui     = np.dot(A2,f2);     
		                              # correction RK2
	return(ui)
'''
ui=RK_CN()
x1=np.cos(np.pi*np.arange(N)/(N-1))
plt.plot(x1,ui);
plt.show();'''
