import numpy as np
import matplotlib.pyplot as plt

######################################
# Definition des conditions initiales
######################################
def init_variable(h, x, h0, L, N, eps, theta):
    init_variable_single(h, x, h0, L, N, eps, theta)
    #init_variable_double(h, x, h0, L, N, eps)
    return

#######################################
# Condition initiale : deux créneaux
#######################################
def init_variable_double(h, x, h0, L, N, eps):
    #
    #--------------------
    # Conditions limites
    # Ghost cells
    #--------------------
    h[0]   = h0
    h[N+1] = h0
    #
    return
''' 
_theta=30
_L=1.
_N=20
_h0=1
eps=1.e-3

#x    = np.linspace(-_L/2*np.cos(_theta), _L/2*np.cos(_theta), _N+1)
x    = np.linspace(0., _L, _N+1)
h    = np.zeros(_N+1)
'''
#######################################
# Condition initiale : un seul creneau
#######################################

def init_variable_single(h, x, h0, L, N, eps, theta):
		
	for i in range(N) :
		if (0.25*L <= x[i] <= 0.5*L) :
			h[i]=4*h0/L* x[i] -h0
					
		elif (0.5*L <= x[i] <= 0.75*L) :
			h[i]=-4*h0/L* x[i] +3*h0
			
		else : 
			h[i]=0			
		
		#h[i]=-4*h0/(L*np.cos(theta))**2*(x[i])**2+h0
	#h[1:N+1]=x[0:N]**2 * ( -4*h0/(L*np.cos(theta))**2 ) + x[0:N] * ( np.sin(theta)/np.cos(theta)+4*(h0+h0*np.sin(theta))/np.cos(theta) )

    #
    #--------------------
    # Conditions limites
    # Ghost cells
    #--------------------
	h[0]   = eps
	h[N]   = eps
    #
	return #h
'''
h=init_variable_single(h, x, _h0, _L, _N, eps, _theta)

plt.plot(x, h)
plt.show()
'''
