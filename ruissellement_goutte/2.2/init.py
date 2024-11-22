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
    h[0]   = 0
    h[N+1] = 0
    #
    return

#######################################
# Condition initiale : un seul creneau
#######################################

def init_variable_single(h, x, h0, L, N, eps, theta):
	
	for i in range(N+2) :
                if (0.25*L <= x[i] <= 0.75*L) :
                        h[i]=h0*(1-((x[i]-0.5*L)/(0.25*L))**2)
                else :
                        h[i]=0

	'''
	for i in range(N) :
		if (0.25*L <= x[i] <= 0.5*L) :
			h[i]=4*h0/L* x[i] -h0

		elif (0.5*L <= x[i] <= 0.75*L) :
			h[i]=-4*h0/L* x[i] +3*h0

		else :
			h[i]=0
        '''
    #
    #--------------------
    # Conditions limites
    # Ghost cells
    #--------------------
	h[0]   = 0
	h[N+1]   = 0
    #
	return
