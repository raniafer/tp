import numpy as np
import output

######################################
# Calcul de la vitesse aux interfaces
######################################
N=100
uhalf= np.zeros(N+1)
def compute_velocity(uhalf, h, x, N, rol, mul, g, theta, L): # en espace
	for i in range(N):
#	uhalf[0:N+1] = rol*h[0:N+1]**2*g/(3*mul)*( (h[0:N+1] - h[1:N+2]) / (x[0:N+1] - x[1:N+2]) )+np.sin(theta)
		uhalf[i] = rol*h[i]**2*g/(3*mul)*( (h[i+1]-h[i]) / (x[i+1]-x[i]) *np.cos(theta) +np.sin(theta) )
	return


##################
# Boucle en temps
##################
def compute_timeloop(h, x, dx, CFL, N, rol, mul, g, itmax, tmax, dtmax, tplot, h0, theta, L, eps):
    #-----------------------
    # Variables temporaires
    #-----------------------
    flux  = np.zeros(N+2)
#    uhalf = np.zeros(N+2)
    #
    #-------------------------------------
    # initialiation de la boucle en temps
    #-------------------------------------
    continuer = True
    it = 1
    t  = 0.0
    tsave = tplot
    #--------------
    # Calcul de dt
    #--------------
    dt = min(dtmax, CFL*3.0*mul*dx*dx/(2.0*rol*g*(h0**3)))
    print('dt=', dt)
    #
    #-------------------
    # Boucle temporelle
    #-------------------
    while (continuer):
        print('it=', it,' t=',t,'s')
        print(h)
        #
        #--------------------------
        # Calcul du flux numerique
        #--------------------------
        compute_velocity(uhalf, h, x, N, rol, mul, g, theta, L)
        
        h[0]   = eps
        h[N] = eps
        
        flux[0:N]=np.maximum(uhalf[0:N],0)*h[0:N]+np.minimum(uhalf[0:N],0)*h[1:N+1]  # cest bon!
           
        #
        #----------------------
        # Schema volumes finis
        #----------------------

        h[1:N+1]=h[0:N]-dt/dx*(flux[0:N]-flux[1:N+1])
        #------------------------
        # Conditions aux limites
        # Ghost cells
        #------------------------

        #        
        print('-------------------')
        #
        it = it + 1
        t  = t  + dt
        #
        #-----------------
        # Test d'ecriture
        #-----------------
        if (t > tsave):
            output.write_output(h, x, N, h0, str(t))
            tsave = tsave + tplot
        #
        #----------------------
        # Test de continuation
        #----------------------
        if (it > itmax):
            continuer = False
        else:
            continuer = True
            if (t > tmax):
                it = itmax
                dt = tmax - t + dt
                t  = tmax
    return
