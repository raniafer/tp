import numpy as np
import output

######################################
# Calcul de la vitesse aux interfaces
######################################
N=20
uhalf= np.zeros(N+2)
k = np.zeros(N+2)
'''
def courbure(k, h, dx, N, x) :
        k[1:N]= (h[2:N+1]+h[0:N-1]-2*h[1:N])/dx**2
        k[0]=0
        k[N+1]=0
        return k
'''
def compute_velocity(uhalf, h, x, N, rol, mul, g, theta, L, dx, k, sigma):
        k[1:N]= (h[2:N+1]+h[0:N-1]-2*h[1:N])/dx**2
        k[0]=0
        k[N+1]=0
        uhalf[1:N+1]= - h[1:N+1]**2/(3*mul)*( g*rol*(h[1:N+1]-h[0:N])/dx *np.cos(theta) + g*rol*np.sin(theta)
		      - sigma*(k[1:N+1]-k[0:N])/dx )
        print ('sigma * (k[1:N+1]-k[0:N]) / dx : ', sigma*(k[1:N+1]-k[0:N])/dx )
        return


##################
# Boucle en temps
##################
def compute_timeloop(h, x, dx, CFL, N, rol, mul, g, itmax, tmax, dtmax, tplot, h0, theta, L, eps, sigma):
    #-----------------------
    # Variables temporaires
    #-----------------------
    flux  = np.zeros(N+2)
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
        print('it=', it,' t=',t,'s', 'dt=', dt)
        print(h)
        print("h[:]:", h[:])
        #
        #--------------------------
        # Calcul du flux numerique
        #--------------------------
        compute_velocity(uhalf, h, x, N, rol, mul, g, theta, L, dx, k, sigma)
        print ("uhalf: ",uhalf[:])
        h[0] = h[1]
        h[N+1] = h[N]
        
        flux[1:]=np.maximum(uhalf[1:],0)*h[:-1]+np.minimum(uhalf[1:],0)*h[1:]
        #
        #----------------------
        # Schema volumes finis
        #----------------------
        #h[1:]=h[:-1]- dt/dx * (flux[:-1]-flux[1:])
        #hn=h.copy()
        h[1:]=h[1:]- dt/dx * (flux[:-1]-flux[1:])
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
