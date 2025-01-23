import numpy as np
import matplotlib.pyplot as plt
import init
import output
import timeloop

#########################
# Principales constantes
#########################
# Numeriques

_N      = 500             # Nombre de cellules par direction
_L      = 1.              # Longueur du domaine
_CFL    = 0.1              # Nombre CFL
_dtmax  = 0.01           # Pas de temps maximum autorisé
_itmax  = 5000         # Nombre d'iterations temporelles max
_tmax   = 0.3            # Temps de simulation max
_tplot  = 0.02             # Temps de sauvegarde
#
# Parametres physiques
_theta= np.pi/2
_g   = 9.81             # Gravite
_rol = 997.0            # Densite de l'eau liquide
_mul = 1.0e-3           # Viscosite de l'eau liquide
_h0  = 0.001         # Hauteur initiale du film
_eps = 1.0e-8           # Hauteur residuelle
#
######################
# Programme principal
######################
# Initialisation des variables
#
x    = np.linspace(0., _L, _N+2)
dx   = x[1]-x[0]
h    = np.zeros(_N+2)                       # Hauteur du film
#--------------------
# Condition initiale
#--------------------
init.init_variable(h, x, _h0, _L, _N, _eps, _theta)
#
#----------------------------------
# Ecriture du pas de temps initial
#----------------------------------
output.write_output(h, x, _N, _h0, '0')
#
#-----------------
# Boucle en temps
#-----------------
timeloop.compute_timeloop(h, x, dx, _CFL, _N, _rol, _mul, _g, _itmax, _tmax, _dtmax, _tplot, _h0, _theta, _L, _eps)
#
#---------------------
# Affichage graphique
#---------------------
#output.write_output(h, x, _N, _h0, 'final')
#
#------------------
# Output to a file
#------------------
#with open("500.txt", 'w') as f:
    #print(h, file=f)
#np.savetxt("500", h, "%.13f")
plt.legend()
plt.show()
