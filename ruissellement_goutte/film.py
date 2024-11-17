import numpy as np
import matplotlib.pyplot as plt
import init
import output
import timeloop

#########################
# Principales constantes
#########################
# Numeriques

_N      = 100             # Nombre de cellules par direction
_L      = 1.              # Longueur du domaine
_CFL    = 0.95              # Nombre CFL
_dtmax  = 0.1           # Pas de temps maximum autorisé
_itmax  = 100         # Nombre d'iterations temporelles max
_tmax   = 0.5            # Temps de simulation max
_tplot  = 0.1             # Temps de sauvegarde
#
# Parametres physiques
_theta=30
_g   = 9.81             # Gravite
_rol = 997.0            # Densite de l'eau liquide
_mul = 1.0e-3           # Viscosite de l'eau liquide
_h0  = 1.0           # Hauteur initiale du film
_eps = 1.0e-8           # Hauteur residuelle
#
######################
# Programme principal
######################
# Initialisation des variables
#
x    = np.linspace(0., _L, _N+1)
dx   = x[1]-x[0]
h    = np.zeros(_N+1)                       # Hauteur du film
#uhalf= np.zeros(_N+2)
#--------------------
# Condition initiale
#--------------------
init.init_variable(h, x, _h0, _L, _N, _eps, _theta)
#
#----------------------------------
# Ecriture du pas de temps initial
#----------------------------------
output.write_output(h, x, _N, _h0, '0') # UNCOMMENT!!!!
#
#-----------------
# Boucle en temps
#-----------------
timeloop.compute_timeloop(h, x, dx, _CFL, _N, _rol, _mul, _g, _itmax, _tmax, _dtmax, _tplot, _h0, _theta, _L, _eps)
#
#---------------------
# Affichage graphique
#---------------------

output.write_output(h, x, _N, _h0, 'final')

#plt.plot(x, h)
plt.legend()
plt.show()
