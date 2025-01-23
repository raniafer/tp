import numpy as np
import os
from numpy import loadtxt
import matplotlib.pyplot as plt

file_n=["16", "32","64","128", "256", "512"]


for i in range(len(file_n)):
	files=open(file_n[i], 'r')

h_16 = np.loadtxt(file_n[0])

h_32 = np.loadtxt(file_n[1])

h_64 = np.loadtxt(file_n[2])

h_128 = np.loadtxt(file_n[3])

h_256 = np.loadtxt(file_n[4])

h_512 = np.loadtxt(file_n[0])

files.close()

def erreur(h_max, h, e, n) :
	h_avr= np.add.reduceat(h_max, np.arange(0, int(len(h_max)), n ))
	e = np.sqrt ( 1/n * np.sum ( (h[0:h_avr.size]-h_avr[:])**2  )  )
	return e
e=0

table = np.array ([erreur(h_1000, h_10, e, 100), erreur(h_1000, h_125, e, 8), erreur(h_1000, h_250, e, 4), erreur(h_1000, h_500, e, 2)  ])
print (table)
plt.plot(file_n[:-1], np.log(table))
plt.xlabel('nombre de cellules')
plt.ylabel('ln($\epsilon$)')
plt.show()
