import numpy as np
import matplotlib.pyplot as plt

###################
# Sortie graphique
###################
def write_output(h, x, N, h0, t):
    xout    = np.zeros(N)
    xout[:] = 0.5*(x[0:N] + x[1:N+1])
    #
    #fig, ax = plt.subplots()
    #ax.plot(xout, h[1:N+1])
    plt.plot(xout, h[1:N+1], label=f"{float(t):.3f}s")
    #
    plt.ylabel('h(x)')
    plt.xlabel('x')
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    #plt.show()
    return
