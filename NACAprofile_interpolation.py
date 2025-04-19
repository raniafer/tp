#script methode des panneaux profil portant méthode HESS AND SMITH
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as lin
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

class panel():
    def __init__(self,N,alpha):
        #Nombre panneaux N
        self.N    =N
        self.X    =np.zeros(N+1)
        self.Y    =np.zeros(N+1)
        self.XC   =np.zeros(N)
        self.YC   =np.zeros(N)
        self.theta=np.zeros(N+1)
        self.BETA =np.zeros((N,N))
        self.PHI  =np.zeros(N)
        self.R    =np.zeros((N,N+1))
        self.S    =np.zeros(N)
        self.SOL  =np.zeros(N+1)
        self.CP   =np.zeros(N)
        self.V    =np.zeros(N)
        #systeme
        self.A=np.zeros((N+1,N+1))
        self.B=np.zeros(N+1)
        #Angle incidence
        self.alpha=alpha
        return
    def geom_NACA(self):
        N=self.N
        dtheta=2.*np.pi/N
        # nombre de points
        Nx=int(N/2) +1
        N_old=98  #98 #68
        X_old=np.zeros(N_old)
        Y_old=np.zeros(N_old)

        # coordonnées intrados et extrados
        X_int=np.zeros(Nx-1)
        X_ext=np.zeros(Nx-2)

        X=np.zeros(N+1)
        Y=np.zeros(N+1)

        # lecture des coordonnées du profil
        X_old, Y_old= read(X_old, Y_old)

        # points sur l'intrados
        for i in range(Nx-1) :
            theta=np.pi-(i+1)/(Nx-1.)*np.pi
            X_int[i]=1./2*(1.-np.cos(theta))
        # points sur l'extrados
        for i in range(Nx-2):
            theta=(i+1)/(Nx-1.)*np.pi
            X_ext[i]=1./2*(1.-np.cos(theta))

        # remplir X + conditions aux limites
        X=np.concatenate(([1.],X_int, X_ext,[1.]))

        # interpolation
        for j in range(0,int(N_old/2)):
            for i in range(0,int(N/2)+1):
                # condition pour trouver l'intervalle ou appartient le point a interpoler
                if ( X_old[j+1]<=X[i]<=X_old[j] ) :
                    pente=(Y_old[j+1]-Y_old[j])/(X_old[j+1]-X_old[j])
                    Y[i]=pente*X[i]+(Y_old[j+1]+Y_old[j])/2-pente*(X_old[j+1]+X_old[j])/2

        for j in range(int(N_old/2)):
            for i in range(1,int(N/2)+1):
                # condition pour trouver l'intervalle ou appartient le point a interpoler
                if ( X_old[int(N_old/2)+j]<=X[i+int(N/2)-1]<=X_old[int(N_old/2)+j+1] ) :
                    pente=(Y_old[int(N_old/2)+j+1]-Y_old[int(N_old/2)+j])/(X_old[int(N_old/2)+j+1]-X_old[int(N_old/2)+j])
                    Y[i+int(N/2)-1]=pente*X[i+int(N/2)-1]+(Y_old[int(N_old/2)+j+1]+Y_old[int(N_old/2)+j])/2-pente*(X_old[int(N_old/2)+j+1]+X_old[int(N_old/2)+j])/2  

        self.X=X
        self.Y=Y

        # conditions aux limites sur Y
        self.Y[0]=0.
        self.Y[N]=0.

        #Calcul des points centraux
        for i in range(N):
            self.XC[i]=(self.X[i]+self.X[i+1])/2.
            self.YC[i]=(self.Y[i]+self.Y[i+1])/2.
           
        #Calcul des angles d orientation des panneaux
        for i in range(N):
            b=self.Y[i+1]-self.Y[i]
            c=self.X[i+1]-self.X[i]
            self.PHI[i]=np.arctan2(b,c)
            self.S[i]=np.sqrt((self.X[i+1]-self.X[i])**2+\
                              (self.Y[i+1]-self.Y[i])**2)
        plt.figure(figsize=(10,2))
        plt.plot(self.X,self.Y, color='k')  #,marker='o',   markersize=8,color='k', linestyle='None')
        plt.plot(self.XC,self.YC,marker='x',   markersize=8,color='red', linestyle='None')
        plt.show()
        for i in range(N):
            for j in range(N):
                if (i==j):
                        self.BETA[i,j]=np.pi
                else:
                        a=(self.YC[i]-self.Y[j+1])*(self.XC[i]-self.X[j])-\
                          (self.XC[i]-self.X[j+1])*(self.YC[i]-self.Y[j])
                        b=(self.XC[i]-self.X[j+1])*(self.XC[i]-self.X[j])+\
                          (self.YC[i]-self.Y[j+1])*(self.YC[i]-self.Y[j])
                        self.BETA[i,j]=np.arctan2(a,b)
        #calcul des distances
        for i in range(N):
           for j in range(N+1):
               self.R[i,j]=np.sqrt((self.X[j]-self.XC[i])**2+\
                                   (self.Y[j]-self.YC[i])**2)

    #Nous remplissons le systeme lineaire pour le profil portant
    def system(self):
        for i in range(N):
            for j in range(N):
                self.A[i,j]=1./(2.*np.pi)*\
                (np.sin(self.PHI[i]-self.PHI[j])*np.log(self.R[i,j+1]/self.R[i,j])+\
                np.cos(self.PHI[i]-self.PHI[j])*self.BETA[i,j])
        for i in range(N):
            for j in range(N):
                self.A[i,N]+=1./(np.pi*2.)*\
                (np.cos(self.PHI[i]-self.PHI[j])*np.log(self.R[i,j+1]/self.R[i,j])-\
                 np.sin(self.PHI[i]-self.PHI[j])*self.BETA[i,j])
            self.B[i]=np.sin(self.PHI[i]-self.alpha)

        #Condition de Kutta
        for j in range(N):
                self.A[N,j]= np.sin(self.PHI[0]-self.PHI[j])*self.BETA[0,j]+\
                             np.sin(self.PHI[N-1]-self.PHI[j])*self.BETA[N-1,j]-\
                             np.cos(self.PHI[0]-self.PHI[j])*np.log(self.R[0,j+1]/self.R[0,j])-\
                             np.cos(self.PHI[N-1]-self.PHI[j])*np.log(self.R[N-1,j+1]/self.R[N-1,j])
                self.A[N,N]+=np.sin(self.PHI[0]-self.PHI[j])*np.log(self.R[0,j+1]/self.R[0,j])+\
                             np.sin(self.PHI[N-1]-self.PHI[j])*np.log(self.R[N-1,j+1]/self.R[N-1,j])+\
                             np.cos(self.PHI[0]-self.PHI[j])*self.BETA[0,j]+\
                             np.cos(self.PHI[N-1]-self.PHI[j])*self.BETA[N-1,j]
        self.A[N,:]=self.A[N,:]/(2.*np.pi)
        self.B[N]=-np.cos(self.PHI[0]-self.alpha)-np.cos(self.PHI[N-1]-self.alpha)
    def calcul(self):
        AI=lin.inv(self.A)
        self.SOL=np.dot(AI,self.B)
    #Calcul de la vitesse tangentielle pour le CP
    def calcul_CP(self):
        for i in range(N):
                self.V[i]=np.cos(self.PHI[i]-self.alpha)
                for j in range(N):
                    self.V[i]+=\
                    self.SOL[j]/(2.*np.pi)*\
                    (np.sin(self.PHI[i]-self.PHI[j])*self.BETA[i,j]-\
                    (np.cos(self.PHI[i]-self.PHI[j])*np.log(self.R[i,j+1]/self.R[i,j])))+\
                    self.SOL[N]/(2.*np.pi)*\
                    (np.sin(self.PHI[i]-self.PHI[j])*np.log(self.R[i,j+1]/self.R[i,j])+\
                    (np.cos(self.PHI[i]-self.PHI[j])*self.BETA[i,j]))
                self.CP[i]=1-self.V[i]**2


# lit les coordonnées du fichier
def read(X, Y):
    with open("data2215.dat", 'r') as f : #data4412.dat
        data=np.loadtxt(f)
    X=data[:,0]
    Y=data[:,1]
    f.close()
    return X, Y


#Nombre de panneaux
N=128
alpha=5.*np.pi/180.
# appeler la classe panel avec les arguments N et l'angle alpha
a=panel(N,alpha)
# appeler les methodes
a.geom_NACA()
a.system()
a.calcul()
a.calcul_CP()

val=0.

for i in range(N):
    val+=a.SOL[N]*a.S[i]
print('coefficient de portance = ',2*val)


N=256
alpha=5.*np.pi/180.
      
b=panel(N,alpha)
b.geom_NACA()
b.system()
b.calcul()
b.calcul_CP()

val=0.

plt.figure()

plt.plot(a.XC,-a.CP,b.XC,-b.CP)
plt.xlabel("x")
plt.ylabel("C_p")
plt.show()

