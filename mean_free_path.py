#Mean free path

import numpy as np
import scipy as sp
import scipy.integrate as integrate
#import matplotlib.pyplot as plt
import random
#from tqdm import tqdm
import time
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

#Initial temperatures
T1 = 300
T2 = 250
T3 = 250
Tmean = (T1+T3)/2
#Dimensions of the surface
L = 1e-9
n1 = 1
n2x = 5
n2y = 5
n3 = 5
#Volum cell
Vcell = L*L*L
#More constants
hbar = 1.05e-34
kb = 1.38e-23
A = 1.5e-44   #scattering constant 
B =2.0e-24    #scattering constant
#Modulo of the phonons' velocity
vg = 2000.  
#Debye's frequency
wD = 1e14 
#Time step
dt = L/vg
#Constant by which we devide the number of phonons
C = 1



Ntemp = np.load('Ntemp250_700.npy')
Temp = np.load('Temp250_700.npy')

#Number of phonons in each cell 
def functionN(w,T):
    return Vcell*1/(np.exp(hbar*w/(kb*T))-1)*w**2/(2*(np.pi**2)*vg**3)

N1 = int(round(integrate.quad(functionN, 0, wD, args=T1)[0]/C))
N2 = int(round(integrate.quad(functionN, 0, wD, args=T2)[0]/C))
N3 = int(round(integrate.quad(functionN, 0, wD, args=T3)[0]/C))
N = N1*n1 + N2*n2x*n2y + N3*n3

#Energy of each phonon
def functionE(w,T,Ncell):
    return Vcell/Ncell*hbar*w/(np.exp(hbar*w/(kb*T))-1)*w**2/(2*(np.pi**2)*vg**3)


N_Tm = int(round(integrate.quad(functionN, 0, wD, args=Tmean)[0]/C))
E_Tm = integrate.quad(functionE, 0, wD, args=(Tmean,N_Tm))[0]

w_Tm = E_Tm/hbar

Lambda3 = vg/(A*w_Tm**4+B*w_Tm**2*Tmean**3)
Lambda = vg/(A*wD**4+B*wD**2*Tmean**3)

vf = open("mean_free_path.txt", "w")
vf.write("w : %f\n"%w_Tm)
vf.write("LambdaMean : %.10f\n"%Lambda3)
vf.write("LambdaDebye : %.10f\n"%Lambda)

vf.close()
