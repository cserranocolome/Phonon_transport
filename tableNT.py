#Code to generate the table N-T.

import numpy as np
import scipy.integrate as integrate
from pynverse import inversefunc

#Initial temperatures
T1 = 300
T2 = 250
T3 = 250
Tmean = (T1+T3)/2
#Dimensions of the surface
L = 1e-8
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
vg = 2000  #Think how to find it correctly
#Debye's frequency
wD = 1e14 #Think about it
#Time step
dt = L/vg
#Constant by which we devide the number of phonons
C = 10



#Number of phonons in each cell 
def functionN(w,T):
    return Vcell*1/(np.exp(hbar*w/(kb*T))-1)*w**2/(2*(np.pi**2)*vg**3)/C

N1 = round(integrate.quad(functionN, 0, wD, args=T1)[0]/C)
N2 = round(integrate.quad(functionN, 0, wD, args=T2)[0]/C)
N3 = round(integrate.quad(functionN, 0, wD, args=T3)[0]/C)
N = N1*n1 + N2*n2x*n2y + N3*n3


def functionNT(T):
    return integrate.quad(functionN, 0, wD, args=T)[0]


inv = inversefunc(functionNT, domain=[200,400])

functionNT(250)
functionNT(300)



#Table of the number of phonons and temperature
Ns = np.arange(30000,60000,1)
Ts = [0]*len(Ns)
for i in range(len(Ns)):
    Ts[i] = inv(Ns[i])
Ts = np.array(Ts)

np.save('Ntemp30000_60000',Ns)
np.save('Temp30000_60000',Ts)




