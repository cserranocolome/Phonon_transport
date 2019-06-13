#Diffusive regime in a two-dimensional domain


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
L = 10e-7
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
B =2.0e-23    #scattering constant
#Modulo of the phonons' velocity
vg = 8000  
#Debye's frequency
wD = 1e14 
#Time step
dt = L/vg
#Constant by which we devide the number of phonons
C = 1e7



Ntemp = np.load('Ntemp350_1300.npy')
Temp = np.load('Temp350_1300.npy')

#Number of phonons in each cell 
def functionN(w,T):
    return Vcell*1/(np.exp(hbar*w/(kb*T))-1)*w**2/(2*(np.pi**2)*vg**3)/C

N1 = int(round(integrate.quad(functionN, 0, wD, args=T1)[0]))
N2 = int(round(integrate.quad(functionN, 0, wD, args=T2)[0]))
N3 = int(round(integrate.quad(functionN, 0, wD, args=T3)[0]))
N = N1*n1 + N2*n2x*n2y + N3*n3

#Energy of each phonon
def functionE(w,T,Ncell):
    return Vcell/Ncell*hbar*w/(np.exp(hbar*w/(kb*T))-1)*w**2/(2*(np.pi**2)*vg**3)/C

Energies = np.array([])
E1 = integrate.quad(functionE, 0, wD, args=(T1,N1))[0]

#Energies = np.full(N1, E1)
Energies1 = [E1]*N1

#Energies2 = np.array([])
E2 = integrate.quad(functionE, 0, wD, args=(T2,N2))[0]
#Energies2 = np.full(N2*n2x*n2y,E2)
Energies2 = [E2]*N2*n2x*n2y

#Energies3 = np.array([])
E3 = integrate.quad(functionE, 0, wD, args=(T3,N3))[0]
Energies3 = [E3]*N3*n3

Energies = np.append(Energies1, Energies2)
Energies = np.append(Energies, Energies3)
#w
wvector = Energies/hbar


#Initialization

#Number of phonons in each cell
Ns1 = N1
Ns2 = np.zeros((n2x, n2y))
for i in range(len(Ns2)):
    for j in range(len(Ns2[1])):
        Ns2[i,j] = N2
Ns3 = np.array([])

for i in range(n3):
    Ns3 = np.append(Ns3, N3)
    
#Coordinates of each phonon
xyz = [[0,0,0]]*int(N)
for i in range(int(Ns1)):
    xyz[i] = [ np.random.uniform(0,1)*L, L*(1+n2y+np.random.uniform(0,1)), np.random.uniform(0,1)*L]
for i in range(n2x):
    for j in range(n2y):
        for k in range(int(N2)):
            xyz[int(N1)+i*n2y*int(N2)+j*int(N2)+k] = [L*(i+np.random.uniform(0,1)) , L*(1+j+np.random.uniform(0,1)), np.random.uniform(0,1)*L]
for i in range(n3):
    for j in range(int(N3)):
        xyz[int(N1)+int(N2)*n2x*n2y+i*int(N3)+j] = [L*(i+np.random.uniform(0,1)) , L*np.random.uniform(0,1), np.random.uniform(0,1)*L]
xyz = np.array(xyz)

#Direction of velocity
v = [[0,0,0]]*int(N)
for i in range(N):
    phi = 2*np.pi*np.random.uniform(0,1)
    theta = np.arccos(2*np.random.uniform(0,1)-1)
    v[i] = [vg*np.cos(phi)*np.sin(theta), vg*np.sin(phi)*np.sin(theta), vg*np.cos(theta)]
v =np.array(v)


#SIMULATION
k = 0 #Counts the number of iterations
times = [0]*int(N)
times = np.array(times)
Et = np.array([]) #Energy vs t
Tp = np.array([np.zeros((n2x, n2y))])

k_max=5001
keq=1999

while k<k_max:
    xyz = xyz + v*dt
    xdelete = np.array([]) #index of the phonons we will delete in the y boundaries
        #boundaries
    for i in range(N):
        f=0
        while(f == 0):
            f = 1
            if(xyz[i][0]>n2x*L):  #Right x-boundary
                v[i][0] = -v[i][0]
                xyz[i][0] = 2*n2x*L-xyz[i][0] 
                f = 0
            if(xyz[i][0]<0):    #left x-boundary
                v[i][0] = -v[i][0]
                xyz[i][0] = -xyz[i][0]
                f = 0
            if(xyz[i][2]>L ):   #top z-boundary
                v[i][2] = -v[i][2]
                xyz[i][2] = 2*L-xyz[i][2] 
                f = 0
            if(xyz[i][2]<0):   # bottom z-boundary
                v[i][2] = -v[i][2]
                xyz[i][2] = -xyz[i][2]
                f = 0
            if(xyz[i][0]>L and xyz[i][1]>(1+n2y)*L):  #y-boundary
                f = 0
                if(xyz[i][1]-v[i][1]*dt>(1+n2y)*L):
                    v[i][0] = -v[i][0]
                    xyz[i][0] = 2*L-xyz[i][0]
                else:
                    v[i][1] = -v[i][1]
                    xyz[i][1] = 2*(1+n2y)*L-xyz[i][1]
        if(xyz[i][1]>(2+n2y)*L or xyz[i][1]<0):   #delete
            xdelete = np.append(xdelete,i)
    xyz = np.delete(xyz,xdelete,0)
    v = np.delete(v,xdelete,0)
    times = np.delete(times,xdelete)
    wvector = np.delete(wvector, xdelete)
    Energies = np.delete(Energies, xdelete)
    N = N - len(xdelete)
    
        #scattering
    Ps = 1-np.exp(-(A*wvector**4+B*wvector**2*Tmean**3)*times)  

    for i in range(N):
        if(np.random.uniform(0,1)<Ps[i]):
            phi = 2*np.pi*np.random.uniform(0,1)
            theta = np.arccos(2*np.random.uniform(0,1)-1)
            v[i] = [vg*np.cos(phi)*np.sin(theta), vg*np.sin(phi)*np.sin(theta), vg*np.cos(theta)]
            times[i] = 0
        #Energy conservation
    
    index = np.array([])
    xdelete = np.array([])
    N1k = 0        #Hot cell
    for i in range(N):
        if(xyz[i][0]<L and xyz[i][1]>(1+n2y)*L):
            index = np.append(index,i)
            N1k = N1k+1

    if(N1k<N1):
        for i in range(N1-N1k):
            xyz = np.append(xyz,[[ np.random.uniform(0,1)*L, L*(1+n2y+np.random.uniform(0,1)), np.random.uniform(0,1)*L]],0)
            phi = 2*np.pi*np.random.uniform(0,1)
            theta = np.arccos(2*np.random.uniform(0,1)-1)
            v = np.append(v,[[vg*np.cos(phi)*np.sin(theta), vg*np.sin(phi)*np.sin(theta), vg*np.cos(theta)]],0)
            times = np.append(times, 0)
            N = N+1
    if(N1k>N1):
        for i in range(N1k-N1): 
            xdelete = np.append(xdelete,index[i])
        xyz = np.delete(xyz,xdelete,0)
        v = np.delete(v,xdelete,0)
        times = np.delete(times,xdelete)
        N = N-len(xdelete)

    for j in range(n3):
        N3k = 0        #Cold cell
        index = np.array([])
        xdelete = np.array([])
        for i in range(N):
            if(xyz[i][1]<L and L*j<xyz[i][0]<(L*(j+1))):
                index = np.append(index,i)
                N3k = N3k+1

        if(N3k<N3):
            for i in range(N3-N3k):
            	xyz = np.append(xyz,[[L*(j+np.random.uniform(0,1)) , L*np.random.uniform(0,1), np.random.uniform(0,1)*L]],0)
                phi = 2*np.pi*np.random.uniform(0,1)
                theta = np.arccos(2*np.random.uniform(0,1)-1)
                v = np.append(v,[[vg*np.cos(phi)*np.sin(theta), vg*np.sin(phi)*np.sin(theta), vg*np.cos(theta)]],0)
                times = np.append(times, 0)
                N = N+1
        if(N3k>N3):
            for i in range(N3k-N3): 
                xdelete = np.append(xdelete,index[i])
            xyz = np.delete(xyz,xdelete,0)
            v = np.delete(v,xdelete,0)
            times = np.delete(times,xdelete)
            N = N-len(xdelete)
        #Temperature of each subcell
    Ts2 = np.zeros((n2x, n2y))
    Ns2 = np.zeros((n2x, n2y))
    for i in range(n2x):
        for j in range(n2y):
            for l in range(N):
                if(L*j<xyz[l][0]<L*(j+1) and ((n2y-i)*L)<xyz[l][1]<((n2y-i+1)*L)):
                    Ns2[i][j] = Ns2[i][j] + 1

            m = abs(Ntemp-Ns2[i][j])
            Ts2[i][j] = Temp[np.where(m == m.min())]
            

        #Energies
    Energies = np.array([])
    E1 = integrate.quad(functionE, 0, wD, args=(T1,N1))[0]
    Energies1 = [E1]*N1
    
    Energies2 = [0]*int(sum(sum(Ns2)))
    s=0
    for i in range(n2x):
        for j in range(n2y):
            E = integrate.quad(functionE, 0, wD, args=(Ts2[i][j],Ns2[i][j]))[0]
            for l in range(int(Ns2[i][j])):
                Energies2[s] = E
                s = s+1

            
        
    E3 = integrate.quad(functionE, 0, wD, args=(T3,N3))[0]
    Energies3 = [E3]*N3*n3
        
    Energies = np.append(Energies1, Energies2)
    Energies = np.append(Energies, Energies3)
    Et = np.append(Et,sum(Energies))
    wvector = Energies/hbar
    
    if(k>keq):
        Tp = np.append(Tp,[Ts2],0)
        
    times = times + dt
    k = k+1
    
Tp = np.delete(Tp,0,0)
np.save('Tp_dif',Tp)


vf = open("Et_2000.txt","w")
for i in range(len(Et)):
  vf.write("%d "%i)
  vf.write("%e\n"%Et[i])
vf.close()