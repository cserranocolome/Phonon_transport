#Ballistic regime in a one-dimensional domain

import sys
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
L = 1.0e-9
n1 = 1
n2x = 1
n2y = 5
n3 = 1
#Volum cell
Vcell = L*L*L
#More constants
hbar = 1.05e-34
kb = 1.38e-23
A = 1.5e-44   #scattering constant 
B =2.0e-24    #scattering constant
#Modulo of the phonons' velocity
vg = 2000.0 
#Debye's frequency
wD = 1e14 
wmin = 0
#Time step
dt = L/vg
#Constant by which we devide the number of phonons
C = 1


Ntemp = np.load('Ntemp250_700.npy')
Temp = np.load('Temp250_700.npy')

#Ntemp = np.load('Ntemp30000_60000.npy')
#Temp = np.load('Temp30000_60000.npy')

#Ntemp = np.load('Ntemp340000_550000.npy')
#Temp = np.load('Temp340000_550000.npy')

#Number of phonons in each cell 
def functionN(w,T):
    return Vcell*1.0/(np.exp(hbar*w/(kb*T))-1.0)*w**2.0/(2.0*(np.pi**2.0)*vg**3.0)/C

N1 = int(round(integrate.quad(functionN, 0, wD, args=T1)[0]))
N2 = int(round(integrate.quad(functionN, 0, wD, args=T2)[0]))
N3 = int(round(integrate.quad(functionN, 0, wD, args=T3)[0]))
N = N1*n1 + N2*n2x*n2y + N3*n3


#Energy of each phonon
def functionE(w,T,Ncell):
    return Vcell/Ncell*hbar*w/(np.exp(hbar*w/(kb*T))-1.0)*w**2.0/(2.0*(np.pi**2.0)*vg**3.0)/C

Energies = np.array([])
E1 = integrate.quad(functionE, wmin, wD, args=(T1,N1))[0]

#Energies = np.full(N1, E1)
Energies1 = [E1]*N1

#Energies2 = np.array([])
E2 = integrate.quad(functionE, wmin, wD, args=(T2,N2))[0]
#Energies2 = np.full(N2*n2x*n2y,E2)
Energies2 = [E2]*N2*n2x*n2y

#Energies3 = np.array([])
E3 = integrate.quad(functionE, wmin, wD, args=(T3,N3))[0]
Energies3 = [E3]*N3*n3

Energies = np.append(Energies1, Energies2)
Energies = np.append(Energies, Energies3)
#w
wvector = Energies/hbar
 

#Initialization

#Number of phonons in each cell
Ns1 = N1
Ns2 = np.array([])
for i in range(n2y):
    Ns2 = np.append(Ns2,N2)
    
Ns3 = np.array([])
for i in range(n3):
    Ns3 = np.append(Ns3, N3)
    
#Coordinates of each phonon and velocity
p=0
xyz = [0.0]*int(N)
xyz = np.array(xyz)
v = [0.0]*int(N)
v = np.array(v)
for i in range(int(Ns1)):
    xyz[i] = L*(1.0+n2y+np.random.uniform(0,1))
    theta = np.arccos(2*np.random.uniform(0,1)-1)
    phi = 2.0*np.pi*np.random.uniform(0,1)
    v[p] = vg*np.sin(theta)*np.cos(phi)
    p=p+1


for j in range(n2y):
    for k in range(int(N2)):
        xyz[p] = L*(1.0+j+np.random.uniform(0,1))
        theta = np.arccos(2*np.random.uniform(0,1)-1)
        phi = 2.0*np.pi*np.random.uniform(0,1)
        v[p] = vg*np.sin(theta)*np.cos(phi)
        p=p+1


for j in range(int(N3)):
    xyz[p] = L*np.random.uniform(0,1)
    theta = np.arccos(2*np.random.uniform(0,1)-1)
    phi = 2.0*np.pi*np.random.uniform(0,1)
    v[p] = vg*np.sin(theta)*np.cos(phi)
    p=p+1


  
#SIMULATION
k = 0 #Counts the number of iterations
times = [0]*int(N)
times = np.array(times)
Et = np.array([]) #Energy vs t

Tcell_list= np.array([[300,250,250,250,250,250,250]])

contar=0
contar_flux=0
if(contar==1):
    vf=open("T_t.txt","w")
if(contar_flux==1):
    vf=open("Flux2_t.txt","w")
    

k_max = 500
while k<k_max:
    xyz = xyz + v*dt
    xdelete = np.array([]) #index of the phonons we will delete in the y boundaries
        #boundaries
    
    if(contar_flux==1):
        
        for i in range(n2y+2):
            flux1=0
            for j in range(N):
                if((n2y+2-i)*L<xyz[j]<(n2y+3-i)*L and (n2y+1-i)*L<(xyz[j]-v[j]*dt)<(n2y+2-i)*L):
                    flux1 = flux1 +1
            vf.write("%i "%flux1)
        vf.write("\n")
        
        for i in range(n2y+2):
            flux1=0
            for j in range(N):
                if((n2y+1-i)*L<xyz[j]<(n2y+2-i)*L and (n2y+1-i)*L<xyz[j]-v[j]*dt<(n2y+2-i)*L):
                    flux1 = flux1 +1
            vf.write("%i "%flux1)
        vf.write("\n")
        
        for i in range(n2y+2):
            flux1=0
            for j in range(N):
                if((n2y-i)*L<xyz[j]<(n2y+1-i)*L and (n2y+1-i)*L<xyz[j]-v[j]*dt<(n2y+2-i)*L):
                    flux1 = flux1 +1
            vf.write("%i "%flux1)
        vf.write("\n")
        vf.write("\n")
    
    
    
    index1 = np.array([])
    index3 = np.array([])
    
    Nvec = [0.0]*(n2y+2)
    Nvec = np.array(Nvec)
    cel = xyz//L
    
    for i in range(N):
        
        if(cel[i]>(n2y+1) or cel[i]<0):
            xdelete = np.append(xdelete,i)
        else:
            Nvec[cel[i]] = Nvec[cel[i]]+1 
            
        
       
    N1k = int(Nvec[n2y+1])
    N3k = int(Nvec[0]) 
    xyz = np.delete(xyz,xdelete)
    v = np.delete(v,xdelete)
    times = np.delete(times,xdelete)
    wvector = np.delete(wvector, xdelete)
    Energies = np.delete(Energies, xdelete)
    N = N - len(xdelete)
    
    cel = xyz//L
    for i in range(N):
        if(cel[i]==(n2y+1)):
            index1 = np.append(index1,i)
        if(cel[i]==0):
            index3 = np.append(index3,i)
            
        #Energy conservation

    if(N1k<N1):
        xyz1 = [0.0]*(N1-N1k)
        v1 = [0.0]*(N1-N1k)
        times1 = [0.0]*(N1-N1k)
        for i in range(N1-N1k):
            xyz1[i]=L*(1+n2y+np.random.uniform(0,1))
            theta = np.arccos(2*np.random.uniform(0,1)-1)
            phi = 2.0*np.pi*np.random.uniform(0,1)
            v1[i] = vg*np.sin(theta)*np.cos(phi)
            N = N+1

        xyz = np.append(xyz,xyz1)
        v = np.append(v,v1)
        times = np.append(times,times1)
    xdelete = np.array([])

    
    if(N1k>N1):
        for i in range(N1k-N1): 
            ind = random.randint(0,len(index1)-1)
            xyz = np.delete(xyz,index1[ind])
            v = np.delete(v,index1[ind])
            times = np.delete(times,index1[ind])
            N = N-1
            index1 = np.delete(index1, ind)
    
    for j in range(n3):

        xdelete = np.array([])

        if(N3k<N3):
            xyz1 = [0.0]*(N3-N3k)
            v1 = [0.0]*(N3-N3k)
            times1 = [0.0]*(N3-N3k)
            for i in range(N3-N3k):
                xyz1[i]=L*np.random.uniform(0,1)
                theta = np.arccos(2*np.random.uniform(0,1)-1)
                phi = 2.0*np.pi*np.random.uniform(0,1)
                v1[i] = vg*np.sin(theta)*np.cos(phi)
                N = N+1
            
            xyz = np.append(xyz,xyz1)
            v = np.append(v,v1)
            times = np.append(times,times1)
            
        if(N3k>N3):
            for i in range(N3k-N3): 
                ind = random.randint(0,len(index3)-1)
                xyz = np.delete(xyz,index3[ind])
                v = np.delete(v,index3[ind])
                times = np.delete(times,index3[ind])
                N = N-1
                index3 = np.delete(index3, ind)
            
            
        #Temperature of each subcell
    Ts2 = [0.0]*(n2y+2)
    Ns2 = [0.0]*n2y
    Ts2 = np.array(Ts2)
    Ts2[0] = T1
    Ns2 = np.array(Ns2)
    
    for j in range(n2y):
        Ns2[j] = Nvec[n2y-j]
        m = abs(Ntemp-Ns2[j])
        Ts2[j+1] = Temp[np.where(m == m.min())]
        
    Ts2[n2y+1]=T3
    
    
        #Energies
    Energies = np.array([])
    E1 = integrate.quad(functionE, wmin, wD, args=(T1,N1))[0]
    Energies1 = [E1]*N1
    
    Energies2 = [0.0]*int(sum(Ns2))
    s=0
    for i in range(n2x):
        for j in range(n2y):
            E = integrate.quad(functionE, wmin, wD, args=(Ts2[j],Ns2[j]))[0]
            for l in range(int(Ns2[j])):
                Energies2[s] = E
                s = s+1

            
        
    E3 = integrate.quad(functionE, wmin, wD, args=(T3,N3))[0]
    Energies3 = [E3]*N3*n3
        
    Energies = np.append(Energies1, Energies2)
    Energies = np.append(Energies, Energies3)
    Et = np.append(Et,sum(Energies))
    wvector = Energies/hbar
    
    if(contar ==1):
        N1k=0
        N3k=0
        for i in range(N):
            if(xyz[i]>(n2y+1)*L):
                N1k=N1k+1
            if(xyz[i]<L):
                N3k=N3k+1
        
        m = abs(Ntemp-N1k)
        N1kT = Temp[np.where(m == m.min())]
        
        vf.write("%f "%N1kT)         
        for j in range(n2y):
        	vf.write("%f "%Ts2[j])
        
        m = abs(Ntemp-N3k)
        N3kT = Temp[np.where(m == m.min())]
        
        vf.write("%f "%N3kT) 
        vf.write("\n")
          
    Tcell_list= np.append(Tcell_list, [Ts2], 0)   
    times = times + dt
    k = k+1

if(contar==1):
    vf.close()

if(contar_flux==1):
    vf.close()



vf = open("Def_bal_9.txt", "w")
for i in range(len(Tcell_list[0])):
    for j in range(len(Tcell_list)):
        vf.write("%f "%Tcell_list[j][i])
    vf.write("\n")
vf.close()    

vf = open("EDef_bal_9.txt","w")
for i in range(len(Et)):
  vf.write("%d "%i)
  vf.write("%e\n"%Et[i])
vf.close()





