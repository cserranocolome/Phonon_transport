#Statistical treatment of the data collected.

import numpy as np


#Ballistic one-dimensional.
Tcell_list = np.loadtxt("Def_dif_8_500.txt", dtype='i')

Tmean_list = np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])
lon = 1000
for i in range(int((len(Tcell_list[0]))/lon)):
    s = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for j in range(lon*i,lon*(i+1)):
        s = s+Tcell_list[:,j]
    Tmean_list = np.append(Tmean_list, [s/lon], 0)
     
Tmean_list = np.delete(Tmean_list,0,0)

vf = open("Tmean_list_barra_balistic.txt","w")
for i in range(7):
    for j in range(len(Tmean_list)):
        vf.write("%f "%Tmean_list[j][i])
    vf.write("\n")
    

Teq = Tcell_list[:,145000:150000]

vf = open("comparacio_bal.txt","w")
for i in range(7):
    x = i*1e-9
    m = np.mean(Teq[i])
    inc = 1.96*(np.var(Teq[i]))**0.5
    vf.write("%e %f %f %f\n"%(x, m, 0, inc))



#Ballistic and Difussive two-dimensional.
Tp = np.load("Tp.npy")

vf = open("comparacio_placa_d.txt","w")
for i in range(5):
    x = (2*i+1)*(1./(2**0.5))*1e-8
    m = np.mean([item[i,i] for item in Tp])
    inc = 1.96*(np.var([item[i,i] for item in Tp]))**0.5
    vf.write("%e %f %f %f\n"%(x, m, 0, inc))
    
    

#Diffusive one-dimensional.
Tcell_list = np.loadtxt("def_dif_8_500_llarg.txt", dtype='i')

Tmean_list = np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])
lon = 1000
for i in range(int((len(Tcell_list[0]))/lon)):
    s = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for j in range(lon*i,lon*(i+1)):
        s = s+Tcell_list[:,j]
    Tmean_list = np.append(Tmean_list, [s/lon], 0)
     
Tmean_list = np.delete(Tmean_list,0,0)

vf = open("Tmean_list_barra_difusiu.txt","w")
for i in range(7):
    for j in range(len(Tmean_list)):
        vf.write("%f "%Tmean_list[j][i])
    vf.write("\n")
    
Tfin=Tcell_list[:,200:500]

vf = open("comparacio2.txt","w")
for i in range(12):
    x = i*1e-8
    m = np.mean(Tfin[i])
    inc = 1.96*(np.var(Tfin[i]))**0.5
    vf.write("%e %f %f %f\n"%(x, m, 0, inc))
    






