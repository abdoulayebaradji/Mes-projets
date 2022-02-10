#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 13:36:52 2020

@author: abdoulayebaradji
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
#%%


#%%
################# exercice 1 ###########################
# ///// coéffcient de la dérivée de la série tronquée

#### quadrature de gauss ###########

# calcul des noeuds et poids 

def GT_noeuds(N) :
    xj=np.zeros(N+1)
    for j in range(N+1):
        xj[j]=np.cos(((2*j+1)/(2*N+2))*np.pi)
        #wj[j]=np.pi/(N+1)
    return xj

def GT_poids(N) :
    wj=np.zeros(N+1)
    for j in range(N+1):
        #xj[j]=np.cos(((2*j+1)/(2*N+2))*np.pi)
        wj[j]=np.pi/(N+1)
    return wj

#### quadrature de gauss lobato ###########

# calcul des noeuds et poids 
def GLT_noeuds(N) :
    xj=np.zeros(N+1)
    for j in range(N+1):
        xj[j]=np.cos(j*np.pi/N)
    return xj

def GLT_poids(N) :
    wj=np.zeros(N+1)
    for j in range(N+1):
        if (j==0 or j==N):
            wj[j]=np.pi/(2*N)
        else:
            wj[j]=np.pi/N
    return wj


# definition des N+1 coéffients de chebychev
### Gauss 

### GaussLobatto

def CoefChebyGL(f) : 
    N=len(f)-1
    k = np.arange(0,N+1,1)
    vect = []
    for i in k :
        vect0=[]
        for j in range(N+1) : 
            if (j==0 or j==N):
                vect0.append(0.5*f[j]*np.cos(i*(np.arccos(GLT_noeuds(N)[j]))))
            else:
                vect0.append(f[j]*np.cos(i*(np.arccos(GLT_noeuds(N)[j]))))
        if (i==0):
            vect.append((1/N)*np.sum(vect0))
        else:
            vect.append((2/N)*np.sum(vect0))
    return vect
    
def f(x) : 
    return np.cos(x)*0.5

def df(x) : 
    return -0.5*np.sin(x)

plt.plot(df(GLT_noeuds(20)))
plt.show()
  
    
#f(GLT_noeuds(4))
    
    
def CoefChebyGauss(f) : 
    N=len(f)-1
    k = np.arange(0,N+1,1)
    vect = []
    for i in k :
        vect0=[]
        for j in range(N+1) : 
            vect0.append(f[j]*GT_poids[j]*np.cos(i*(np.arccos(GT_noeuds(N)[j]))))
        if (i==0):
            vect.append((1/np.pi)*np.sum(vect0))
        else:
            vect.append((2/np.pi)*np.sum(vect0))
    return vect

                                
        
def CoefDerivéeCheby(f) :
    N=len(f)-1
    #k=np.arange(N-1,-1,-1)
    vect=np.zeros(N)
    vect[N-1]=2*N*CoefChebyGL(f)[N]
    vect[N-2]=2*(N-1)*CoefChebyGL(f)[N-1]
    i=1
    while (2*i+2 < N ) :
        vect[N-(2*i+1)] = 2*((N-(2*i+1))+1)*CoefChebyGL(f)[N-(2*i+1)+1]+ vect[N-(2*i-1)]
        vect[N-(2*i+2)] = 2*((N-(2*i+2))+1)*CoefChebyGL(f)[N-(2*i+2)+1]+ vect[N-2*i]
        i+=1
    vect[1] = 2*2*CoefChebyGL(f)[2]+vect[3]
    vect[0]= CoefChebyGL(f)[1]+vect[2]
    return vect

def DerivéeCheby(f,x) : 
    N= len(f)-1
    k=np.arange(0,N,1)
    vect=[]
    for j in k : 
        vect.append(CoefDerivéeCheby(f)[j]*np.cos(j*(np.arccos(x))))
    return np.sum(vect)
    
vv=[]
for i in range(20) : 
    vv.append(DerivéeCheby(f(GLT_noeuds(20)),GLT_noeuds(20)[i]))
    

plt.plot(vv)
plt.show()

pp1= plt.plot(df(GLT_noeuds(20)),marker='o',label="vraie solution")
pp2= plt.plot(vv,marker='v',label="solution approchée")
plt.title("Rreprésentation graphique de l'approximation de la derivée ")
plt.legend()
plt.show()




########## Test de l'exercie 1 #################
    


    
    
    
    
    
    

    

#%%


#%%
# exercice 2 

## definition du coéfficient an

def an(f) : 
    N=len(f)-1
    k=np.arange(0,N+1,1)
    vect=[]
    for i in k : 
        vect0=[]
        for j in range(N+1) :
            if (j==0 or j==N) : 
                vect0.append(0.5*f[j]*np.cos(j*i*np.pi/N))
            else:
                vect0.append(f[j]*np.cos(j*i*np.pi/N))
        if (i==0):
            vect.append((1/N)*np.sum(vect0))
        else:
            vect.append((2/N)*np.sum(vect0))
    return vect

an(f(GLT_noeuds(20)))

# coéffcient de la derivée 

def d_an(f) : 
    N=len(f)-1
    vect=np.zeros(N)
    vect[N-1] = 2*N*an(f)[N]
    vect[N-2] = 2*(N-1)*an(f)[N-1]
    i=1
    while 2*i +2 < N : 
        vect[N-(2*i+1)] = 2*(vect[N-(2*i-1)] + (N-2*i)*an(f)[N-2*i])
        vect[N-(2*i+2)] = 2*(vect[N-2*i] + (N-(2*i+1))*an(f)[N-(2*i+1)])
        i+=1
    vect[1] = 2*(2*an(f)[2]+vect[3])
    vect[0]= vect[2] +an(f)[1]
    return vect

d_an(f(GLT_noeuds(20)))


# calcul de la derivée par rapport à x 

def d_u(f) : 
    N=len(f)-1
    k = np.arange(0,N,1)
    vect=[]
    for j in range(N) :
        vect0=[]
        for n in k : 
            vect0.append(d_an(f)[n]*np.cos(np.pi*n*j/N))
        vect.append(np.sum(vect0))
    return vect


plt.plot(d_u(f(GLT_noeuds(20))))
plt.show()


# calcul du second membre de G

def U0(x) : 
    return x**2+x

def U0bis(x) : 
    if (x==-1):
        return -np.sin(5*np.pi*(-x-1))
    else:
        return 0


def ft(x,t) : 
    return x+1 +(2*x+1+t)*np.exp((x+1)*(x+t)+x)

def G(u,t,N) : 
    #N=len(f)-1
    vect=[]
    for j in range(N) : 
        vect.append(np.exp(u[j] + GLT_noeuds(N)[j])*d_u(ft(GLT_noeuds(N),t))[j] )
    return ft(GLT_noeuds(N-1),t) - vect

G(U0(GLT_noeuds(10)),0,10)


# Résolution de la méthode d'intégrations numérique Euler modifié  

def Euler_modifié(N,U0,T,dt) : 
    #N=len(f)-1
    k1 = U0 + np.dot(G(U0,0,N),dt/2)
    k2 = G(k1,0+(dt/2),N)
    U1 = U0 + np.dot(k2,dt)
    no = int(T/(dt))
    #res=np.zeros(len(U1))
    for j in range(1,no) : 
        U0=U1
        k1 = U0 + np.dot(G(U0,j,N),dt/2)
        k2 = G(k1,j+(dt/2),N)
        res = U0 + np.dot(k2,dt)
        U1 = res
    return res



ooo = Euler_modifié(10,U0(GLT_noeuds(9)),0.5,0.1)   

plt.plot(ooo)
plt.show()



# solution exacte


# RK4 

def RK3Pas_1(N,U0,dt,t) : 
    U = U0
    k1 = G(U,t,N)
    U = U + np.dot(k1,(1/3)*dt)
    k2 = np.dot(k1,(-5/9)) + G(U,t+(dt/3),N)
    U = U + np.dot(k2,dt*(15/16))
    k3 = np.dot(k2,(-153/128)) + G(U,t+(3*dt/4),N)
    return U + np.dot(k3,dt*(8/15))


def RK3Pas(N,U0,T,dt) :
    no = int(T/(dt))
    gal = RK3Pas_1(N,U0,dt,0) 
    li = []
    li.append(gal)
    for n in range(1,no) :
        U0 = gal
        u = RK3Pas_1(N,U0,dt,n)
        li.append(u)
        gal = u
    return li

oo = RK3Pas(14,U0(GLT_noeuds(13)),0.5,0.1)   
o= oo[0] + oo[1] + oo[2] +  oo[3] + oo[4] 


plt.plot(o)
plt.show()




#%%



#%%

# exercice 3 




#%%








