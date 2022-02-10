#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 14:09:35 2020

@author: abdoulayebaradji
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

#%%

#%%

# fonction qui calcule la derivé d'une fonction en utilisant la dérivation de fourier 

def Dn(N) : 
    M = np.zeros((N,N))
    for n in range(N) : 
        for j in range(N) : 
            if (j!=n) :
                M[n,j] = 0.5*(-1)**(n-j)*np.cos((n-j)*np.pi/N)/np.sin((n-j)*np.pi/N)
            else : 
                M[n,j] = 0
    return M


                
def approx_dev(f,N) : 
    return np.dot(Dn(N),f)

# question 1 

def CollocFourierDerTemps(phi,N,nu) : 
    F = nu*approx_dev(phi,N)
    u = F - phi
    G = approx_dev(u,N)
    return G

# question 2

# Méthode RK3 
def CollocFourierRK3Pas(dt,phi,N,nu) : 
    U = phi
    k1 = CollocFourierDerTemps(U,N,nu)
    U = U + (1/3)*dt*k1
    k2 = (-5/9)*k1 + CollocFourierDerTemps(U,N,nu)
    U = U + (15/16)*dt*k2
    k3 = (-153/128)*k2 +CollocFourierDerTemps(U,N,nu)
    return U + (8/15)*dt*k3





# definition de la fonction h_j(x) pour tout x


def alpha(k,N) : 
    if (k == N/2 or k == -N/2 ):
        alpha = 2
    else :
        alpha =1
    return alpha



def h(x, N) : 
    n0 = np.arange(-N/2,1+(N/2),1)
    xj = np.arange(0,np.pi*2,2*np.pi/N)
    lij = []
    for j in range(N) : 
        li = []
        for k in range(N+1) :
            li.append(np.exp(complex(0,n0[k])*(x - xj[j]))/(alpha(n0[k],N)))
        lij.append(1/N*(np.sum(li)))
    return lij

#hh = h(1,16)

# definition de la foncion PHI en tout point 
 
def CollocFourierCalculApprox(x,vect0,dt,N,nu,T) :
    li = []
    for o in range(N) :
        li.append(CollocFourierRK3Pas(dt,vect0,N,nu)[o]*h(x,N)[o])
    return np.sum(li)


# question 3 

def CollocFourier(nu,T,N,dt,phi0,M) :
    xj_nout = np.arange(0,np.pi*2,2*np.pi/M)
    col = CollocFourierRK3Pas(dt,phi0,N,nu)
    no = int(T/(dt))
    li = []
    li.append(col)
    for o in range(1,no) :
        phi0 = col
        u = CollocFourierRK3Pas(dt,phi0,N,nu)
        li.append(u)
        col = u
    val_co = []
    for j in xj_nout :
        val_co.append(CollocFourierCalculApprox(j,li[no-1],dt,N,nu,T))
    return val_co
        
# question 4
# a.) interpolation : (erreur d'aliasing ) 

def phio(x) : 
    return 3/(5 - 4*np.cos(x))

# b.) 


# test 
popo = CollocFourier(0.2,2,16,1.25*10**(-3),phio(np.arange(0,np.pi*2,2*np.pi/16)),100)

# t = 0, 1 , 2
plt.plot(popo)
plt.show()




######### vraie solution ############################ 

def vrai_sol(t,nu,N) : 
    l = np.arange(-N/2,1 + N/2,1)
    xj = np.arange(0,np.pi*2,2*np.pi/N)
    u = []
    for i in xj : 
        uu = []
        for k in l : 
            uu.append(2**(-np.abs(k))*np.exp(complex(0,k*(i - t)) - nu*k**2*t))
        u.append(sum(uu))
    return u

# t=0
plt.plot(vrai_sol(0,1,50))
plt.show()

# t=1
plt.plot(vrai_sol(1,1,16))
plt.show()

# t=2
plt.plot(vrai_sol(2,1,16))
plt.show()

# erreur dans une norme discrète

# definition de la fonction alpha_k



#hhh = CollocFourierCalculApprox(1,phio(np.arange(0,np.pi*2,2*np.pi/16)),1.25*10**(-3),16,1,0)


def erreur_colloc(T,M,N,vect0,dt,nu) :
    err = (2*np.pi/M)*np.sum((np.array(vrai_sol(T, nu, M)) - np.array(CollocFourier(nu,T,N,dt,vect0,M)))**2)
    return err


er_colloc = erreur_colloc(2,100,50,phio(np.arange(0,np.pi*2,2*np.pi/50)),5*10**(-3),0.2)

Nnn = np.arange(6,51,2)
liii = []
for p in Nnn :
    liii.append(erreur_colloc(2,100,p,phio(np.arange(0,np.pi*2,2*np.pi/p)),5*10**(-3),0.2))



plt.plot(liii)
plt.show()


































##### Partie 2 ###################################

# création matrice de projection
# question 1 

# question a.) 

# 
def chap_init(N) :
    vect=[]
    li=np.arange(-N/2,1+(N/2),1)
    for k in li :
        vect.append(2**(-np.abs(k)))
    return vect




def DerConvectionDiffusionTemos(nu,vect0,N) : 
    li = np.arange(-N/2,1+(N/2),1)
    vect = []
    for k in range(N+1) :
        vect.append(-(complex(nu*li[k]**2,li[k]))*vect0[k])
    return vect

DerConvectionDiffusionTemos(0.2,chap_init(16),16)

# question b.) 

# un pas 
def PasMethodeGalerkinFourier(nu,dt,vect0,N) : 
    U = vect0
    k1 = np.array(DerConvectionDiffusionTemos(nu,U,N))
    U = U + k1*(1/3)*dt
    k2 = k1*(-5/9) + np.array(DerConvectionDiffusionTemos(nu,U,N))
    U = U + k2*(15/16)*dt
    k3 = np.dot(k2,(-153/128)) + np.array(DerConvectionDiffusionTemos(nu,U,N))
    return U + np.dot(k3,(8/15)*dt)

PasMethodeGalerkinFourier(0.2,5*10**(-3),chap_init(16),16)

# question c.)

# calcul la valeur dans un point x

def CalculApproxGalerkinFourier(nu,x,vect0,N,dt) : 
    n0 = np.arange(-N/2,1+(N/2),1)
    li = []
    for n in range(N+1) :
        li.append(PasMethodeGalerkinFourier(nu,dt,vect0,N)[n]*np.exp(complex(0,n0[n]*x)))
    return np.sum(li)


     
# question 2 


def FourierGalerkin(nu,N,dt,T,vect0,N_out) :
    no = int(T/(dt))
    xj_nout = np.arange(0,np.pi*2,2*np.pi/N_out)
    #xj_nout = np.arange(-np.pi,np.pi+(2*np.pi/N_out),2*np.pi/N_out)
    gal = PasMethodeGalerkinFourier(nu,dt,vect0,N) 
    li = []
    li.append(gal)
    for n in range(1,no) :
        vect0 = gal
        u = PasMethodeGalerkinFourier(nu,dt,vect0,N)
        li.append(u)
        gal = u
    val = []
    for j in xj_nout :
        val.append(CalculApproxGalerkinFourier(nu,j,li[int(T/dt)-1],N,dt))
    return val
                   


# t=2

test0 = FourierGalerkin(0.2,16,5*10**(-3),2,chap_init(16),16)


plt.plot(test0)
plt.show()

# vrai solution graphe galerkin

# ne change pas comme pour collocation

# erreur galerkin


def erreur_galerkin(T,M,N,vect0,dt,nu) :
    err = (2*np.pi/M)*np.sum((np.array(vrai_sol(T, nu, M)) - np.array(FourierGalerkin(nu,N,dt,T,vect0,M)))**2)
    return err


er = erreur_galerkin(2,50,16,chap_init(16),2.5*10**(-3),0.2)

# dt1
Nn = np.arange(6,30,2)
lii = []
for p in Nn :
    lii.append(erreur_galerkin(2,50,p,chap_init(p),5*10**(-3),0.2))

#dt2
Nn = np.arange(6,30,2)
liii = []
for p in Nn :
    liii.append(erreur_galerkin(2,50,p,chap_init(p),2.5*10**(-3),0.2))

plt.plot(liii)
plt.show()

#dt3
Nn = np.arange(6,30,2)
luu = []
for p in Nn :
    luu.append(erreur_galerkin(2,50,p,chap_init(p),1*10**(-4),0.2))

plt.plot(luu)
plt.show()

# graphe erreur pour t=2 avec N=16......

plt.plot(np.log(Nn),np.log(lii))
plt.show()


# superposition des graphes N=16, N_out=50,
# t=0
p1= plt.plot(vrai_sol(0,0.2,50),marker='o',label="vraie solution")
p2= plt.plot(FourierGalerkin(0.2,16,5*10**(-3),0,chap_init(16),50),marker='v',label="solution approchée")
plt.title("Rreprésentation graphique fouriergalerkin t=0")
plt.legend()
plt.show()

# t=1
p12= plt.plot(vrai_sol(1,0.2,50),marker='o',label="vraie solution")
p23= plt.plot(FourierGalerkin(0.2,16,5*10**(-3),1,chap_init(16),50),marker='v',label="solution approchée")
plt.title("Rreprésentation graphique fouriergalerkin t=1")
plt.legend()
plt.show()


# t=2
p11= plt.plot(vrai_sol(2,0.2,50),marker='o',label="vraie solution")
p22= plt.plot(FourierGalerkin(0.2,16,5*10**(-3),2,chap_init(16),50),marker='v',label="solution approchée")
plt.title("Rreprésentation graphique fouriergalerkin t=2")
plt.legend()
plt.show()


# erreur galerkin en fonction de N et differentes valeurs dt
ep= plt.plot(lii,marker='o',label="dt = 5*10**(-3)")
epp=  plt.plot(liii,marker='v',label="dt = 2.5*10**(-3)")

eppp= plt.plot(luu,marker='+',label="dt = 1.25*10**(-3)")



plt.title("Rreprésentation graphique de l'erreur fouriergalerkin t=2")
plt.legend()
plt.show()







#%%