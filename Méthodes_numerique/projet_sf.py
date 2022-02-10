#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 21:34:51 2020

@author: abdoulayebaradji
"""
#%%
import numpy as np
import os
import pandas as pd

spam = pd.read_table("/Users/abdoulayebaradji/Desktop/Stage_DIFUSE/data/commentaires.txt", sep="", header=0)

#%%

#%%
# interpolant des fonctions
import numpy as np
import matplotlib.pyplot as plt
def f1(x,h) : 
    return h*np.sin(np.pi*x/h)/(np.pi*x)
L=10
h=1
x = np.arange(-L,L+1,h)
plt.plot(x,f1(x,h))
plt.grid()
plt.show()

#%%


#%%
y = np.arange(-L/4,L/4+1,1)
q=[]
def f2(x,h) :
    for m in y :
        #print(m)
        q.append(f1(x-m*h,h))
        #print (q)
    return sum(q)
z = f2(x,h)
    
plt.plot(x,z)
plt.show()

#%%

#%%

k1 = np.arange(-L/3+1,1,1)
k2 = np.arange(1,L/3+1,1)
b=[]
def f3_bis1(x,h) : 
    for m in k1 :
        b.append( 9*f1(x-m*h,h) / (3*L*m*h + L**2))
    return sum(b)
i1 = f3_bis1(x,h)

v=[]
def f3_bis2(x,h) : 
    for m in k2 : 
        v.append(f1(x-m*h,h)*(L*m*h-3)/(L*m*h))
    return sum(v)
i2 = f3_bis2(x,h)    


i3 = i1 + i2

plt.plot(x,i3)
plt.show()
#%%


#%%
# DFT Algo 
x_no= np.arange(0 , 2 * np.pi, 2*np.pi/24)  
def DFT(g,N) :
    v_new=[]
    h = 2*np.pi/N
    k = np.arange(-N/2+1,N/2+1,1)
    for i in k : 
        v=[]
        for l in range(1,N+1) :
            v.append(h*np.exp(complex(0,-i*l*h))*g(l*h))
            #print (v)
        v_new.append(h*sum(v))
        #print(i)
        #print(v_new)
    return v_new


def DFT_bis(g,N) : 
    k = np.arange(-N/2+1,N/2,1)
    lala=[]
    for m in range(N-1) : 
        s = complex(k[m]*(DFT(g,N)[m].imag),k[m]*(DFT(g,N)[m].real))
        lala.append(s)
    lala.append(0)
    return lala


def DerSpecPerDFT (g,N) : 
    k=np.arange(-N/2+1,N/2+1,1)
    h = 2*np.pi/N
    u_new=[]
    for l in range(1,N+1) :
        u=[]
        for o in range(N) : 
            u.append(np.exp(complex(0,k[o]*l*h))*DFT_bis(g,N)[o])
        u_new.append(1/(np.pi*2)*sum(u))
    return u_new




# test avec les fonctions g1 et g2  : 

def g1(x) : 
    return max(0,abs(x-np.pi)/2)


def g2(x) : 
    return np.exp(np.sin(x))

zaza=DerSpecPerDFT(g1,24)
plt.plot(x_no,zaza)
plt.show()

zozo= DerSpecPerDFT(g2,24)
plt.plot(x_no,zozo)
plt.show()

#%%

#%%


#  MATRICE DE DERIVÉE

def mat_deriv(N) :
    h = 2*np.pi/N
    M = np.zeros((N,N))
    for i in range(0, N):              
      for j in range(0, N):            
         if j!=i:
             M[i,j] = 0.5 * pow(-1,i-j) * np.cos((i-j) * h / 2)/ np.sin((i-j) * h / 2)
    return M




def DerSpecPer(vecteur) :
    M = mat_deriv(len(vecteur))
    return np.dot(M,vecteur)


def g3(x):
    return np.cos(x) * np.exp(np.sin(x))

  
plt.plot(DerSpecPer(g1(x_no)))
plt.show()



#%%


#%%
# illustration graphique des erreurs 

# definitions des fonctions 
def q1(x) : 
    return np.abs(np.sin(x))**3

def q2(x) : 
    return np.exp(-(np.sin(x/2))**(-2))

def q3(x) : 
    return 1/(1 + (np.sin(x/2))**2)

def q4(x) : 
    return np.sin(10*x)

# derivé des fonctions 

def dq1(x) : 
    return -3*np.cos(x)*(np.sin(x))**2


def dq2(x) : 
    return np.cos(x/2)*(np.sin(x/2))**(-3)*np.exp(-(np.sin(x/2))**(-2))


def dq3(x) : 
    return -(np.cos(x/2)*np.sin(x/2))/(1+(np.sin(x/2))**2)**2

def dq4(x) : 
    return 10*np.cos(10*x)

# norme infinie d'un vecteur 

def norm_inf(vecteur) :
    nr= - 1.0
    for i in range(len(vecteur)) :
        if ( np.abs(vecteur[i]) > nr ) : 
            nr = np.abs(vecteur[i])
    return nr  
# calcul de l'erreur
#import math  
def erreur(g,dg,iter):
    vect = []
    for n in range(1,iter+1) : 
        h=2*np.pi/n
        x=np.arange(0,2*np.pi,h)
        s= DerSpecPer(g(x))
        u= dg(x)
        t=u-s
        U=norm_inf(t)
        vect.append(U)
    return vect
 
 
# graphe erreur fonction q1
plt.plot(erreur(q1,dq1,50))
plt.show()

# graphe erreur fonction q2
plt.plot(erreur(q2,dq2,50))
plt.show()


# graphe erreur fonction q3
plt.plot(erreur(q3,dq3,50))
plt.show()


# graphe erreur fonction q4
plt.plot(erreur(q4,dq4,50))
plt.show()




#%%


#%%
# valeurs propres d'un oscillateur harmonique :
    
# construction de la matrice de dérivation spectrale d'ordre 2 
L1=8
N1=6

def mat_deriv_2(N) : 
    M=np.zeros((N,N))
    h = 2*np.pi/N
    for n in range(N) : 
        for j in range(N) : 
            if (j==n) : 
                M[n,j] = -(np.pi)**2/(3*(h)**2) - 1/6 
            else :
                M[n,j] = -(-1)**(j+1)/(2*(np.sin((j+1)*h/2))**2)
    return M

mat_deriv_2(6)

x1 = np.arange(0,2*np.pi,2*np.pi/N1)
    

S = np.diag(x1)
# valeur propre 
import numpy.linalg as lng
B =  S - mat_deriv_2(N1)
val_p = lng.eig(B)[0]




#%%




















