# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:14:42 2024

@author: imtey
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate as intg
import math as m
from scipy.optimize import fsolve
from mpl_toolkits.mplot3d import Axes3D



## essai 1D

"fcts pratiques/test"
#1D
def trivial(x):
    if x<=1 and x>=0:
        return 1
    else:
        return 0
def iden(x):  #id
    return x
def arghh(x): #e^-x
    return m.exp(-x)
def indicplus(x): #indicatrice de R+ * exp(-x)
    if x < 0 : 
        return 0
    else : 
        return arghh(x)
def triangle(x):
    if x>=0 and x<=1:
        return x
    elif x>1 and x<=2: 
        return 2-x
    else: 
        return 0
sigma=0.5
mu=3
def normal(bins):
    return 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu)**2 / (2 * sigma**2))


def dégueu(x):   #les pics en 1/n^2
    n = int(x)
    if x <= 0.5 : 
        return 0
    elif x >= n+1 - 1/(2*(n+1)**3):
        return (2*((n+1)**4)*x +(n+1) - 2*((n+1)**5))*12/(np.pi)**2
    elif x <= n + 1/(2*n**3):
        return (-2*(n**4)*x +n+2*n**5 )*12/(np.pi)**2
    else : return 0

def coscos(x):
    if x<1:
        return 0
    else:
        return (1+np.cos(x))/x**2*1/(0.915702)
#2D
def indic(x,y):
    if 0<x<1 and 0<y<1:
        return 1
    else:
        return 0
def normal_superieur(bins,gims):
    return 1/(2 * np.pi) *np.exp( - (np.sqrt(bins**2+gims**2) -mu)**2 / (2))
def indicplus2(x,y):
    return np.exp(-np.sqrt(x**2+y**2)) *1/(6.283185307171164)
            
"fcts du modèle"
#1D
def P(p,phi):
    def f(x):  #le truc ds l'integrale
        if p(x) == 0 or phi(x) == 0 : 
            return 0
        else : 
            return p(x)*(1-m.exp(-phi(x)))  #retourne aussi l'erreur
    return intg.quad(f,-m.inf,m.inf)
def P2(p,phi):
    def f(x,y):  #le truc ds l'integrale
        if p(x,y) == 0 or phi(x,y) == 0 : 
            return 0
        else : 
            return p(x,y)*(1-m.exp(-phi(x,y)))  #retourne aussi l'erreur
    return intg.dblquad(f,-m.inf,m.inf,lambda x:-m.inf,lambda x:m.inf)
def Popti(p,Gphi):  #pour phi optimal
    lambd=iwillfindyou(p, Gphi)
    def phi(t):
        if p(t)>lambd:
            return m.log(p(t))-m.log(lambd)
        else:
            return 0
    return P(p,phi)
def Popti2(p,Gphi):  #pour phi optimal en 2D
    lambd=iwillfindyou2(p, Gphi)
    def phi(s,t):
        if p(s,t)>lambd:
            return m.log(p(s,t))-m.log(lambd)
        else:
            return 0
    return P2(p,phi)

# résolution 
# def iwillfindyou(distrib_cible,Gphi):  #distrib:p, Gphi:qte effort
#     def searchlambda(Gphi):
#         def searchf(u): #integrale de phiu
#             def phiu(u): #log(p(t)/u) 
#                 def lafonction(t):
#                     condition= (indicplus(t) >u)
#                     # if distrib_cible(t)-u>0:
#                     #     return m.log(distrib_cible(t))-m.log(u)
#                     # else:
#                     #     return 0
#                     np.where(condition,np.log(indicplus(t))-np.log(u),0)
#                 return lafonction
#             return intg.quad(phiu(u),-np.Inf,np.Inf)[0]-Gphi
#         borneinf=0.1
#         bornesup=1 # pour le fsolve
#         while searchf(borneinf)<0:
#             borneinf=borneinf/10
#         while searchf(bornesup)>0:
#             bornesup=bornesup*10
#         return fsolve(searchf,[borneinf])[0]
def iwillfindyou(distrib_cible,Gphi):  #distrib:p, Gphi:qte effort ILMARCHEEEE?
    def phiu(u): #log(p(t)/u) 
        def lafonction(t):
            if distrib_cible(t)>u:
                return m.log(distrib_cible(t))-m.log(u)
            else:
                return 0
        return lafonction
    def searchf(u): #integrale de phiu
        return intg.quad(phiu(u),-m.inf,m.inf)[0]-Gphi
    def searchlambda(Gphi):
        borneinf=0.1
        bornesup=1 # pour le fsolve
        while searchf(borneinf)<0:
            borneinf=borneinf/10
        while searchf(bornesup)>0:
            bornesup=bornesup*10
        return float(fsolve(searchf,[borneinf])[0])
    return searchlambda(Gphi)
    

def trace_phi(distrib_cible,Gphi):
    X=np.linspace(-0.5,2.5,1000)
    lambd=iwillfindyou(distrib_cible,Gphi)
    def phi(t):
        if distrib_cible(t)>lambd:
            return m.log(distrib_cible(t))-m.log(lambd)
        else:
            return 0
    Y=[phi(x) for x in X]
    Z=[distrib_cible(x) for x in X]
    Lambda=[lambd for x in X]
    abel="lambda="+str(0.2231)
    plt.clf()
    #plt.figure().set_figwidth(14)
    #plt.rcParams['font.size'] = 25
    plt.plot(X,Y,label="phi")  #bleu
    plt.plot(X,Z,label="p")  #orange
    #plt.plot(X,Lambda,label=abel) #vert
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend(loc="upper left")
    plt.show()
    return lambd

    
#2D
def iwillfindyou2(distrib_cible,Gphi):  #distrib:p, Gphi:qte effort ILMARCHEEEE?
    def phiu(u): #log(p(t)/u) 
        def lafonction(s,t):
            if distrib_cible(s,t)>u:
                return m.log(distrib_cible(s,t))-m.log(u)
            else:
                return 0
        return lafonction
    def searchf(u): #integrale de phiu
        return intg.dblquad(phiu(u),-m.inf,m.inf,lambda x:-m.inf,lambda x:m.inf)[0]-Gphi
    def searchlambda(Gphi):
        borneinf=0.1
        bornesup=1 # pour le fsolve
        while searchf(borneinf)<0:
            borneinf=borneinf/10
        while searchf(bornesup)>0:
            bornesup=bornesup*10
        return float(fsolve(searchf,[borneinf])[0])
    return searchlambda(Gphi)

def trace_phi2(distrib_cible,Gphi):
    def phi(s,t):
        if distrib_cible(s,t)>lambd:
            return m.log(distrib_cible(s,t))-m.log(lambd)
        else:
            return 0
    lambd=iwillfindyou2(distrib_cible,Gphi)
    x = np.linspace(-5, 5, 100)
    y = np.linspace(-5, 5, 100)
    X, Y = np.meshgrid(x, y)
    #Z = phi(X, Y)
    Z = np.array([[phi(a,b) for a in x] for b in y])
    #Proba=np.array([[distrib_cible(a,b) for a in x] for b in y])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Graphique de phi(x, y)')
    #ax.plot_surface(X, Y, Proba, cmap='viridis')
    ax.mouse_init()
    plt.show()

X = np.linspace(-0.5, 2.5, 1000)
p=[triangle(x) for x in X]
plt.clf()
plt.plot(X,p,label="p",color="orange")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc="upper left")












