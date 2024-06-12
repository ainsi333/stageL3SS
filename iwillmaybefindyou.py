import cv2
import numpy as np
import matplotlib.pyplot as plt
import random as rd

A=np.zeros((100,100))
A[95,40]=100
B=np.zeros((10,5))
B[2,4]=100

iwillmaybefindyou(A,B,20,100,D)

def D(A,pos,P2):  #distance entre 2 patch
    n,m=np.shape(P2)
    P1=A[pos[0]:pos[0]+n,pos[1]:pos[1]+m]
    return np.sqrt(np.sum(np.abs((P1-P2)**2)))

def recherche(A,hb,wb,a,pos_visite,pos,vois):
    ha,wa=A.shape
    n=ha-hb+1
    m=wa-wb+1
    c=0 #compter
    T=[] #patch a explorer, un patch est representé par sa case en haut a gauche
    if pos[0]<=ha-2*a-1 and pos[1]<=wa-2*a-1:
        i=0
        j=0
        while pos[0]+i<ha and pos[1]+i<wa and i<2*a+1 and j<2*a+1:
            if pos_visite[pos[0]+i,pos[1]+j]==0:
                c+=1
                T.append([pos[0]+i,pos[1]+j])
    if c<(a+1/2)**2: #bcp de cases non explorées, patchmatch
        return "patchmatch"
    else:  #on va chercher a la main
        while c<(2*a+1)**2:
            for v in vois:
                if pos_visite[v[0],v[1]]==1 and v not in T:
                    T.append(v)
                    c+=1
        return T
def fait_la_recherche(resultat,A,B,D,S,P,seuil):
    #if type(resultat)==str:
     #   "a faire" #patchmatch jsp
    if 1==2:
        return 8
    else:
        T=resultat
        if len(T)!=0:
            n=len(T)
            min=D(A,T[0],B)
            meilleur=0
            for i in range(n):
                d=D(A,T[i],B)
                if d<min:
                    meilleur=i
                    min=d
                if d<seuil:
                    S.append([d,T[i]])
                elif len(S)==0 and d>=seuil and d<P[1]: #mettre a jour P
                    P=[T[i],d]




def iwillmaybefindyou(A,B,a,seuil,D):
    ha,wa=A.shape
    hb,wb=B.shape
    n=ha-hb+1
    m=wa-wb+1
    S=[]   #Positions sous le seuil
    P=[[0,0],D(A,[0,0],B)]  #Position qui minimise D si aucune n'est sous le seuil
    T=100
    pos_visite=np.zeros((n,m))  #patch deja visite, representes par le pixel en haut a gauche du patch
    # début position initiale
    init=[[rd.randint(0,n-1),rd.randint(0,m-1)]]
    while len(init)<3:
        pos=[rd.randint(0,n-1),rd.randint(0,m-1)]
        if pos not in init:
            init.append(pos)
    phase=1 #phase 1: recherche, phase 2: déplacement
    pos=init[0]
    for i in range(3):
        if D(A,init[i],B)<D(A,pos,B):
            pos=init[i]
     # fin position initiale
    for l in range(T):
        # phase 1 : on fait la recherche dans un domaine de rayon a
        x0,y0 = pos
        vois = []
        for j in range(-a,a+1):
            for i in range(-a,a+1):
                vois += [[x0+i,y0+j] for i in range(-a,a+1)]
        vois = np.array(vois)
        fait_la_recherche(recherche(A,hb,wb,a,pos_visite,pos,vois),A,B,D)
        # phase 2 : maintenant on tire la nouvelle position et on change la phase
        newpos=[]
        while len(newpos)==0:
            test=[rd.randint(0,n-1),rd.randint(0,m-1)]
            if pos_visite[test[0],test[1]]==0 :
                newpos = test
        pos = newpos

    S = sorted(S)
    return S,P







# np.where