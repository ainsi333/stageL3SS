import numpy as np
import matplotlib.pyplot as plt
from random import *
from math import *


def tc(rayon,x0,y0):
	theta = np.linspace(0, 2*np.pi, 100)  # Angles de 0 à 2*pi

	xc =x0+ rayon * np.cos(theta)
	yc =y0+ rayon * np.sin(theta)

	return xc,yc


def intermittentSearch(x_i,y_i,phasei,x_0,y_0,tau1,tau2,D,V,a,b,N) : # fonction qui va créer la zone de recherche, et effectuer la recherche en temps réel (mettre b=0 si on ne veut chercher qu'une seule cible (correspond à un espace de recherche infini))
	# création zone de recherche avec cible
	# pas&temps
	k = 0.1
	T = 100

	# liste des temps
	vt = np.arange(0,T,k)

	# liste des tau tirés
	t1 = [0]+[np.random.exponential(tau1) for i in range(len(vt))]
	t2 = [0]+[np.random.exponential(tau2) for i in range(len(vt))]



	# lancement
	x_init = x_i
	y_init = y_i
	x_sol = [x_init]
	y_sol = [y_init]
	c1 = 0
	c2 = 0
	t=0
	if phasei == 1 :
		direc = np.random.uniform(-np.pi,np.pi)
		for i in range(1,int(t1[1])+1) :
			x_sol.append(x_i + k*V*i*cos(direc))
			y_sol.append(y_i + k*V*i*sin(direc))
		for i in range(len(x_sol)):
			plt.clf()
			#plt.axis([-20+x_i,20+x_i,-20,20]) # le remettre partout pour que l'écran bouge avec le chercheur
			if b == 0 :
				plt.plot(x_0,y_0,marker='o',label="cible")
			else :
				N=[]
				kx = floor((x_sol[i]-x_0)/(2*b))
				ky = floor((y_sol[i]-y_0)/(2*b))
				N.append([2*kx*b+x_0,2*ky*b+y_0])
				N.append([2*(kx+1)*b+x_0,2*ky*b+y_0])
				N.append([2*(kx+1)*b+x_0,2*(ky+1)*b+y_0])
				N.append([2*kx*b+x_0,2*(ky+1)*b+y_0])
				N = np.array(N)
				plt.scatter(N[:,0],N[:,1],marker='o',label="cible")
			plt.plot(x_init,y_init,marker='x',label="position de départ")
			plt.plot(x_sol[i],y_sol[i], marker="x",label="chercheur")

			plt.plot(x_sol[:i],y_sol[:i])
			plt.title("recherche stochastique au temps t = {}".format(i)+" dans la phase 1")
			# Tracer le cercle en pointillé
			for j in range(len(N)) :
				plt.plot(tc(a,N[j][0],N[j][1])[0], tc(a,N[j][0],N[j][1])[1], linestyle='dashed', linewidth=2, color='blue', label='Cercle en pointillé')

			plt.pause(10**(-3))
		c1 += 1 #  nombre d'apparitions de la phase 1
		x_i,y_i = x_sol[-1],y_sol[-1]
		t+=int(t1[1])+1
		#phasei = 2
	else :
		x_sol = x_sol + diffusion(x_i,y_i,t2[1],k,D)[0]
		y_sol = y_sol + diffusion(x_i,y_i,t2[1],k,D)[1]
		for i in range(len(x_sol)):
			plt.clf()
			#plt.axis([-20,20,-20,20])
			if b == 0 :
				plt.plot(x_0,y_0,marker='o',label="cible")
			else :
				N=[]
				kx = floor((x_sol[i]-x_0)/(2*b))
				ky = floor((y_sol[i]-y_0)/(2*b))
				N.append([2*kx*b+x_0,2*ky*b+y_0])
				N.append([2*(kx+1)*b+x_0,2*ky*b+y_0])
				N.append([2*(kx+1)*b+x_0,2*(ky+1)*b+y_0])
				N.append([2*kx*b+x_0,2*(ky+1)*b+y_0])
				N = np.array(N)
				plt.scatter(N[:,0],N[:,1],marker='o',label="cible")
			plt.plot(x_init,y_init,marker='x',label="position de départ")
			plt.plot(x_sol[i],y_sol[i], marker="x",label="chercheur")
			# Tracer le cercle en pointillé
			for j in range(len(N)) :
				plt.plot(tc(a,N[j][0],N[j][1])[0], tc(a,N[j][0],N[j][1])[1], linestyle='dashed', linewidth=2, color='blue', label='Cercle en pointillé')
			for x in N:
				if np.linalg.norm(np.array([x_sol[i],y_sol[i]])-np.array(x))<=a :
					return "cible trouvée en {} étapes !".format(i)
			plt.plot(x_sol[:i],y_sol[:i])
			plt.title("recherche stochastique au temps t = {}".format(i)+" dans la phase 2")

			plt.pause(10**(-3))
		c2 += 1 #compte nombre d'apparitions de la phase 2
		x_i,y_i = x_sol[-1],y_sol[-1]
		t+= int(t2[1])+1
		#phasei = 1

	# afficher graphes successivement
	for l in range(len(vt)):

		if phasei == 2 :
			c1+=1
			direc = np.random.uniform(-np.pi,np.pi)
			x_sol = x_sol + [x_i + k*V*i*cos(direc) for i in range(int(t2[c2]),int(t2[c2])+int(t1[c1]))]
			y_sol =	y_sol + [y_i + k*V*i*sin(direc) for i in range(int(t2[c2]),int(t2[c2])+int(t1[c1]))]
			x_i,y_i = x_sol[-1],y_sol[-1]
			for i in range(t,t+int(t1[c1])) :

				# condition
				# if np.linalg.norm(np.array([x_sol[i],y_sol[i]])-np.array([x_0,y_0]))<=a :
				# 	return "cible trouvée en {} étapes !".format(i)
				plt.clf()
				#plt.axis([-20,20,-20,20])
				if b == 0 :
					plt.plot(x_0,y_0,marker='o',label="cible")
				else :
					N=[]
					kx = floor((x_sol[i]-x_0)/(2*b))
					ky = floor((y_sol[i]-y_0)/(2*b))
					N.append([2*kx*b+x_0,2*ky*b+y_0])
					N.append([2*(kx+1)*b+x_0,2*ky*b+y_0])
					N.append([2*(kx+1)*b+x_0,2*(ky+1)*b+y_0])
					N.append([2*kx*b+x_0,2*(ky+1)*b+y_0])
					N = np.array(N)
					plt.scatter(N[:,0],N[:,1],marker='o',label="cible")
				plt.plot(x_init,y_init,marker='x',label="position de départ")
				plt.plot(x_sol[i],y_sol[i], marker="x",label="chercheur")

				plt.plot(x_sol[:i],y_sol[:i])
				plt.title("recherche stochastique au temps t = {}".format(i)+" dans la phase 1")
				# Tracer le cercle en pointillé
				for j in range(len(N)) :
					plt.plot(tc(a,N[j][0],N[j][1])[0], tc(a,N[j][0],N[j][1])[1], linestyle='dashed', linewidth=2, color='blue', label='Cercle en pointillé')

				plt.pause(10**(-3))
				phasei = 1
				end = i
			t+=int(t1[c1])
			print("1-->2")



		else :
			c2+=1
			x_sol = x_sol + diffusion(x_i,y_i,t2[c2],k,D)[0]
			y_sol =	y_sol + diffusion(x_i,y_i,t2[c2],k,D)[1]
			x_i,y_i = x_sol[-1],y_sol[-1]
			for i in range(t,t+int(t2[c2])) :
				plt.clf()
				#plt.axis([-20,20,-20,20])
				if b == 0 :
					plt.plot(x_0,y_0,marker='o',label="cible")
				else :
					N=[]
					kx = floor((x_sol[i]-x_0)/(2*b))
					ky = floor((y_sol[i]-y_0)/(2*b))
					N.append([2*kx*b+x_0,2*ky*b+y_0])
					N.append([2*(kx+1)*b+x_0,2*ky*b+y_0])
					N.append([2*(kx+1)*b+x_0,2*(ky+1)*b+y_0])
					N.append([2*kx*b+x_0,2*(ky+1)*b+y_0])
					N = np.array(N)

					plt.scatter(N[:,0],N[:,1],marker='o',label="cible")
				plt.plot(x_init,y_init,marker='x',label="position de départ")
				plt.plot(x_sol[i],y_sol[i], marker="x",label="chercheur")
				plt.plot(x_sol[:i],y_sol[:i])
				plt.title("recherche stochastique au temps t = {}".format(i)+" dans la phase 2")
				plt.legend()

				# Tracer le cercle en pointillé
				for j in range(len(N)) :
					plt.plot(tc(a,N[j][0],N[j][1])[0], tc(a,N[j][0],N[j][1])[1], linestyle='dashed', linewidth=2, color='blue', label='Cercle en pointillé')
				# condition
				for w in N:
					if np.linalg.norm(np.array([x_sol[i],y_sol[i]])-np.array(w))<=a :
						return "cible trouvée en {} étapes !".format(i)

				plt.pause(10**(-3))
				end = i
				phasei = 2
			t+= int(t2[c2])
			print("2-->1")




	return "cible non trouvée en {} étapes...".format(end)


# Diffusion
# Paramètres de la simulation

def diffusion(x_start,y_start,T,k,D) :
	num_steps = int(T/k)  # Nombre d'étapes de la simulation
	dt = k  # Intervalle de temps entre chaque étape

	# Initialisation de la position
	x = x_start
	y = y_start

	# Listes pour stocker les positions au fil du temps
	x_positions = []
	y_positions = []

	# Simulation du mouvement diffusif en 2D
	for step in range(num_steps):
		# Calcul des déplacements aléatoires selon une distribution gaussienne
		dx = np.sqrt(2 * D * dt) * np.random.randn()  # Déplacement selon x
		dy = np.sqrt(2 * D * dt) * np.random.randn()  # Déplacement selon y

		# Mise à jour de la position
		x += dx
		y += dy

		# Enregistrement des nouvelles positions
		x_positions.append(x)
		y_positions.append(y)

	return x_positions,y_positions





