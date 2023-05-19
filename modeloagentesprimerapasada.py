# -*- coding: utf-8 -*-
"""
@author: Manuel
"""


#Programa para el estudio de primera pasada.

#Importamos las librerías que vamos a necesitar.
import random, math 
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate


#Esta es la función principal.
def fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Para que la división temporal sea lo suficientemente pequeña, debemos considerar al menos dos órdenes de magnitud 
    #respecto a los rates característicos.
    delta_t=abs(10**(-2)*min(1/kappa_0_to_1,1/kappa_1_to_0,1/k_A[0],1/k_A[1],1/k_B[0],1/k_B[1],1/pi_A_to_B,1/pi_B_to_A))
    #El tiempo define el número de iteraciones.
    tiempo=np.arange(delta_t,1000,delta_t)
    #Logaritmo de la población inicial.
    n_total=[math.log(n_A+n_B)]
    extincion=False #Como introducimos un número positivo de bacterias no hay extinción.
    tiempo_extincion=False #Ni tampoco tiempo de extinción.
    #Generemos la semilla para los números aleatorios.
    semilla=np.random.default_rng()
    #Elegimos un entorno inicial al azar.
    entorno_actual=semilla.choice([0,1],p=[0.5,0.5]) 
    for j in range (len(tiempo)):
        #¿Cambia el ambiente?
        entorno_actual=cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0, delta_t)
        #Estudiamos la reproducción y muerte de las bacterias.
        n_B=reproduccion_muerte_bacteria(n_B, k_B[entorno_actual], delta_t, semilla)
        n_A=reproduccion_muerte_bacteria(n_A, k_A[entorno_actual], delta_t, semilla)
        #Recorremos todas las bacterias que han sobrevivido y vemos si cambian el fenotipo.
        n_A, n_B=cambio_fenotipo(n_A,n_B,pi_A_to_B,pi_B_to_A,delta_t,semilla)
        #Estudiamos si hay extinción.
        if (n_A+n_B)<=1: #Como estamos con logaritmos, log 1=0 y hay extinción con una sola bacteria.
            extincion=True
            tiempo_extincion=tiempo[j] #Guardamos el tiempo de extinción.
            #Rellenamos el resto del vector y salimos del programa para ahorrar tiempo.
            for i in range(j, len(tiempo)):
                n_total+=[(0)]
            break
        n_total+=[math.log(n_A+n_B)]
        if (n_A+n_B)>100000000: #Cota superior.
            #Rellenamos el resto del vector y salimos del programa para ahorrar tiempo.
            for i in range(j, len(tiempo)-1):
                n_total+=[math.log(n_A+n_B)]
            break
    return n_A, n_B, n_total, np.insert(tiempo,0,0), extincion, tiempo_extincion #El insert simplemente es para añadir el origen de tiempos.
    
def cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0,delta_t):
    r=random.random() #Devuelve un número aleatorio entre 0 y 1.
    #Si estamos en el primer entorno y se cumple la condición...
    if entorno_actual==0 and r<=(kappa_0_to_1*delta_t):
        entorno_actual=1 #... pasamos al entorno 1.
    #Lo mismo para el otro.
    elif entorno_actual==1 and r<=(kappa_1_to_0*delta_t):
        entorno_actual=0 
    return entorno_actual

#Recordemos que n es el número de bacterias, k es la tasa de crecimiento en el entorno.
def reproduccion_muerte_bacteria(bacterias,k,delta_t,semilla):
    #Generamos un número aleatorio que siga una distribución binomial, que nos indicará el número de seres que se interactuan.
    r=semilla.binomial(n=bacterias,p=abs(k*delta_t))
    if k<0:
        #Se mueren ese número de seres.
        return bacterias-int(r) 
    if k>0:
        #Se reproducen ese número de seres.
        return bacterias+int(r)
        
def cambio_fenotipo(n_A,n_B,pi_A_to_B,pi_B_to_A,delta_t,semilla):
    #Para aseguarnos de que al recorrer un fenotipo no modifiquemos el otro instaneamente, los cambios de fenotipos
    #se guardan en unas variables extra.
    n_A_extra=0; n_B_extra=0
    #Generamos un número aleatorio siguiendo una distribución binomial. 
    r1=semilla.binomial(n=n_A,p=pi_A_to_B*delta_t)
    #Ese número nos da el número de bacterias del tipo A que han cambiado al B.
    n_B_extra+=int(r1)
    #Generamos un número aleatorio siguiendo una distribución binomial. 
    r2=semilla.binomial(n=n_B,p=pi_B_to_A*delta_t)
    #Ese número nos da el número de bacterias del tipo B que han cambiado al A.
    n_A_extra+=int(r2)
    #El número final de cada fenotipo, es el original más los que han cambiado del fenotipo contrario y los que se 
    #pierden porque han cambiado a su opuesto.
    return n_A+n_A_extra-n_B_extra, n_B+n_B_extra-n_A_extra


#La siguiente función devuelve el tiempo de primera pasada normalizada.
def tiempo_primera_pasada(t,x_0,mu,sigma,dominio):
    g=((x_0)/(sigma*math.sqrt(2*math.pi)*np.power(t,3/2))*np.exp(-(x_0+mu*t)**2/(2*sigma**2*t)))
    integral=integrate.trapz(g, dominio)
    return g/integral

#Esta función devuelve la extinción y nos permite comparar el tiempo de primera pasada teórico con el de las simulaciones.
def prob_extincion(numero_simulaciones,entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Creamos listas vacías para guardar los tiempos donde se produce la extinción.
    tiempo_extincion_vector=[]
    for i in range(numero_simulaciones): #Bucle en el número de simulaciones.
        #En este caso solo nos interesa la extinción.
        _, _, _, _, extincion,tiempo_extincion=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
        if extincion==True:
            #Si hay extinción, guardemos el tiempo.
            tiempo_extincion_vector+=[tiempo_extincion]
    #Calculamos la probabilidad de extinción.
    prob_extincion=len(tiempo_extincion_vector)/numero_simulaciones
    #En caso de que haya más de 1 un caso con extinción, tiene sentido representar un histograma.
    if len(tiempo_extincion_vector)>=2:
        intervalo=np.linspace(min(tiempo_extincion_vector),max(tiempo_extincion_vector))
        plt.hist(x=tiempo_extincion_vector, bins=intervalo,color = "lightblue", ec="red", density='True')
        #Introducimos los datos del punto obtenidos con el modelo colectivo.
        sigma=1.13
        mu=0.239
        #Representamos el tiempo de primera pasada.
        intervalo2=np.arange(0.01,max(tiempo_extincion_vector),0.05)
        primera_pasada=tiempo_primera_pasada(intervalo2,math.log(10), mu, sigma, intervalo2)
        plt.plot(intervalo2,primera_pasada, label='Tiempo de primera pasada teórico')
        #Personalizamos la gráfica.
        plt.title('Histograma de tiempos de extinción con $π_{A\_to\_B}=$'+str(pi_A_to_B)+',$π_{B\_to\_A}=$'+str(pi_B_to_A),fontsize = 13)
        plt.xlabel('Tiempos de extinción',fontsize = 14); plt.xticks(fontsize=14)
        plt.ylabel('Frecuencia',fontsize = 14); plt.yticks(fontsize=14)
        plt.legend(fontsize=12)
        #Y mostramos.
        plt.show()
        #Comprobamos que el tiempo de primera pasada integra aproximadamente 1.
        integral=integrate.trapz(primera_pasada, intervalo2)
    return prob_extincion, integral

#Definimos el entorno. El entorno 0 se corresponde con un ambiente normal mientras que el 1 uno con antibióticos.
entornos=[0,1]
#La tasa de cambio entre entornos viene dada por las kappas; kappa_0_to_1 para cambiar del entorno 0 al 1 y kappa_1_to_0
#para el proceso inverso.
kappa_0_to_1=1; kappa_1_to_0=1
#El número de bacterias las contabilizamos en las siguientes variables:
n_A=5; n_B=5
#Necesitamos también las tasas de crecimiento de cada fenotipo en cada ambiente. Los fenotipos son A y B.
k_A0=2; k_A1=-2; k_A=[k_A0,k_A1]; k_B0=0.2; k_B1=-0.2; k_B=[k_B0,k_B1]

#Este programa es solo para estudiar el tiempo de primera pasada.
pi_A_to_B=0.569; pi_B_to_A=0.207
prob_extincion, integral=prob_extincion(1000, entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A)
print('La probabilidad de extinción es:', prob_extincion)
print('El tiempo de primera pasada integra:', integral)