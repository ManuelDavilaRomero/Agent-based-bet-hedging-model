# -*- coding: utf-8 -*-
"""
@author: Manuel
"""


#Program for the stuyding of the first passage time.

import random, math 
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy.special import erf 
import pickle

#This is the main function.
def fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Definition of delta t.
    delta_t=abs(1000**(-1)*min(1/kappa_0_to_1,1/kappa_1_to_0,1/k_A[0],1/k_A[1],1/k_B[0],1/k_B[1],1/pi_A_to_B,1/pi_B_to_A))
    #Time defines the number of iterations.
    tiempo=np.arange(delta_t,1000,delta_t)
    #Logarithm of the initial population. 
    n_total=[math.log(n_A+n_B)]
    extincion=False #A positive number in the time origin: no extinction.
    tiempo_extincion=False #Nether there is an extinction time.
    #Seed for random generator (123456) for reproducibility if only 1 simulation is needed
    semilla=np.random.default_rng()
    #We choose an initial environment randomly.
    entorno_actual=semilla.choice([0,1],p=[0.5,0.5]) 
    random.seed()
    for j in range (len(tiempo)):
        #Does the environmet change?
        entorno_actual=cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0, delta_t)
        #Death and reproduction of bacteria.
        n_B=reproduccion_muerte_bacteria(n_B, k_B[entorno_actual], delta_t, semilla)
        n_A=reproduccion_muerte_bacteria(n_A, k_A[entorno_actual], delta_t, semilla)
        #Change of phenotype.
        n_A, n_B=cambio_fenotipo(n_A,n_B,pi_A_to_B,pi_B_to_A,delta_t,semilla)
        #Is there extinction?
        if (n_A+n_B)<=1: #Using logarithms, log 1=0 and one bacteria is equal to extinction.
            extincion=True
            tiempo_extincion=tiempo[j] 
            for i in range(j, len(tiempo)):
                n_total+=[(0)]
            break
        n_total+=[math.log(n_A+n_B)]
        if (n_A+n_B)>100000000: #Upper bound.
            for i in range(j, len(tiempo)-1):
                n_total+=[math.log(n_A+n_B)]
            break
    return n_A, n_B, n_total, np.insert(tiempo,0,0), extincion, tiempo_extincion 
    
def cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0,delta_t):
    r=random.random()
    if entorno_actual==0 and r<=(kappa_0_to_1*delta_t):
        entorno_actual=1 
    elif entorno_actual==1 and r<=(kappa_1_to_0*delta_t):
        entorno_actual=0 
    return entorno_actual

def reproduccion_muerte_bacteria(bacterias,k,delta_t,semilla):
    #Random number following a binomial distribution.
    r=semilla.binomial(n=bacterias,p=abs(k*delta_t))
    if k<0:
        #r bacteria die
        return bacterias-int(r) 
    if k>0:
        #r bacteria reproduce.
        return bacterias+int(r)
        
def cambio_fenotipo(n_A,n_B,pi_A_to_B,pi_B_to_A,delta_t,semilla):
    #We use aux variable to avoid modifying the original.
    n_A_extra=0; n_B_extra=0
    #Random number following a binomial distribution.
    r1=semilla.binomial(n=n_A,p=pi_A_to_B*delta_t)
    n_B_extra+=int(r1)
    #Random number following a binomial distribution.
    r2=semilla.binomial(n=n_B,p=pi_B_to_A*delta_t)
    n_A_extra+=int(r2)
    return n_A+n_A_extra-n_B_extra, n_B+n_B_extra-n_A_extra

#Theorical survial function
def survival_function(t,x_0,mu,sigma):
    #f=(np.exp(-(2*mu*x_0)/sigma**2)*(np.exp((2*mu*x_0)/sigma**2)*erf((math.sqrt(2)*x_0+math.sqrt(2)*mu*t)/(2*sigma*np.sqrt(t)))+erf((math.sqrt(2)*x_0-math.sqrt(2)*mu*t)/(2*sigma*np.sqrt(t)))+np.exp((2*mu*x_0)/sigma**2)-1))/2 
    f2=(1/2)*(1-np.exp(-(2*mu*x_0)/(sigma**2))*(1+erf((mu*t-x_0)/(math.sqrt(2)*sigma*np.sqrt(t))))+erf((mu*t+x_0)/(math.sqrt(2)*sigma*np.sqrt(t))))
    return f2

def tiempo_primera_pasada(t,x_0,mu,sigma,dominio):
    ext_t=1-survival_function(1000, x_0, mu, sigma)
    g=((x_0)/(sigma*math.sqrt(2*math.pi)*np.power(t,3/2))*np.exp(-(x_0+mu*t)**2/(2*sigma**2*t)))
    integral=integrate.trapz(g, dominio)
    return ext_t, g/integral*ext_t

#To compare the theorical first passage time with the simulations.
def prob_extincion(numero_simulaciones,entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Empty list to save time of extinction.
    tiempo_extincion_vector=[]
    for i in range(numero_simulaciones): 
        #In this case we are only interested in extinction.
        _, _, _, _, extincion,tiempo_extincion=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
        if extincion==True:
            #We save the time if there extinction.
            tiempo_extincion_vector+=[tiempo_extincion]
        print(i)
    #Probability of exinction.
    prob_extincion=len(tiempo_extincion_vector)/numero_simulaciones
    #In case of more than two extiction cases, we can make a histogram.
    plt.figure(figsize=(2.75,1.96),dpi=300)
    if len(tiempo_extincion_vector)>=2:
        intervalo=np.linspace(min(tiempo_extincion_vector),max(tiempo_extincion_vector))
        #plt.hist(x=tiempo_extincion_vector, bins=intervalo/2,color = "lightblue", ec="red", density='True')
        n, bins = np.histogram(tiempo_extincion_vector, bins=30, density=True)
        integral_2 = np.trapz(n, x=bins[:-1])
        n*=prob_extincion/integral_2
        integral_3 = np.trapz(n, x=bins[:-1])
        plt.bar(bins[:-1], n, width=np.diff(bins), color="lightblue", edgecolor="red", alpha=0.75)
        print("Integral (simulated pdf aproximation):", integral_3)
        #Collective model data.
        sigma=0.20
        mu=0.000072
        #Graphic representation of the first passage time.
        intervalo2=np.arange(0.01,max(tiempo_extincion_vector),0.05)
        (ext_t,primera_pasada)=tiempo_primera_pasada(intervalo2,math.log(10), mu, sigma, intervalo2)
        plt.plot(intervalo2,primera_pasada, label='Theorical first passage time')
        #Histogram.
        plt.title('Low growth and variance (PF)',fontsize = 10)
        plt.xlabel('Time extinction',fontsize = 10); plt.xticks(fontsize=8)
        plt.ylabel('PDF',fontsize = 10); plt.yticks(fontsize=8)
        plt.xlim(-3,200)
        plt.legend(fontsize=6)
        #Y mostramos.
        plt.show()
        #Checking if the first passage time integrates 1.
        integral=integrate.trapz(primera_pasada, intervalo2)
    return ext_t,prob_extincion, integral,tiempo_extincion_vector,intervalo2, primera_pasada,n,bins

def guardar_datos(datos, nombre_archivo):
    with open(nombre_archivo, 'wb') as archivo:
        pickle.dump(datos, archivo)

#Definition of the environment.
entornos=[0,1]
#Environment change rates.
kappa_0_to_1=1; kappa_1_to_0=1
#Number of bacteria of each phenotype
n_A=5; n_B=5.
#Growth rate in each environment.
k_A0=2; k_A1=-2; k_A=[k_A0,k_A1]; k_B0=0.2; k_B1=-0.2; k_B=[k_B0,k_B1]

#First passage time study
pi_A_to_B=6.894; pi_B_to_A=0.001
ext_t,prob_extincion, integral,tiempo_extincion_vector,intervalo2, primera_pasada,n,bins=prob_extincion(10000, entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A)
print('Simulated probability of extinction: ', prob_extincion)
print('The theorical probability of extinction is: ', ext_t)
print('Integral of the first passage time: ', integral)
#Save data 
guardar_datos((ext_t,prob_extincion, integral,tiempo_extincion_vector,intervalo2, primera_pasada,n,bins), "primera_pasada_lg.pkl")
