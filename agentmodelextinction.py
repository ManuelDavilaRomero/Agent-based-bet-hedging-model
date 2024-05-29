# -*- coding: utf-8 -*-
"""
@author: Manuel
"""

#Agent-based model, analysis of extinction.

#Libraries.
import random
import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib.lines import Line2D

#This is the main function.
def fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Definition of delta t.
    delta_t=100**(-1)*min(1/kappa_0_to_1,1/kappa_1_to_0,1/abs(k_A[0]),1/abs(k_A[1]),1/abs(k_B[0]),1/abs(k_B[1]),1/pi_A_to_B,1/pi_B_to_A)
    #Time defines the number of iterations.
    tiempo=np.arange(delta_t,1000,delta_t)
    #Initial population.
    n_total=np.zeros(len(tiempo))
    n_total[0]=(n_A+n_B)
    extincion=False #A positive number in the time origin: no extinction.
    tiempo_extincion=False #Nether there is an extinction time.
    tiempo_cota_superior=False #And the upper bound is still yet unsed.
    #Seed for random generator (123456) for reproducibility if only 1 simulation is needed
    semilla=np.random.default_rng()
    #We choose an initial environment randomly.
    entorno_actual=semilla.choice([0,1],p=[0.5,0.5])
    for j in range (1,len(tiempo)):
        #Does the environmet change?
        entorno_actual=cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0, delta_t)
        #Death and reproduction of bacteria.
        n_B=reproduccion_muerte_bacteria(n_B, k_B[entorno_actual], delta_t,semilla)
        n_A=reproduccion_muerte_bacteria(n_A, k_A[entorno_actual], delta_t,semilla)
        #Change of phenotype.
        n_A, n_B=cambio_fenotipo(n_A,n_B,pi_A_to_B,pi_B_to_A,delta_t,semilla)
        #We save the total number of bacteria to represent them.
        if (n_A+n_B)==0:
            extincion=True
            tiempo_extincion=tiempo[j]
            break
        n_total[j]=(n_A+n_B)
        if (n_A+n_B)>100000000: #Upper bound.
            tiempo_cota_superior=tiempo[j]
            for i in range(j, len(tiempo)):
                n_total[i]=(n_A+n_B)
            break
    return n_A, n_B,n_total,tiempo, extincion, tiempo_extincion, tiempo_cota_superior 
    
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
        #r bacteria die.
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

#Probability of extinction with active cases and upper bond cases.
def prob_extincion(numero_simulaciones,entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Empty list to save time of extinction or time of reaching the upper bond.
    tiempo_extincion_vector=[]
    tiempo_cota_superior_vector=[]
    #Aux variable to calculate the mean of the active populations.
    n_A_media=0
    n_B_media=0
    n_plot=[]
    plt.figure(figsize=(2.75,1.96),dpi=300)
    for i in range(numero_simulaciones): 
        n_Abis, n_Bbis, n_totalbis, tiempo, extincion,tiempo_extincion,tiempo_cota_superior=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
        if extincion==True:
            #If there is extinction, we save the time...
            tiempo_extincion_vector+=[tiempo_extincion]
        if extincion==False:
            #or the upper bond time...
            if max(n_totalbis)>100000000:
                tiempo_cota_superior_vector+=[tiempo_cota_superior]
            else:
                #... or the population.
                n_A_media+=n_Abis
                n_B_media+=n_Bbis
        n_plot.append(n_totalbis)
    print(tiempo)
    #Plotting 30 for better plot
    for i in range(30): 
        plt.plot(tiempo,n_plot[i],'-',markersize=0.1)
    #Personalized graphic.
    plt.xlabel('t', fontsize = 10); plt.ylabel('n(t)',fontsize = 10)
    plt.yticks(fontsize=8); plt.xticks(fontsize=8); plt.yscale('log')
    plt.ylim(top=10**7)
    plt.show()
    #Probability of extinction.
    prob_extincion=len(tiempo_extincion_vector)/numero_simulaciones
    #In case of more than two extiction cases, we can make a histogram.
    plt.figure(figsize=(2.75,1.96),dpi=300)
    if len(tiempo_extincion_vector)>=2:
        intervalo=np.linspace(min(tiempo_extincion_vector),max(tiempo_extincion_vector))
        plt.hist(x=tiempo_extincion_vector, bins=intervalo,color = "lightblue", ec="red")
        plt.xlabel('Extinction time',fontsize = 10); plt.xticks(fontsize=8)
        plt.ylabel('Frequency',fontsize = 10); plt.yticks(fontsize=8)
        plt.show()
    return prob_extincion, len(tiempo_extincion_vector), len(tiempo_cota_superior_vector), n_A_media,n_B_media

#If we only want the probablity of extinction without studying active populations or other aspects.
def prob_extincion_simple(numero_simulaciones,entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Counting the extinction cases.
    veces_extintas=0
    for i in range(numero_simulaciones):
        #We want the same population for every case. 
        _, _, _,_, extincion,_,_=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
        if extincion==True:
            veces_extintas+=1
    prob_extincion=veces_extintas/numero_simulaciones
    return prob_extincion

#Function for the Pareto front.
def frente_pareto(pi_A_to_B_vector):
    pi_B_to_A_vector=[]
    pi_A_to_B_vector_grafica=[]
    for i in range(len(pi_A_to_B_vector)):
        T=1/(2*3.50-1/pi_A_to_B_vector[i])
        #we are only intereseted in positive values and the region we are representing.
        if T>0 and T<1.5:
            pi_B_to_A_vector+=[T]
            pi_A_to_B_vector_grafica+=[pi_A_to_B_vector[i]]
    return pi_B_to_A_vector, pi_A_to_B_vector_grafica

#Definition of the environment.
entornos=[0,1]
#Environment change rates.
kappa_0_to_1=1; kappa_1_to_0=1
#Number of bacteria of each phenotype.
n_A=5; n_B=5
#Growth rate in each environment.
k_A0=2; k_A1=-2; k_A=[k_A0,k_A1]; k_B0=0.2; k_B1=-0.2; k_B=[k_B0,k_B1]

'''
#Trayectories.
#Case 1:
pi_A_to_B=0.263; pi_B_to_A=0.246
plt.figure(figsize=(2.75,1.96),dpi=300)
for _ in range(1):
    #Compute the trajetory...
    n_total1, tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    #and plot it.
    plt.plot(tiempo,n_total1,'.',color='green',markersize=0.3);
#Case 2:
pi_A_to_B=6.894; pi_B_to_A=0.001
for _ in range(1):
    n_total3, tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    plt.plot(tiempo,n_total3,'.',color='brown',markersize=0.1);
'''
'''
#Case 3:
pi_A_to_B=0.569; pi_B_to_A=0.207
for _ in range(100):
    n_total2, tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    plt.plot(tiempo,n_total2,'.',color='magenta',markersize=0.8)
#Case 4:
pi_A_to_B=0.74; pi_B_to_A=0.1
for _ in range(100):
    n_total4, tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    plt.plot(tiempo,n_total4,'.',color='orange',markersize=0.8);
'''
'''
#Personalize the plot.
#plt.title('Trajectories', fontsize=15)
plt.xlabel('t', fontsize=10); plt.ylabel('n(t)', fontsize=10); plt.xticks(fontsize=8); plt.yticks(fontsize=8)
custom_lines = [Line2D([0], [0], color='green', lw=4),
                #Line2D([0], [0], color='magenta', lw=4),
                #Line2D([0], [0], color='orange', lw=4),
                Line2D([0], [0], color='brown', lw=4)]
plt.legend(custom_lines, ['Maximum growth PF',
                          #'Minimal extinction in the PF',
                          #'Minimal extinction',
                          'Low growth and variance PF'],
           fontsize=6,
           bbox_to_anchor=(0.50, 1.05),
           loc='upper center',
           fancybox=True, shadow=True)
plt.yscale('log')
plt.ylim(top=10**7)
plt.show()
'''

#Simple extinction if we are interested in one point.
pi_A_to_B=6.894; pi_B_to_A=0.001
'''
numero_extincion=prob_extincion_simple(10000, entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A)
'''
'''
ext=[]
for i in range(10):
    print(i)
    ext.append(prob_extincion_simple(1000, entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A))
print('Extinction probability:', np.mean(ext))
print('with std ', np.std(ext)/np.sqrt(10))
'''


#Probability of extinction with active cases study.
'''
pi_A_to_B=0.263; pi_B_to_A=0.249
#We use the function created for this objective.
prob_extincion, numero_extincion, numero_cota_superior,n_A,n_B=prob_extincion(1000, entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A)
#Show the results
print('Probability of extinction: ', prob_extincion)
print('Number of cases with extinction: ', numero_extincion)
print('Upper bond cases: ', numero_cota_superior)
casos_finales=1000-int(numero_extincion)-int(numero_cota_superior) #Active cases: no extinction and no upper bond
print('Active cases: ',casos_finales)
#If there is active cases...
if casos_finales!=0:
    #...we can calculate an approximate stacionary case.
    print('Expected stacionary case: ', pi_B_to_A/pi_A_to_B)
    print('Reach stacionary case: ',(n_A/casos_finales)/(n_B/casos_finales))
'''


#Map of extinction.
'''
start_time = time.time()
#Vectors of the zone we want to study.
pi_A_to_B_vector=np.linspace(0.1,0.3,3)
pi_B_to_A_vector=np.linspace(0.5,1.5,11)
#Empty variables to save the point with least extinction and a matrix to save extinction.
extincion_minima=1
valores_extincion=[]
des_est=np.zeros((len(pi_B_to_A_vector),len(pi_A_to_B_vector)))
f=np.zeros((len(pi_B_to_A_vector),len(pi_A_to_B_vector)))
for i in range(len(pi_A_to_B_vector)):
    print(i)
    for j in range (len(pi_B_to_A_vector)):
        print(j)
        ext=[]
        for k in range(0,1000):
            ext.append(prob_extincion_simple(1000,entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B_vector[i],pi_B_to_A_vector[j]))
        f[j,i]=np.mean(ext)
        des_est[j,i]=np.std(ext)/np.sqrt(1000)
        print(f[j,i],'pi constants', pi_A_to_B_vector[i],pi_B_to_A_vector[j],'std',des_est[j,i])
        if f[j,i]<extincion_minima:
            extincion_minima=f[j,i]
            valores_extincion=[pi_A_to_B_vector[i],pi_B_to_A_vector[j]]    
print('The minimum extinction is ',extincion_minima)
print('For ',valores_extincion)
np.save('matriz_guardada_zoom2.npy', f)
np.savez('minima_extincion_zoom2.npy',valores_extincion)
np.save('des_est.npy',des_est)

#Vector for the Pareto front.
plt.figure(figsize=(2.75,1.96),dpi=300)
pi_A_to_B_vector_pareto=np.linspace(0.1,3,1000)
pi_B_to_A_vector_pareto, pi_A_to_B_vector_pareto=frente_pareto(pi_A_to_B_vector_pareto)
plt.plot(pi_A_to_B_vector_pareto,pi_B_to_A_vector_pareto,'-',color='blue',label='Pareto front',linewidth=2.5)
plt.xlim([0,3]); plt.ylim([0,1.5]); plt.legend(loc="upper right")
plt.plot(0.263, 0.246, marker="o", markersize=10, markeredgecolor="green", markerfacecolor="green")
plt.plot(0.569, 0.207, marker="o", markersize=10, markeredgecolor="magenta", markerfacecolor="magenta")
plt.legend(loc="upper right",fontsize=6)
plt.pcolormesh(pi_A_to_B_vector, pi_B_to_A_vector,f, cmap='OrRd',shading='nearest')
cbar = plt.colorbar()
cbar.set_label('Extinción', fontsize=10);
plt.xlabel(r'$π_{A\to B}$',fontsize=10); plt.ylabel(r'$π_{B\to A}$',fontsize=10); 
plt.xticks(fontsize=8); plt.yticks(fontsize=8)

e = int(time.time() - start_time)
print('{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))
'''
