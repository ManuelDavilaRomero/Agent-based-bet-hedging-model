# -*- coding: utf-8 -*-
"""
@author: Manuel
"""

#Agent-based model focus on the logarithm of the population.

from scipy.special import erf 
from scipy.stats import norm
import random, math 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.integrate import quad 
import pickle

#This is the main function.
def fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Definition of delta t.
    delta_t=abs(1000**(-1)*min(1/kappa_0_to_1,1/kappa_1_to_0,1/k_A[0],1/k_A[1],1/k_B[0],1/k_B[1],1/pi_A_to_B,1/pi_B_to_A))
    #0.005
    #Time defines the number of iterations.
    tiempo=np.arange(delta_t,50,delta_t)
    #Initial population.
    n_Alista=np.zeros(len(tiempo))
    n_Alista[0]=math.log(n_A)
    n_Blista=np.zeros(len(tiempo))
    n_Blista[0]=math.log(n_B)
    n_total=np.zeros(len(tiempo))
    n_total[0]=(n_A+n_B)
    #Seed for random generator (123456) for reproducibility if only 1 simulation is needed.
    semilla=np.random.default_rng()
    #We choose an initial environment randomly.
    entornos=np.zeros(len(tiempo))
    entorno_actual=semilla.choice([0,1],p=[0.5,0.5])
    entornos[0]=entorno_actual
    for j in range (1,len(tiempo)):
        #Does the environmet change?
        entorno_actual=cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0, delta_t)
        #Death and reproduction of bacteria.
        n_B=reproduccion_muerte_bacteria(n_B, k_B[entorno_actual], delta_t, semilla)
        n_A=reproduccion_muerte_bacteria(n_A, k_A[entorno_actual], delta_t, semilla)
        #Change of phenotype.
        n_A, n_B=cambio_fenotipo(n_A, n_B, pi_A_to_B, pi_B_to_A, delta_t, semilla)
        if n_A+n_B>=1:    
            n_total[j]=(n_A+n_B)
            if n_A>=1:
                n_Alista[j]=math.log(n_A)
            if n_B>1:
                n_Blista[j]=math.log(n_B)
            entornos[j]=entorno_actual
        else:
            #Extinction
            break
    return n_Alista, n_Blista,n_total, tiempo, entornos 
    
def cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0,delta_t):
    #r=semilla.random() 
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

def normal_pdf(x, mu, sigma):
    return norm.pdf(x, loc=mu, scale=sigma)

def analisis_log(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A, posicion_analisis_1):
    #Empty list to study several positions at once.
    n_analizar_1=[]
    for i in range(10000):
        n,tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
        if n[posicion_analisis_1]>0:
            n_analizar_1+=[math.log(n[posicion_analisis_1])]
            print(i)
    t1=round(tiempo[posicion_analisis_1],1)
    
    #Survival with the logarithim
    supervivencia=len(n_analizar_1)/10000
    print(supervivencia)
    mu_fit=np.mean(n_analizar_1)
    sigma_fit=np.std(n_analizar_1) 
    x_range=np.linspace(min(n_analizar_1), max(n_analizar_1), len(n_analizar_1))
    y_fit=normal_pdf(x_range, mu_fit, sigma_fit)   
    area_bajo_curva=np.trapz(y_fit, x_range)
    factor_normalizacion=1/area_bajo_curva
    y_fit_normalizada=y_fit*factor_normalizacion*supervivencia
    integral=np.trapz(y_fit_normalizada,x_range)
    print('Truncated gaussian integral:', integral)
    plt.plot(x_range, y_fit_normalizada, 'r-')
    
    #Plotting the pdf aproximation
    n, bins = np.histogram(n_analizar_1,bins=25,density=True)
    integral_2 = np.trapz(n, x=bins[:-1])
    n*=supervivencia/integral_2
    integral_3 = np.trapz(n, x=bins[:-1])
    print('PDF aproximation integral:', integral_3)
    plt.bar(bins[:-1], n, width=np.diff(bins), color="lightblue", edgecolor="red", alpha=0.75)
    
    #Fokker-planck resolution
    x_range_2=np.linspace(0, max(n_analizar_1), len(n_analizar_1))
    fokker_planck=density_probability(x_range_2,math.log(10),0.000072,np.sqrt(0.02*2),t1)
    super_nor=survival_function(t1,math.log(10),0.000072,np.sqrt(0.02*2)) 
    print('Theorical survival: ', super_nor)
    #area_bajo_curva_2=np.trapz(fokker_planck, x_range_2)
    plt.plot(x_range_2,fokker_planck,'g')
    vector_equiespaciado=np.linspace(-5,40,100000)
    no_extincion=gaussiana(0.000072,0.2/2, t1, vector_equiespaciado, np.log(10))
    plt.plot(vector_equiespaciado,no_extincion,'--',color='g')
    
    plt.title('Low growth and variance at $t=$'+str(t1), fontsize=14)
    #legend_elements_1 = [Line2D([0], [0], color='b', label='$\mu=%.2f,\ \sigma=%.2f$' %(mean, std), markersize=15)]
    legend_elements_1 = [Line2D([0], [0], linestyle='--', color='g', label='Collective without extinction', markersize=15),
                        Line2D([0], [0], color='g', label='Collective model', markersize=15),
                         Line2D([0], [0], color='r', label='$\mu=%.2f,\ \sigma=%.2f$' %(mu_fit, sigma_fit), markersize=15)]
    plt.legend(handles=legend_elements_1, loc='best')
    plt.xlabel('Ln n(t)', fontsize=13); plt.xticks(fontsize=14) 
    plt.ylabel('Frequency', fontsize=13);  plt.yticks(fontsize=13)
    plt.show()
    
    return fokker_planck,y_fit_normalizada,n_analizar_1,x_range,x_range_2,supervivencia,super_nor,t1,mu_fit,sigma_fit,n, bins


def gaussiana(lambdad, D, t, x, x_0):
    f = lambda x_var: (1/np.sqrt(4*np.pi*D*t))*np.exp(-(x_var-x_0-lambdad*t)**2/(4*D*t))
    integral, _ = quad(f, -np.inf, np.inf)
    return f(x) / integral

def media_varianza(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Empty matrix to save all cases. Functional but not optimal
    n10,tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    n40,_=fluacting_environments(entornos,20,20,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    matriz10=np.zeros((1000, len(tiempo)))
    matriz40=np.zeros((1000, len(tiempo)))
    matriz10[0,:]=np.log(n10/10)
    matriz40[0,:]=np.log(n40/40)
    for i in range(1,1000):
        n10=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:3]
        n40=fluacting_environments(entornos,20,20,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:3]
        matriz10[i,:]=np.log(np.array(n10)/10)
        matriz40[i,:]=np.log(np.array(n40)/40)
        print(i)
    #Transform -inf to nan to not count them
    matriz10[matriz10 == -np.inf] = np.nan
    matriz40[matriz40 == -np.inf] = np.nan
    matriz10 = matriz10/tiempo
    media_columnas_10= np.nanmean(matriz10, axis=0)
    varianza_columnas_10 = (np.nanvar(matriz10, axis=0)*tiempo)/2
    matriz40 = matriz40/tiempo
    media_columnas_40= np.nanmean(matriz40, axis=0)
    varianza_columnas_40 = (np.nanvar(matriz40, axis=0)*tiempo)/2
    
    proporcion_nan_10 = np.sum(np.isnan(matriz10[:, -1])) / matriz10.shape[0]
    proporcion_nan_40 = np.sum(np.isnan(matriz40[:, -1])) / matriz40.shape[0]
    print("Extinction 10:", proporcion_nan_10)
    print("Extinction 40:", proporcion_nan_40)

    #Figures. Change the label and values for each case
    '''
    plt.figure(dpi=600)
    plt.xlabel('Time', fontsize=13); plt.xticks(fontsize=14) 
    plt.ylabel('$\Lambda_t$ (agent)', fontsize=11);  plt.yticks(fontsize=13)
    legend_elements = [
    Line2D([0], [0], color='brown', linestyle='-', label='Agent ($n_0=10$)'),
    Line2D([0], [0], color='brown', linestyle='--', label='Agent ($n_0=40$)'),
    Line2D([0], [0], color='brown', linestyle='dotted', label='$D$ (collective)')]
    # Agregar título a la leyenda
    plt.legend(handles=legend_elements, loc='best', title='Low growth and variance (PF)')

    plt.plot(tiempo[10:-1],media_columnas_10[10:-1], color='brown')
    plt.plot(tiempo[10:-1],media_columnas_40[10:-1], color='brown',linestyle='--')
    plt.axhline(y=0.000072, color='brown', linestyle='dotted')
    plt.ylim(-0.01, 0.04)
    plt.show()
    
    plt.figure(dpi=600)
    plt.plot(tiempo[10:-1],varianza_columnas_10[10:-1], color='brown')
    plt.plot(tiempo[10:-1],varianza_columnas_40[10:-1], color='brown',linestyle='--')
    plt.xlabel('Time', fontsize=13); plt.xticks(fontsize=14) 
    plt.ylabel('$Var(\Lambda_t)t/2$ (agent)', fontsize=10);  plt.yticks(fontsize=13)
    legend_elements = [
    Line2D([0], [0], color='brown', linestyle='-', label='Agent ($n_0=10$)'),
    Line2D([0], [0], color='brown', linestyle='--', label='Agent ($n_0=40$)'),
    Line2D([0], [0], color='brown', linestyle='dotted', label='$D$ (collective)')]
    # Agregar título a la leyenda
    plt.legend(handles=legend_elements, loc='best', title='Low growth and variance (PF)')
    plt.axhline(y=0.020, color='brown', linestyle='dotted')
    plt.ylim(-0.01, 0.08)
    plt.show()
    '''
    
    '''
    plt.figure(dpi=600)
    plt.xlabel('Time', fontsize=13); plt.xticks(fontsize=14) 
    plt.ylabel('$\Lambda_t$ (agent)', fontsize=11);  plt.yticks(fontsize=13)
    legend_elements = [
    Line2D([0], [0], color='green', linestyle='-', label='Agent ($n_0=10$)'),
    Line2D([0], [0], color='green', linestyle='--', label='Agent ($n_0=40$)'),
    Line2D([0], [0], color='green', linestyle='dotted', label='$\Lambda$ (collective)')]
    # Agregar título a la leyenda
    plt.legend(handles=legend_elements, loc='best', title='Maximum growth (PF)')

    plt.plot(tiempo,media_columnas_10, color='green')
    plt.plot(tiempo,media_columnas_40, color='green',linestyle='--')
    plt.axhline(y=0.239, color='green', linestyle='dotted')
    plt.ylim(-0.01,0.5)
    plt.show()
    
    plt.figure(dpi=600)
    plt.plot(tiempo,varianza_columnas_10, color='green')
    plt.plot(tiempo,varianza_columnas_40, color='green',linestyle='--')
    plt.xlabel('Time', fontsize=13); plt.xticks(fontsize=14) 
    plt.ylabel('$Var(\Lambda_t)t/2$ (agent)', fontsize=10);  plt.yticks(fontsize=13)
    legend_elements = [
    Line2D([0], [0], color='green', linestyle='-', label='Agent ($n_0=10$)'),
    Line2D([0], [0], color='green', linestyle='--', label='Agent ($n_0=40$)'),
    Line2D([0], [0], color='green', linestyle='dotted', label='$D$ (collective)')]
    # Agregar título a la leyenda
    plt.legend(handles=legend_elements, loc='best', title='Maximum growth (PF)')
    plt.axhline(y=0.638, color='green', linestyle='dotted')
    plt.ylim(0.25,0.7)
    plt.show()
    '''
    return matriz10, media_columnas_10, varianza_columnas_10,matriz40, media_columnas_40, varianza_columnas_40

def guardar_datos(datos, nombre_archivo):
    with open(nombre_archivo, 'wb') as archivo:
        pickle.dump(datos, archivo)

#Studying the logarithm of the population.
def trayectorias_log(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #We do a simulation.
    n_A_1,n_B_1,n,tiempo,entornos=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
    #We represente the logarithm of n_A versus n_B using colours depeding on the environment.
    figure, axis = plt.subplots(2, 2,figsize=(8.5, 6),dpi=300)
    figure.tight_layout(pad=2.0)
    #axis[0,0].set_title('Trajectories of the logarithm of the population'); 
    axis[0,0].set_xlabel('Ln $n_A(t)$', fontsize=12)  
    axis[0,0].set_ylabel('Ln $n_B(t)$', fontsize=12)
    axis[0,0].tick_params(labelsize=10)
    axis[0,0].text(0.92, 0.08, '(a)', transform=axis[0, 0].transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='right')
    legend_elements_1 = [Line2D([0], [0], marker='.', color='w', label='Environment 0',markerfacecolor='b', markersize=10),
                   Line2D([0], [0], marker='.', color='w', label='Environment 1',markerfacecolor='r', markersize=10),
                   Line2D([0], [0], color='g', label='$\phi_{0}$', markersize=10),
                   Line2D([0], [0], color='k', label='$\phi_{1}$', markersize=10)]
    #Equilibrium.
    Sigma_0=k_A[0]-k_B[0]
    Sigma_1=k_A[1]-k_B[1]
    pi_total=pi_A_to_B+pi_B_to_A
    phi_1=(Sigma_0-pi_total+math.sqrt((Sigma_0-pi_total)**2+4*Sigma_0*pi_B_to_A))/(2*Sigma_0)
    phi_2=(Sigma_1-pi_total+math.sqrt((Sigma_1-pi_total)**2+4*Sigma_1*pi_B_to_A))/(2*Sigma_1)
    x=np.linspace(min(n_A_1),max(n_A_1),10000)
    axis[0,0].plot(x,x+math.log((1/phi_1)-1),'g', linewidth=1)
    axis[0,0].plot(x,x+math.log((1/phi_2)-1),'k', linewidth=1)
    for j in range(len(n_B_1)):
        if entornos[j]==0:
            axis[0,0].plot(n_A_1[j],n_B_1[j],'.', color='blue',markersize=4)
        else:
            axis[0,0].plot(n_A_1[j],n_B_1[j],'.', color='red',markersize=4)
    axis[0,0].legend(handles=legend_elements_1, loc='best')
    #axis[0,1].set_title('Trajectory of ln $n_B(t)$'); 
    axis[0,1].set_xlabel('Time',fontsize=12)
    axis[0,1].set_ylabel('Ln $n_B(t)$',fontsize=12)
    axis[0,1].tick_params(labelsize=10)
    axis[0,1].text(0.92, 0.08, '(b)', transform=axis[0, 1].transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='right')
    legend_elements_2 = [Line2D([0], [0], marker='.', color='w', label='Ambiente 0',markerfacecolor='b', markersize=15),
                   Line2D([0], [0], marker='.', color='w', label='Ambiente 1',markerfacecolor='r', markersize=15)]
    #axis[0,1].legend(handles=legend_elements_2, loc='best')
    axis[0,1].plot(tiempo,n_A_1+math.log((1/phi_1)-1),'g', linewidth=1)
    axis[0,1].plot(tiempo,n_A_1+math.log((1/phi_2)-1),'k', linewidth=1)
    #Temporal evolution of ln n_B.
    for j in range(len(n_B_1)):
        print(entornos[j])
        if entornos[j]==0:
            axis[0,1].plot(tiempo[j],n_B_1[j],'.', color='blue',markersize=4)
        else:
            axis[0,1].plot(tiempo[j],n_B_1[j],'.', color='red',markersize=4)
    axis[0,0].get_shared_y_axes().join(axis[0,0], axis[0,1])
    axis[1,0].set_ylabel('Ln $n_A(t)$',fontsize=12); axis[1,0].set_xlabel('Time',fontsize=12); 
    axis[1,0].tick_params(labelsize=10)
    #axis[1,0].set_title('Trajectory of ln $n_A(t)$')
    #axis[1,0].legend(handles=legend_elements_2, loc='best')
    axis[1,0].plot(tiempo,n_B_1-math.log((1/phi_1)-1),'g', linewidth=1)
    axis[1,0].plot(tiempo,n_B_1-math.log((1/phi_2)-1),'k', linewidth=1)
    axis[1,0].text(0.92, 0.08, '(c)', transform=axis[1, 0].transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='right')
    #Temporal evolution of ln n_A.
    for j in range(len(n_A_1)):
        if entornos[j]==0:
            axis[1,0].plot(tiempo[j],n_A_1[j],'.', color='blue',markersize=4)
        else:
            axis[1,0].plot(tiempo[j],n_A_1[j],'.', color='red',markersize=4)
    #Temporal evolution of total population.
    for j in range(len(n_B_1)):
        if entornos[j]==0:
            axis[1,1].plot(tiempo[j],(np.exp(n_B_1[j])+np.exp(n_A_1[j])),'.', color='blue',markersize=2) #La exponencial es para eliminar logaritmos.
        else:
            axis[1,1].plot(tiempo[j],(np.exp(n_B_1[j])+np.exp(n_A_1[j])),'.', color='red',markersize=2)
    #axis[1,1].set_title('Trajectory  of $n_B(t)+n_A(t)$'); 
    axis[1,1].set_xlabel('Time',fontsize=12) 
    axis[1,1].set_ylabel('$n_B(t)+n_A(t)$',fontsize=12)
    axis[1,1].tick_params(labelsize=10)
    axis[1,1].text(0.92, 0.08, '(d)', transform=axis[1, 1].transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='right')
    #axis[1,1].legend(handles=legend_elements_2, loc='best')
    axis[1,1].set(yscale='log')
    figure.tight_layout()
    plt.show()
    return n_A_1, n_B_1

#Theorical survial function
def survival_function(t,x_0,mu,sigma):
    f2=(1/2)*(1-np.exp(-(2*mu*x_0)/(sigma**2))*(1+erf((mu*t-x_0)/(math.sqrt(2)*sigma*np.sqrt(t))))+erf((mu*t+x_0)/(math.sqrt(2)*sigma*np.sqrt(t))))
    return f2

def ext_function_approx(t,x_0,mu,sigma):
    return math.sqrt(1/(2*np.pi))*(sigma/x_0)*np.exp(-(mu*t+x_0)**2/(2*sigma**2*t))*np.sqrt(t)

def ext_function_approx_long(mu,x_0,sigma):
    return np.exp(-(2*mu*x_0)/sigma**2)

def density_probability(x,x_0,mu,sigma,t):
   x=np.array(x)
   x_0=np.array(x_0)
   mu=np.array(mu)
   sigma=np.array(sigma)
   t=np.array(t)
   result = (1/np.sqrt(2*np.pi*sigma**2*t))*(-np.exp(-2*x_0*mu/sigma**2)*
              np.exp(-(x+x_0-mu*t)**2/(2*sigma**2*t))+
              np.exp(-(x_0-x+mu*t)**2 /(2*sigma**2 *t)))
   return result  
 
#Definition of the environment.
entornos=[0,1]
#Environment change rates.
kappa_0_to_1=1; kappa_1_to_0=1
#Number of bacteria of each phenotype.
n_A=5; n_B=5
#Growth rate in each environment.
k_A0=2; k_A1=-2; k_A=[k_A0,k_A1]; k_B0=0.2; k_B1=-0.2; k_B=[k_B0,k_B1]

#Analysis of the evolution of the logarithm of the population.
'''
pi_A_to_B=0.7; pi_B_to_A=0.1
fokker_planck,y_fit_normalizada,n_analizar_1,x_range,x_range_2,supervivencia,super_nor,t1,mu_fit,sigma_fit,n,bins=analisis_log(entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A,9000)
guardar_datos((fokker_planck, y_fit_normalizada, n_analizar_1, x_range, x_range_2, supervivencia, super_nor, t1, mu_fit, sigma_fit,n,bins), "datos_guardados_4.pkl")
'''
'''
matriz10, media_columnas_10, varianza_columnas_10,matriz40, media_columnas_40, varianza_columnas_40=media_varianza(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
'''
#Theoric extinction probability


'''
#Time intervale. We can also choose directly the time t as an argument of the function.
t=np.linspace(0,1000,20000)
#Initial condition (logarithm of the population). 
#Growth rate and variance for the collective model.
x_0=math.log(10)
mu=0.210
sigma=0.91
#Survival function and its representation.
supervivencia=survival_function(t,x_0,mu,sigma)
plt.figure(figsize=(2.75,1.96),dpi=700)
plt.plot(t,supervivencia,label='$F(t)$ survival function')
#The extinction funtion is the complementary.
extincion=1-supervivencia 
plt.plot(t,extincion,color='red',label='$E(t)$ extinction function' )
#Approximation
ext_1=ext_function_approx(t,x_0,mu,sigma)
ext_2=ext_function_approx_long(mu,x_0,sigma)*np.ones(len(t))
plt.plot(t[:800], ext_1[:800], ':', color='red')  # Para ext_1 de t=0 a t=35 (36 elementos)
plt.plot(t[200:], ext_2[200:], '--', color='red')  # Para ext_2 de t=10 hacia adelante
plt.plot(t[:800], 1 - ext_1[:800], ':', color='blue')  # Para ext_1 de t=0 a t=35 (36 elementos)
plt.plot(t[200:], 1 - ext_2[200:], '--', color='blue')  # Para ext_2 de t=10 hacia adelante
#Personalising the figure.
plt.xlabel('Time', fontsize=10); plt.yticks(fontsize=8); plt.xticks(fontsize=8)
plt.title('Minimun extinction PF',fontsize=10)
plt.xlim(0,50)
plt.legend(fontsize=7)
plt.grid()
plt.show()
#Showing the last value of the array.
extincion=1-supervivencia[-1]
print('At time ',t[-1],'the probablity of extinction is ',extincion,'.')
#Studying the behavoir depending on x_0.
x_0_2=np.arange(1,100000,1)
log_x_0=np.log(x_0_2) 
t_2=1000
#Survival function.
plt.figure(figsize=(2.75,1.96),dpi=700)
supervivencia_2=survival_function(t_2,log_x_0,mu,sigma)
plt.plot(log_x_0,supervivencia_2, label='$F(t)$ survival function')
#Extinction function.
extincion_2=1-supervivencia_2
plt.plot(log_x_0,extincion_2, label='$E(t)$ extinction function',color='red')
ext_1=ext_function_approx(t_2,x_0_2,mu,sigma)
ext_2=ext_function_approx_long(mu,x_0_2,sigma)
plt.plot(log_x_0,ext_2,'--',color='red')
plt.plot(log_x_0,1-ext_2,'--',color='blue')
plt.title('Minimun extinction PF',fontsize=10)
plt.xlabel('Initial population ln $n(t=0)$', fontsize=10)
plt.yticks(fontsize=8); plt.xticks(fontsize=8)
plt.legend(fontsize=7)
plt.grid()
plt.show()
'''
#Trajectories.
'''
#We only need to specify the phenotypic change rate.
pi_A_to_B=0.569; pi_B_to_A=0.207
n_A_1,n_B_1=trayectorias_log(entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A)
#Several executions may be need (in case of extinction we only plot a point several times).
'''