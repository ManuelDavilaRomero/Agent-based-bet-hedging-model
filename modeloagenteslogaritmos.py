# -*- coding: utf-8 -*-
"""
@author: Manuel
"""

#Modelo de agentes centrado en el logaritmo de la población.

#Importamos las librerías que vamos a necesitar.
from sklearn import preprocessing
from scipy.special import erf 
from scipy.stats import norm
import random, math 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

#Esta es la función principal.
def fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Para que la división temporal sea lo suficientemente pequeña, debemos considerar al menos dos órdenes de magnitud 
    #respecto a los rates característicos.
    delta_t=abs(10**(-2)*min(1/kappa_0_to_1,1/kappa_1_to_0,1/k_A[0],1/k_A[1],1/k_B[0],1/k_B[1],1/pi_A_to_B,1/pi_B_to_A))
    #El tiempo define el número de iteraciones.
    tiempo=np.arange(delta_t,50,delta_t)
    #Población inicial para la representación. En este caso nos interesa el logaritmo de la población.
    n_Alista=[math.log(n_A)]
    n_Blista=[math.log(n_B)]
    #En las variables totales no usamos el logaritmo de momento.
    n_total=[(n_A+n_B)]
    n_total_2=[((n_A+n_B)**2)]
    #Generemos la semilla para los números aleatorios.
    semilla=np.random.default_rng()
    #Elegimos un entorno inicial al azar.
    entorno_actual=semilla.choice([0,1],p=[0.5,0.5])
    #Lo guardamos en una lista.
    entornos=[entorno_actual]
    for j in range (len(tiempo)):
        #¿Cambia el ambiente?
        entorno_actual=cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0, delta_t)
        #Estudiamos la reproducción y muerte de las bacterias.
        n_B=reproduccion_muerte_bacteria(n_B, k_B[entorno_actual], delta_t, semilla)
        n_A=reproduccion_muerte_bacteria(n_A, k_A[entorno_actual], delta_t, semilla)
        #Recorremos todas las bacterias que han sobrevivido y vemos si cambian el fenotipo.
        n_A, n_B=cambio_fenotipo(n_A, n_B, pi_A_to_B, pi_B_to_A, delta_t, semilla)
        if n_A>0: #Aunque no estudiemos extinción, solo podemos tomar el logaritmo de cantidades postivas.
            n_Alista+=[math.log(n_A)]
        else: #Como usamos la binomial, si n_A<=0, vamos a acabar devolviendo solo 0.
            n_Alista+=[n_A]
        if n_B>0:
            n_Blista+=[math.log(n_B)]
        else:
            n_Blista+=[n_B]
        n_total+=[(n_A+n_B)]
        n_total_2+=[(n_A+n_B)**2]
        entornos+=[entorno_actual]
    #Convertimos en arrays para así poder operar elemento a elemento a sin problemas.
    n_Alista=np.array(n_Alista)
    n_Blista=np.array(n_Blista)
    n_total=np.array(n_total)
    n_total_2=np.array(n_total_2)
    return n_Alista, n_Blista,n_total, n_total_2, np.insert(tiempo,0,0), entornos #El insert simplemente es para añadir el origen de tiempos.
    #La variable n_total_2 no se ha usado finalmente pero se conserva porque puede ser interesante para calcular momentos de orden 2 en un futuro.
    
#Las siguientes tres funciones son exactamente iguales que en el programa que incluye simulaciones para la extinción.   
def cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0,delta_t):
    r=random.random() #Devuelve un número aleatorio entre 0 y 1.
    #Si estamos en el primer entorno y se cumple la condición...
    if entorno_actual==0 and r<=(kappa_0_to_1*delta_t):
        entorno_actual=1 #... pasamos al entorno 1
    #Lo mismo para el otro.
    elif entorno_actual==1 and r<=(kappa_1_to_0*delta_t):
        entorno_actual=0 
    return entorno_actual

#Recordemos que n es el número de bacterias, k es la tasa de crecimiento en el entorno.
def reproduccion_muerte_bacteria(bacterias,k,delta_t, semilla):
    #Generamos un número aleatorio que siga una distribución binomial, que nos indicará el número que actúan.
    r=semilla.binomial(n=bacterias,p=abs(k*delta_t))
    if k<0:
        #Se mueren ese número de individuos.
        return bacterias-int(r) #El int no hace falta pues la binomial devuelve un entero pero nos curamos en salud.
    if k>0:
        #Se reproducen ese número de individuos.
        return bacterias+int(r)
        
def cambio_fenotipo(n_A,n_B,pi_A_to_B,pi_B_to_A,delta_t, semilla):
    #Para aseguarnos de que al recorrer un fenotipo no modifiquemos el otro instaneamente, los cambios de fenotipos
    #se guardan en unas variables extra.
    n_A_extra=0; n_B_extra=0
    #Generamos un número aleatorio siguiendo una distribución binomial. 
    r1=semilla.binomial(n=n_A,p=pi_A_to_B*delta_t)
    #Ese número nos da el número de bacterias del tipo A que han cambiado al B.
    n_B_extra+=int(r1)
    #Generamos un número aleatorio siguiendo una distribución binomial. 
    r2=semilla.binomial(n=n_B,p=pi_B_to_A*delta_t)
    #Este otro nos indica el número del tipo B que han cambiado al A.
    n_A_extra+=int(r2)
    #El número final de cada fenotipo, es el original más los que han cambiado del fenotipo contrario y los que se 
    #pierden porque han cambiado a su opuesto.
    return n_A+n_A_extra-n_B_extra, n_B+n_B_extra-n_A_extra

#La función de análisis del logaritmo de la población nos permite estudiar en distintas posiciones temporales su distribución.
def analisis_log(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A, posicion_analisis_1,pos_ana_2, pos_ana_3):
    #Creamos listas vacías para analizar la posición que nos interese:
    n_analizar_1=[]
    n_analizar_2=[]
    n_analizar_3=[]
    #Hacemos diez mil simulaciones:
    for i in range(10000):
        n, n_2,tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:5]
        #Como vamos a estudiar el logaritmo de la población no nos interesa los valores iguales a 0.
        if n[posicion_analisis_1]>0:
            n_analizar_1+=[math.log(n[posicion_analisis_1])]
        if n[pos_ana_2]>0:
            n_analizar_2+=[math.log(n[pos_ana_2])]
        if n[pos_ana_3]>0:
            n_analizar_3+=[math.log(n[pos_ana_3])]
    #El por qué no incluimos directamente el logaritmo en el programa principal únicamente se debe a que conviene 
    #guardar los 0 de la población, que no tienen sentido al tomar logaritmos.
    #Numpy nos permite obtener la media y varianza y representar una gaussiana con esas características.
    t1=round(tiempo[posicion_analisis_1],1)
    t2=round(tiempo[pos_ana_2],1)
    t3=round(tiempo[pos_ana_3],1)
    std=np.std(n_analizar_1,ddof=1)
    mean=np.mean(n_analizar_1)
    dominio=np.linspace(np.min(n_analizar_1),np.max(n_analizar_1))
    plt.plot(dominio, norm.pdf(dominio,mean,std))
    plt.hist(n_analizar_1,color = "lightblue", bins=30,ec='red',density=True, stacked=True)
    plt.title('Histograma del ln $n(t)$ a tiempo '+str(t1)+' con $π_{A\_to\_B}=$'+str(pi_A_to_B)+',$π_{B\_to\_A}=$'+str(pi_B_to_A), fontsize=13)
    legend_elements_1 = [Line2D([0], [0], color='b', label='Gaussiana con $\mu=%.3f,\ \sigma=%.3f$' %(mean, std), markersize=15)]
    plt.legend(handles=legend_elements_1, loc='best')
    plt.xlabel('Ln n(t)', fontsize=13); plt.xticks(fontsize=14) 
    plt.ylabel('Frecuencia', fontsize=13);  plt.yticks(fontsize=13)
    plt.show()
    std=np.std(n_analizar_2,ddof=1)
    mean=np.mean(n_analizar_2)
    dominio=np.linspace(np.min(n_analizar_2),np.max(n_analizar_2))
    plt.plot(dominio, norm.pdf(dominio,mean,std))
    plt.hist(n_analizar_2,color = "lightblue", bins=30, ec='red',density=True, stacked=True)
    plt.title('Histograma del ln $n(t)$ a tiempo '+str(t2)+' con $π_{A\_to\_B}=$'+str(pi_A_to_B)+',$π_{B\_to\_A}=$'+str(pi_B_to_A), fontsize=13)
    legend_elements_2 = [Line2D([0], [0], color='b', label='Gaussiana con $\mu=%.3f,\ \sigma=%.3f$' %(mean, std), markersize=15)]
    plt.legend(handles=legend_elements_2, loc='best')
    plt.xlabel('Ln n(t)', fontsize=13); plt.xticks(fontsize=14) 
    plt.ylabel('Frecuencia', fontsize=13);  plt.yticks(fontsize=13)
    plt.show()
    std=np.std(n_analizar_3,ddof=1)
    mean=np.mean(n_analizar_3)
    dominio=np.linspace(np.min(n_analizar_3),np.max(n_analizar_3))
    plt.plot(dominio, norm.pdf(dominio,mean,std))
    plt.hist(n_analizar_3, color = "lightblue", bins=30, ec='red',density=True, stacked=True)
    plt.title('Histograma del ln $n(t)$ a tiempo '+str(t3)+' con $π_{A\_to\_B}=$'+str(pi_A_to_B)+',$π_{B\_to\_A}=$'+str(pi_B_to_A), fontsize=13)
    legend_elements_3 = [Line2D([0], [0], color='b', label='Gaussiana con $\mu=%.3f,\ \sigma=%.3f$' %(mean, std), markersize=15)]
    plt.legend(handles=legend_elements_3, loc='best')
    plt.xlabel('Ln n(t)', fontsize=13); plt.xticks(fontsize=14) 
    plt.ylabel('Frecuencia', fontsize=13);  plt.yticks(fontsize=13)
    plt.show()
    return tiempo[posicion_analisis_1]

#Esta función es para estudiar las trayectorias de los logaritmos de la población.
def trayectorias_log(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Generamos una simulación.
    n_A_1,n_B_1,n, n_2, tiempo, entornos=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
    #Representamos el logartimo de n_A frente al de n_B usando colores según el entorno.
    figure, axis = plt.subplots(2, 2,figsize=(10, 7))
    for j in range(len(n_B_1)):
        if entornos[j]==0:
            axis[0,0].plot(n_A_1[j],n_B_1[j],'.', color='blue')
        else:
            axis[0,0].plot(n_A_1[j],n_B_1[j],'.', color='red')
    axis[0,0].set_title('Trayectorias de los logaritmos de la población'); axis[0,0].set(xlabel='Ln $n_A(t)$', ylabel='Ln $n_B(t)$')
    legend_elements_1 = [Line2D([0], [0], marker='.', color='w', label='Entorno 0',markerfacecolor='b', markersize=15),
                   Line2D([0], [0], marker='.', color='w', label='Entorno 1',markerfacecolor='r', markersize=15),
                   Line2D([0], [0], color='g', label='Línea $\phi_{0}$', markersize=15),
                   Line2D([0], [0], color='k', label='Línea $\phi_{1}$', markersize=15)]
    #Nos faltan las líneas de equilibrio.
    Sigma_0=k_A[0]-k_B[0]
    Sigma_1=k_A[1]-k_B[1]
    pi_total=pi_A_to_B+pi_B_to_A
    phi_1=(Sigma_0-pi_total+math.sqrt((Sigma_0-pi_total)**2+4*Sigma_0*pi_B_to_A))/(2*Sigma_0)
    phi_2=(Sigma_1-pi_total+math.sqrt((Sigma_1-pi_total)**2+4*Sigma_1*pi_B_to_A))/(2*Sigma_1)
    x=np.linspace(min(n_A_1),max(n_A_1),10000)
    axis[0,0].plot(x,x+math.log((1/phi_1)-1),'g')
    axis[0,0].plot(x,x+math.log((1/phi_2)-1),'k')
    axis[0,0].legend(handles=legend_elements_1, loc='best')
    #Y también representamos la evolución temporaral del logaritmo de n_B(t).
    for j in range(len(n_B_1)):
        print(entornos[j])
        if entornos[j]==0:
            axis[0,1].plot(tiempo[j],n_B_1[j],'.', color='blue')
        else:
            axis[0,1].plot(tiempo[j],n_B_1[j],'.', color='red')
    axis[0,1].set_title('Trayectoria del ln $n_B(t)$'); axis[0,1].set(xlabel='Tiempo',ylabel='Ln $n_B(t)$')
    legend_elements_2 = [Line2D([0], [0], marker='.', color='w', label='Entorno 0',markerfacecolor='b', markersize=15),
                   Line2D([0], [0], marker='.', color='w', label='Entorno 1',markerfacecolor='r', markersize=15)]
    axis[0,1].legend(handles=legend_elements_2, loc='best')
    axis[0,0].get_shared_y_axes().join(axis[0,0], axis[0,1])
    #Represtnamos también la evolución de ln n_A.
    for j in range(len(n_A_1)):
        if entornos[j]==0:
            axis[1,0].plot(tiempo[j],n_A_1[j],'.', color='blue')
        else:
            axis[1,0].plot(tiempo[j],n_A_1[j],'.', color='red')
    axis[1,0].set(ylabel='Ln $n_A(t)$'); axis[1,0].set(xlabel='Tiempo'); axis[1,0].set_title('Trayectoria del ln $n_A(t)$')
    axis[1,0].legend(handles=legend_elements_2, loc='best')
    #Y la evoluación de la población total.
    for j in range(len(n_B_1)):
        if entornos[j]==0:
            axis[1,1].plot(tiempo[j],(np.exp(n_B_1[j])+np.exp(n_A_1[j])),'.', color='blue') #La exponencial es para eliminar logaritmos.
        else:
            axis[1,1].plot(tiempo[j],(np.exp(n_B_1[j])+np.exp(n_A_1[j])),'.', color='red')
    axis[1,1].set_title('Trayectoria  de $n_B(t)+n_A(t)$'); axis[1,1].set(xlabel='Tiempo',ylabel='$n_B(t)+n_A(t)$')
    axis[1,1].legend(handles=legend_elements_2, loc='best')
    axis[1,1].set(yscale='log')
    figure.tight_layout()
    plt.show()
    return 1 #Devolvemos 1 simplemente para que la función tenga un retorno.

#La siguiente función nos devuelve la función de supervivencia obtenida suponiendo el movimiento browiano con deriva y difusión.
def survival_function(t,x_0,mu,sigma):
    f=(np.exp(-(2*mu*x_0)/sigma**2)*(np.exp((2*mu*x_0)/sigma**2)*erf((math.sqrt(2)*x_0+math.sqrt(2)*mu*t)/(2*sigma*np.sqrt(t)))+erf((math.sqrt(2)*x_0-math.sqrt(2)*mu*t)/(2*sigma*np.sqrt(t)))+np.exp((2*mu*x_0)/sigma**2)-1))/2 
    return f


#Estas son las tasas y condiciones que tomamos de base en nuestro análisis.
#Definimos el entorno. El entorno 0 se corresponde con un ambiente normal mientras que el 1 uno con antibióticos.
entornos=[0,1]
#La tasa de cambio entre entornos viene dada por las kappas; kappa_0_to_1 para cambiar del entorno 0 al 1 y kappa_1_to_0
#para el proceso inverso.
kappa_0_to_1=1; kappa_1_to_0=1
#El número de bacterias las contabilizamos en las siguientes variables. Los fenotipos son A y B.
#Iniciamos la simulación en un total de 10 bacterias.
n_A=5; n_B=5
#Necesitamos también las tasas de crecimiento de cada fenotipo en cada ambiente.
#Por comodidad, se guardan las de cada fenotipo en una lista. 
k_A0=2; k_A1=-2; k_A=[k_A0,k_A1]; k_B0=0.2; k_B1=-0.2; k_B=[k_B0,k_B1]

'''#Análisis de la evolución del logaritmo de la población.
pi_A_to_B=6.894; pi_B_to_A=0.001
_=analisis_log(entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A,250,500,750)'''

#Cálculo de la probabilidad de extinción teórica.

'''#Comenzamos definiendo el intervalo de tiempos. También podemos calcular directamente la probabilidad a un tiempo t dado
#introduciéndolo como argumento de la función. Se elige tiempo 1000 para comparar con la extinción de la simulación.
t=np.linspace(0,1000,20000)
#Necesitamos la condición inicial, que el logartimo de la población inicial así como la media y la varianza que describen 
#el par de valores \pi y se obtienen del modelo colectivo.
x_0=math.log(10)
mu=0.210
sigma=0.91
#Calculamos la función de supervivencia y la representamos gráficamente.
supervivencia=survival_function(t,x_0,mu,sigma)
plt.plot(t,supervivencia,label='Función de supervivencia')
#La función de extinción es simplemente el complementario de la función de supervivencia.
extincion=1-supervivencia 
plt.plot(t,extincion,label='Función de extinción' )
#Personalizamos la figura.
plt.xlabel('Tiempo', fontsize=14); plt.yticks(fontsize=14); plt.xticks(fontsize=14)
plt.xlim(0,100)
plt.legend(fontsize=14)
plt.grid()
plt.show()
#Finalmente, mostramos en pantalla el último valor que toma.
extincion=1-supervivencia[-1]
print('A tiempo ',t[-1],'la probabilidad de extinción es ',extincion,'.')
#También puede ser interesante estudiar cómo varían estas funciones en función de la población inicial x_0.
#Creamos un vector de poblaciones iniciales.
x_0_2=np.arange(1,100000,1)
log_x_0=np.log(x_0_2)
#En este caso, se mantiene el tiempo fijo. 
t_2=1000
#Calculamos la función de supervivencia y la representamos gráficamente.
supervivencia_2=survival_function(t_2,log_x_0,mu,sigma)
plt.plot(log_x_0,supervivencia_2, label='Función de supervivencia')
#La función de extinción es simplemente el complementario de la función de supervivencia.
extincion_2=1-supervivencia_2
plt.plot(log_x_0,extincion_2, label='Función de extinción' )
plt.xlabel('Población inicial ln $n(t=0)$', fontsize=14)
plt.yticks(fontsize=14); plt.xticks(fontsize=14)
plt.legend(fontsize=14)
plt.grid()
plt.show()'''

#Nota: Aunque se muestre un error de división por 0 para t=0, como el denominador que se anula se encuentra dentro 
#de la función error, que tiende a 1 en infinito, la libraría nos devuelve el valor 1 en este caso.'''

#Trayectorias.

'''#Para estudiar las trayectorias de los logaritmos simplmente tenemos que indicar los valores de \pi deseados.
pi_A_to_B=0.569; pi_B_to_A=0.207
_=trayectorias_log(entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A)
#Observamos que si llegamos a una población nula dibujaremos un punto el origen repetidas veces.
#Puede ser necesario ejecutar el programa varias veces para obtener resultados que podamos interpretar debido al factor estocástico.'''

