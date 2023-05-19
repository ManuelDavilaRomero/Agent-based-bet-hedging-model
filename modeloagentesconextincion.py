# -*- coding: utf-8 -*-
"""
@author: Manuel
"""

#Modelo de agentes con análisis de la extinción.

#Importamos las librerías que vamos a necesitar.
import random
import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib.lines import Line2D

#Esta es la función principal.
def fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Para que la división temporal sea lo suficientemente pequeña, debemos considerar al menos dos órdenes de magnitud 
    #respecto a las tasas característicos.
    delta_t=abs(10**(-2)*min(1/kappa_0_to_1,1/kappa_1_to_0,1/k_A[0],1/k_A[1],1/k_B[0],1/k_B[1],1/pi_A_to_B,1/pi_B_to_A))
    #El tiempo define el número de iteraciones.
    tiempo=np.arange(delta_t,1000,delta_t)
    #Introducimos la población inicial.
    n_total=[(n_A+n_B)]
    extincion=False #Como introducimos un número positivo de bacterias en el origen temporal no hay extinción.
    tiempo_extincion=False #Tampoco habrá aún un tiempo de extinción.
    tiempo_cota_superior=False #Ni un tiempo donde se haya alcanzado la cota superior.
    #Generemos la semilla para los números aleatorios.
    semilla=np.random.default_rng()
    #Elegimos un entorno inicial al azar.
    entorno_actual=semilla.choice([0,1],p=[0.5,0.5])
    for j in range (len(tiempo)):
        #¿Cambia el ambiente?
        entorno_actual=cambio_ambiente(entorno_actual,kappa_0_to_1,kappa_1_to_0, delta_t)
        #Estudiamos la reproducción y muerte de las bacterias.
        n_B=reproduccion_muerte_bacteria(n_B, k_B[entorno_actual], delta_t,semilla)
        n_A=reproduccion_muerte_bacteria(n_A, k_A[entorno_actual], delta_t,semilla)
        #Recorremos todas las bacterias resultantes y vemos si cambian su fenotipo.
        n_A, n_B=cambio_fenotipo(n_A,n_B,pi_A_to_B,pi_B_to_A,delta_t,semilla)
        #Finalmente guardadamos el número total de individuos para luego representarlos.
        if (n_A+n_B)==0:
            extincion=True
            tiempo_extincion=tiempo[j]
            #Rellenamos el resto del vector y salimos del programa para ahorrar tiempo.
            for i in range(j, len(tiempo)):
                n_total+=[(0)]
            break
        n_total+=[(n_A+n_B)]
        if (n_A+n_B)>100000000: #Cota superior.
            tiempo_cota_superior=tiempo[j]
            #Rellenamos el resto del vector y salimos del programa para ahorrar tiempo.
            for i in range(j, len(tiempo)-1):
                n_total+=[(n_A+n_B)]
            break
    return n_A, n_B,n_total, np.insert(tiempo,0,0), extincion, tiempo_extincion, tiempo_cota_superior #El insert simplemente es para añadir el origen de tiempos.
    
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
    #Generamos un número aleatorio que siga una distribución binomial, que nos indicará el número de bacterias que interactuan.
    r=semilla.binomial(n=bacterias,p=abs(k*delta_t))
    if k<0:
        #Se mueren ese número de bacterias.
        return bacterias-int(r) #El int no hace falta pues la binomial devuelve un entero pero nos curamos en salud.
    if k>0:
        #Se reproducen ese número de seres.
        return bacterias+int(r)
        
def cambio_fenotipo(n_A,n_B,pi_A_to_B,pi_B_to_A,delta_t,semilla):
    #Para aseguarnos de que al recorrer un fenotipo no modifiquemos el otro instáneamente, los cambios de fenotipos
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

#La siguiente función calcula probabilidad de extinción pero además permite conocer la población activa así como
#los casos que alcanzan la cota superior.
def prob_extincion(numero_simulaciones,entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Creamos listas vacías para guardar los tiempos donde se produce la extinción o se alcanza la cota.
    tiempo_extincion_vector=[]
    tiempo_cota_superior_vector=[]
    #Creamos unas variables auxiliares para calcular la media de las poblaciones activas.
    n_A_media=0
    n_B_media=0
    for i in range(numero_simulaciones): #Bucle en el número de simulaciones.
        #Queremos que todos los intentos sean con la misma población por eso usamos las variables bis.
        n_Abis, n_Bbis, n_totalbis, tiempo, extincion,tiempo_extincion,tiempo_cota_superior=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
        if extincion==True:
            #Si hay extinción, guardemos el tiempo.
            tiempo_extincion_vector+=[tiempo_extincion]
        if extincion==False:
            #Si no hay, guardamos el tiempo de la cota...
            if max(n_totalbis)>100000000:
                tiempo_cota_superior_vector+=[tiempo_cota_superior]
            else:
                #... o la población de cada fenotipo para calcular la media.
                n_A_media+=n_Abis
                n_B_media+=n_Bbis
        #Representamos la trayectoria.
        plt.plot(tiempo,n_totalbis,'-',markersize=0.8)
    #Personalizamos la gráfica a nuestro gusto.
    plt.title('1000 simulaciones para las trayectorias',fontsize = 16)
    plt.xlabel('t', fontsize = 16); plt.ylabel('n(t)',fontsize = 14)
    plt.yticks(fontsize=14); plt.xticks(fontsize=14); plt.yscale('log')
    #Y la mostramos en pantalla.
    plt.show()
    #Calculamos la probabilidad de extinción.
    prob_extincion=len(tiempo_extincion_vector)/numero_simulaciones
    #En caso de que haya más de 1 un caso con extinción, tiene sentido representar un histograma.
    if len(tiempo_extincion_vector)>=2:
        intervalo=np.linspace(min(tiempo_extincion_vector),max(tiempo_extincion_vector))
        plt.hist(x=tiempo_extincion_vector, bins=intervalo,color = "lightblue", ec="red")
        #Personalizamos.
        plt.title('Histograma de tiempos de extinción',fontsize = 16)
        plt.xlabel('Tiempos de extinción',fontsize = 14); plt.xticks(fontsize=14)
        plt.ylabel('Frecuencia',fontsize = 14); plt.yticks(fontsize=14)
        #Y mostramos.
        plt.show()
    return prob_extincion, len(tiempo_extincion_vector), len(tiempo_cota_superior_vector), n_A_media,n_B_media

#Si solo nos interesa la probabilidad de extinción y no la cota superior.
def prob_extincion_simple(numero_simulaciones,entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A):
    #Creamos un contador para el número de extinciones.
    veces_extintas=0
    for i in range(numero_simulaciones):
        #Queremos que todos los intentos sean con la misma población. 
        n_Abis, n_Bbis, n_totalbis, tiempo, extincion,tiempo_extincion,tiempo_cota=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)
        if extincion==True:
            veces_extintas+=1
    prob_extincion=veces_extintas/numero_simulaciones
    return prob_extincion

#Función para representar el frente de Pareto.
def frente_pareto(pi_A_to_B_vector):
    pi_B_to_A_vector=[]
    pi_A_to_B_vector_grafica=[]
    for i in range(len(pi_A_to_B_vector)):
        T=1/(2*3.33-1/pi_A_to_B_vector[i])
        #Solo nos interesan los valores postivos y en el rango que vamos a representar.
        if T>0 and T<1.5:
            pi_B_to_A_vector+=[T]
            pi_A_to_B_vector_grafica+=[pi_A_to_B_vector[i]]
    return pi_B_to_A_vector, pi_A_to_B_vector_grafica

#Definimos el entorno. El entorno 0 se corresponde con un ambiente normal mientras que el 1 uno con antibióticos.
entornos=[0,1]
#La tasa de cambio entre entornos viene dada por las kappas; kappa_0_to_1 para cambiar del entorno 0 al 1 y kappa_1_to_0
#para el proceso inverso.
kappa_0_to_1=1; kappa_1_to_0=1
#El número de bacterias las contabilizamos en las siguientes variables:
n_A=5; n_B=5
#Necesitamos también las tasas de crecimiento de cada fenotipo en cada ambiente. Los fenotipos son A y B.
k_A0=2; k_A1=-2; k_A=[k_A0,k_A1]; k_B0=0.2; k_B1=-0.2; k_B=[k_B0,k_B1]

#Trayectorias.

'''#Caso 1:
#Mantendremos el resto de variables fijas y variaremos la probabilidad de cambio de fenotipo.
pi_A_to_B=0.569; pi_B_to_A=0.207
for _ in range(100):
    #Calculamos una trayectoria...
    n_total1, tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    #y la representamos.
    plt.plot(tiempo,n_total1,'.',color='blue',markersize=0.8);
#Caso 2:
pi_A_to_B=1; pi_B_to_A=0.1
for _ in range(100):
    n_total2, tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    plt.plot(tiempo,n_total2,'.',color='red',markersize=0.8);
#Caso 3:
pi_A_to_B=6; pi_B_to_A=0.001
for _ in range(100):
    n_total3, tiempo=fluacting_environments(entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B,pi_B_to_A)[2:4]
    plt.plot(tiempo,n_total3,'.',color='green',markersize=0.8);
#Personalizamos la gráfica.
plt.title('Trayectorias', fontsize=15)
plt.xlabel('t', fontsize=14); plt.ylabel('n(t)', fontsize=14); plt.xticks(fontsize=14); plt.yticks(fontsize=14)
custom_lines = [Line2D([0], [0], color='blue', lw=4),
                Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='green', lw=4)]
plt.legend(custom_lines, ['$π_{B\_to\_A}=0.569,π_{B\_to\_A}=0.207$', '$π_{B\_to\_A}=1,π_{B\_to\_A}=0.1$', '$π_{B\_to\_A}=6,π_{B\_to\_A}=0.001$'])
plt.yscale('log')
plt.show()'''

#Extinción simple. Si solo nos interesa un punto concreto, es lo más rápido.

'''#Definimos las constantes que nos faltan.
pi_A_to_B=6.894; pi_B_to_A=0.001
numero_extincion=prob_extincion_simple(10000, entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A)
print('La probabilidad de extinción es:', numero_extincion)'''


#Probabildidad de extinción con estudio de la cota superior y las poblaciones activas.

'''#Definimos las constantes que nos falta. 
pi_A_to_B=1; pi_B_to_A=0.5
#Usamos la función que hemos creado para ello. 
prob_extincion, numero_extincion, numero_cota_superior,n_A,n_B=prob_extincion(1000, entornos, n_A, n_B, kappa_0_to_1, kappa_1_to_0, k_A, k_B, pi_A_to_B, pi_B_to_A)
#Mostramos los resultados.
print('La probabilidad de extinción es:', prob_extincion)
print('El número casos donde hay extinción es:', numero_extincion)
print('El número de casos que llegan a la cota es:', numero_cota_superior)
casos_finales=1000-int(numero_extincion)-int(numero_cota_superior) #Las activas son las que no están ni extintas ni en la cota.
print('Los casos que llegan al tiempo final sin llegar a la cota:',casos_finales)
#Si hay casos activos...
if casos_finales!=0:
    #...podemos aproximar el estado estacionario.
    print('Estado estacionario esperado:', pi_B_to_A/pi_A_to_B)
    print('Estado estacionario alcanzado:',(n_A/casos_finales)/(n_B/casos_finales))'''




#Mapa de extincion.
'''start_time = time.time()
#Creamos los vectores de la zona que vamos a estudiar.
pi_A_to_B_vector=np.linspace(0.1,3,10)
pi_B_to_A_vector=np.linspace(0.1,1.5,10)
#Las siguientes variables son para guardar la extinción mínima.
extincion_minima=1
valores_extincion=[]
#Y esta matriz es para los valores de la extinción de todos los puntos.
f=np.zeros((len(pi_B_to_A_vector),len(pi_A_to_B_vector)))
#Hay que rellenarla de esta forma para que la representación sea la correcta.
for i in range(len(pi_A_to_B_vector)):
    for j in range (len(pi_B_to_A_vector)):
        f[j,i]=prob_extincion_simple(10000,entornos,n_A,n_B,kappa_0_to_1,kappa_1_to_0,k_A,k_B,pi_A_to_B_vector[i],pi_B_to_A_vector[j])
        print(f[j,i],'valores', pi_A_to_B_vector[i],pi_B_to_A_vector[j])
        if f[j,i]<extincion_minima:
            extincion_minima=f[j,i]
            valores_extincion=[pi_A_to_B_vector[i],pi_B_to_A_vector[j]]    
print('La extinción mínima es',extincion_minima)
print('Para los valores',valores_extincion)

#Creamos un vector para representar el frente de Pareto.
pi_A_to_B_vector_pareto=np.linspace(0.1,3,1000)
pi_B_to_A_vector_pareto, pi_A_to_B_vector_pareto=frente_pareto(pi_A_to_B_vector_pareto)
#Lo representamos.
plt.plot(pi_A_to_B_vector_pareto,pi_B_to_A_vector_pareto,'.',color='blue',label='Frente de Pareto')
#Personalizamos la gráfica y representamos el mapa de extinción.
plt.xlim([0,3]); plt.ylim([0,1.5]); plt.legend(loc="upper right")
plt.pcolormesh(pi_A_to_B_vector, pi_B_to_A_vector,f, cmap='OrRd',shading='nearest')
cbar = plt.colorbar()
cbar.set_label('Extinción', fontsize=16);
plt.legend(loc="upper right")
plt.xlabel('$π_{A\_to\_B}$'); plt.ylabel('$π_{B\_to\_A}$'); plt.title('Extinción')
e = int(time.time() - start_time)
print('{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))'''