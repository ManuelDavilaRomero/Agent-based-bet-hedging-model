# TFG-Creciendo-y-evolucionando-en-ambientes-fluctuantes
Los programas incluidos en este repositorio han sido desarrollados como parte del TFG Creciendo y evolucionado en ambientes fluctuantes para el grado de Físicas de la UCM. El autor íntegro ha sido el alumno Manuel Dávila Romero, ayudado por los consejos de los dos supervisores del TFG Luis Dinis y Francisco J. Cao.

Se ha usado el lenguaje de programación Python en su versión 3.7.9, con ayuda del entorno de Spyder en un ordenador MSI GL63 8RD con procesador Intel(R) Core(TM) i7-8750.

Aunque los códigos están comentados para facilitar su compresión, conviene hacer algunas anotaciones generales. El modeloagentesconextincion es el programa con cálculo de probabilidad de extinción, mientras que modeloagenteslogaritmos se centra en el estudio de propiedades del logaritmo de la población. Finalmente, el programa modeloagentesprimerapasada se ha utilizado para calcular el tiempo de primera pasada. Los códigos están formados por una serie de funciones, que modelizan la situación presentada en la sección 3.1 de la memoria del TFG. La función principal en todos ellos se denomina fluctuating_environments y calcula los cambios de ambiente y fenotipo, así como el crecimiento de la población. Al final del código, aparecen comentados los comandos que permiten reproducir los resultados presentados. 

-Para obtener las trayectorias de la sección 4.1 descomentamos en modeloagentesconextincion Trayectorias, que ejecuta fluctuating_environments 100 veces para distintos valores de las constantes de cambio de fenotipo.

-Para el análisis de la distribución del logaritmo de la población de la sección 4.2, descomentamos en modeloagenteslogaritmos Análisis del logaritmo de la población. Se calcula entonces la distribución del logaritmo de la población a 3 tiempos diferentes y, con el uso de una librería, se ajusta una gaussiana para cada uno de ellos.

-Para la extinción de la sección 4.3 descomentamos en el primer programa Extinción simple, si solo nos interesa un punto concreto, o Mapa de extinción, si lo que deseamos es reproducir figuras como la 5 o la 6. Como se desarrolla en la memoria, la probabilidad de extinción es el cociente entre el número de poblaciones que se extinguen y el número total de simulaciones, por lo que estos comandos básicamente ejecutan varias veces fluctuating_environments y almacenan cuántas trayectorias se absorben.

-Para el análisis teórico de la sección 4.4 ya hemos mencionado que modeloagentesprimerapasada es el programa que debemos usar para el estudio del tiempo de primera pasada. Si lo que nos interesa es la función de supervivencia o extinción teóricas recurrimos a modeloagenteslogaritmos y descomentamos Cálculo de la probabilidad de extinción teórica.

-Para estudiar la trayectoria en el plano (ln(n_A),ln(n_B)) como se hace en la sección 4.5, basta descomentar Trayectorias en modeloagenteslogaritmos.
