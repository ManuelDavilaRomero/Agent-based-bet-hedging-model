# Estrategias óptimas en ambientes fluctuantes/ Extinction in agent-based and collective models of bet-heding 

Español 

Las simulaciones se han desarrollado usando el lenguaje de programación de *Python* en su versión 3.7.9 y 3.8.10, con ayuda del entorno de *Spyder* usando principalmente un Legion 5 15IAH7H Laptop (Lenovo) - Type 82RB con procesador 12th Gen Intel(R) Core(TM) i7-12700H   2.30 GHz. 

Los tres programas tienen una parte común. Incluyen tres funciones encargadas de generar  los números aleatorios necesarios para los procesos de cambio de ambiente (**cambio_ambiente**), reproducción/muerte  (**reproduccion_muerte_bacteria**) y el cambio de fenotipo (**cambio_fenotipo**). Para calcular el paso temporal de la simulación y ejecutar estos programas en cada paso, se crea la función principal **fluacting_environments**. Los resultados presentados en el trabajo se incluyen como *scripts* que basta descomentar para reproducibilidad.
- **agentmodelextinction**. Este programa utiliza la cota superior. Se utiliza para calcular la probabilidad de extinción, tanto de puntos aislados como mapas de extinción. 
- **agentmodelfirstpassagetime**. Se centra en el estudio de los tiempos de primera pasada, por lo que utiliza el logaritmo de la población.
- **agentmodellogarithm**. El análisis de las trayectorias de la Sección 4.1 se realiza con este programa. No incluye la cota superior, pues en caso contrario, la distribución de las trayectorias no sería correcta. 

English 

The simulations have been developed using the Python programming language in its versions 3.7.9 and 3.8.10, with the assistance of the Spyder environment, primarily running on a Legion 5 15IAH7H Laptop (Lenovo) - Type 82RB with a 12th Gen Intel(R) Core(TM) i7-12700H 2.30 GHz processor.

All three programs share a common part. They include three functions responsible for generating the random numbers necessary for the processes of environmental change (**environment_change**), reproduction/death (**reproduction_death_bacteria**), and phenotype change (**phenotype_change**). To calculate the time step of the simulation and execute these programs at each step, the main function **fluacting_environments** is created. The results presented in the paper are included as scripts that only need to be uncommented for reproducibility.

- **agentmodelextinction**. This program utilizes the upper bound. It is used to calculate the probability of extinction, both for isolated points and extinction maps.
- **agentmodelfirstpassagetime**. It focuses on studying the first passage times, thus utilizing the logarithm of the population.
- **agentmodellogarithm**. The analysis of trajectories is performed with this program. It does not include the upper bound because otherwise, the distribution of trajectories would not be correct.

