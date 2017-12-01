[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org

Instalation
-----------

./waf configure --with-optim  --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex

./waf install

export PKG_CONFIG_PATH=/home/directorio_ibex/ibex-2.3.4/share/pkgconfig   -> Esto se tiene que hacer cada vez que se conecta a la maquina

export PKG_CONFIG_PATH=/home/iaraya/github/ibex/ibex-dev-dag/ibex-lib/share/pkgconfig


In plugins/optim/examples:
make optimizer-mop


Para compilar todo en una sola linea yo hago lo siguiente:
1. ingreso al directorio raiz de ibex (solo la primera vez)
2. cd plugins/optim/examples (solo la primera vez)
3. cd -; sudo ./waf install; cd -; rm optimizer-mop; make optimizer-mop (cada vez que quiero re-compilar todo)

Y luego resolver un problema de ejemplo:
./optimizer-mop ../benchs/MOP/mop-7.txt --cy-contract --eps 0.02 -b lsmear --nb_ub_sols 10 --plot --w2 0.00


TODO
----

Graficar resultados on-the-fly apretando tecla para avanzar **(Matias)**:
  - [ ] Tener la opcion de mostrar UB como funcion escalonada/puntos
  - [x] Graficar recta lb dentro de cajas: z1 + a*z2=w_lb
  - [ ] Agregar parametros al ejecutable: mostrar plot (paso a paso)

Técnicas de selección de nodo:
  - [x] [OC](http://ben-martin.fr/files/publications/2016/EJOR_2016.pdf): min (z1.lb-z1_init.lb)/wid(z1_init) +  (z2.lb-z2_init.lb)/wid(z2_init)
  - [x] SR1, [OC1](https://tel.archives-ouvertes.fr/tel-01146856/document): Min lb1
  - [x] [OC2](https://tel.archives-ouvertes.fr/tel-01146856/document): Min lb2
  - [x] [OC3](https://tel.archives-ouvertes.fr/tel-01146856/document): Min lb1 + lb2
  - [x] [OC4](https://tel.archives-ouvertes.fr/tel-01146856/document): Decreasing value of
  hypervolume of the point y (considering the initial values for z1_ub and z2_ub)
  - [ ] [OC5](https://tel.archives-ouvertes.fr/tel-01146856/document): Decreasing box size
  of boxes such that lb is not dominated **(Damir's CellNSSet)**
  - [ ] Escoger caja random del Nondominated Set **(Damir)**
  - [x] Escoger caja que maximiza la distancia a UB (Optimizada usando Pqueue)
  - [x] Diving compatible con los metodos anteriores
  - [ ] Que hacer cuando aun no hay upperbounds?

Calculo de distancia:
  - [x] Modificar funcion de calculo de distancia usando ideas de Damir

Criterio de termino:
  - [x] Definir criterio relativo: abs_prec = rel_prec * min(wid(z1), wid(z2)) --> **repensar**
  - [x] Calcular hipervolumen relativo de la solucion
  - [ ] Definicion del lowerbound tiene errores aun (limitarse a reparar errores graves por ahora)


Biseccion:
  - [x] Adaptar LSmear (tecnica de biseccion)

Discarding boxes:
  - [x] Lowerbounding usando restriccion auxiliar z1+a*z2=w
  - [x] w_lb delimitado por puntos UB
  - [x] pendiente igual a pendiente entre puntos extremos
  - [x] Implementar monotonicity test (FT) from [here](https://link.springer.com/content/pdf/10.1007%2Fs10589-007-9135-8.pdf)
  - [ ] Filtrar cajas del buffer

**(Ignacio)** Upperbounding:
  - [x] Criterio dinamico para establecer cantidad de puntos que se generan
  - [x] Adaptar para trabajar con ecuaciones (Hansen-Sengupta)
  - [ ] Bug: Problema con PdcHansenFeasibility (LoupFinderMOP) para problema test.bch
  - [ ] Encontrar recta factible en x usando simplex,
  para luego obtener segmentos upperbound de la curva asociada en y
  - [ ] Implementar metodos para manejar el set de segmentos no dominados (nuevo UB)

Definicion del lowerbound (y eventualmente UB):
  - [x] Algoritmo para definir segmentos LB o UB **reparar bugs**
  


[Algoritmo para encontrar interseccion entre 2 segmentos](https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect)

Ver Vincent et al. (2013) Multiple objective branch and bound for mixed 0-1 linear programming:
Corrections and improvements for the biobjective case

Preparar experimentos:
  - [x] Agregar todos los benchmarks de este [paper](http://ben-martin.fr/files/publications/2016/EJOR_2016.pdf)
  - [x] Agregar opcion de setear cantidad de soluciones del upperbounding (nb_ub_sols)
  - [x] Agregar opcion para mostrar/no mostrar plot (plot?)
  - [x] Agregar opcion para modificar distancia minima aceptada entre soluciones factibles (min_ub_dist)

**Estructura de papers.**

*Paper 1. Nonlinear biobjective optimization. Improvements to interval Branch&Bounds algorithms*  (contribuciones):
  - Propiedad de las soluciones (ub):
  any feasible solution x is eps-dominated by at least one solution x' in the ub set, i.e.,
  for all x feasible, there exists at least one x' in ub_set, such that: f1(x') <= f1(x) + eps  and f2(x') <= f2(x) + eps
  - Crierios de seleccion del siguiente nodo (max_distance + diving?)
    - Explicar algoritmo ub_distance en detalle (usando pqueue y recalculo de distancias)
  - Upperbounding usando simplex (min f1 + min f2 + midpoint)
  - Discarding boxes by using an auxiliary constraint: z1+a*z2=w


**Experiments**
  - determinar buen valor para la precision para obtener resultados en todas las instancias en tiempo razonable
  (precision relativa al rango de las soluciones no dominadas)
  - tiempo, nodos, soluciones
  - Estrategia basica (std): -f hc4 -b largestfirst -s weighted_sum --nb_ub_sols=1
  - Estrategia full contractor (fullctc): -f acidhc4 --lr=compo
  - upperbounding simplex (politopo) + (acidhc4 + compo)
  	- std
  	- std + --nb_ub_sols=X, X in {3, 5, 10, 50, 100}
  	- fullctc + --nb_ub_sols=X, X in {3, 5, 10, 50, 100}
  - Metodo de caja box + cy (lo que mejora la precision w_lowerbound) 
    - fullctc + --nb_ub_sols=best_X
    - fullctc + --nb_ub_sols=best_X + --cy-contract
  - Metodo de caja box + cy (lo que mejora el filtrado w_upperbound) 
    - fullctc + --nb_ub_sols=best_X + --cy-contract-full
  - Comparar estrategias de seleccion de nodo (NDSdist - diving-NDSdist)
    - -s weighted_sum, -s NDSdist, -s diving-NDSdist
  - Comparar estrategias de biseccion (lsmear - largestfirst)
    - -b largestfirst, -b lsmear


*Paper 2. Nonlinear biobjective optimization. Improving the precision of the nondominated set by using edges.* (contribuciones):
  - Definicion del ub_set usando segmentos factibles
  - Proponer algoritmo eficiente (busqueda binaria) para encontrar soluciones factibles asociadas a puntos en los segmentos
