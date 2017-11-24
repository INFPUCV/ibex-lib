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
./optimizer-mop ../benchs/MOP/osy.txt acidhc4 compo lsmear NDSdist 1e-6 10000 0 10 1e-7 0.0 1

Parametros:
./optimizer-mop instance filtering linear_relaxation bisector node_selection eps maxtime plot? nb_ub_sols min_ub_dist weight_finder rand_seed

weight_finder: finder minimize f1+weight*f2 and weight*f1+f2
nb_ub_sols: max number of solutions returned by the upperbounding method on each tree node
min_ub_dist: minimal distance between two ub points



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
  - tiempo, nodos, soluciones
  - upperbounding simplex (politopo)
  - Selecccion de nodos (maxdistance (mas lejos del upperbound))
  - Selecccion de nodos (diving decir que era buena para un objetivo y no para dos objetivos)
  - Metodo de caja box + cy (lo que mejora la precision w_lowerbound)
  - Metodo de caja box + cy (lo que mejora el filtrado w_upperbound) 

  - Compare hypervolumes and times between using the line z1+a*z2=w for improving the lower bound and the traditional
  method (using just the boxes lb)
  - Compare the new heuristic (max_distance) with other approaches. Show also the anytime behaviour.
  - Compare the upperbounding using the midpoint, simplex with n points and the dynamic version of simplex
  - Compare different approaches for selecting variable to bisect
  - Compare different contractor strategies (hc4, acidhc4, acidhc4+compo)

*Paper 2. Nonlinear biobjective optimization. Improving the precision of the nondominated set by using edges.* (contribuciones):
  - Definicion del ub_set usando segmentos factibles
  - Proponer algoritmo eficiente (busqueda binaria) para encontrar soluciones factibles asociadas a puntos en los segmentos
