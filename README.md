[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org

Instalation
-----------

./waf configure --with-optim  --with-ampl --with-affine --prefix=. --gaol-dir= --optim-lib=soplex

./waf install

export PKG_CONFIG_PATH=/home/directorio_ibex/ibex-2.3.4/share/pkgconfig   -> Esto se tiene que hacer cada vez que se conecta a la maquina

export PKG_CONFIG_PATH=/home/iaraya/github/ibex/ibex-dev-dag/ibex-lib/share/pkgconfig

In plugins/optim/examples:
make optimizer-mop

Then, to run an example:
./optimizer-mop test2.txt acidhc4 compo smearsumrel diving-minlb 1e-1 100 1

Para compilar todo en una sola linea yo hago lo siguiente:
1. ingreso al directorio raiz de ibex (solo la primera vez)
2. cd plugins/optim/examples (solo la primera vez)
3. cd -; sudo ./waf install; cd -; rm optimizer-mop; make optimizer-mop (cada vez que quiero re-compilar todo)

Y luego resolver un problema de ejemplo:
./optimizer-mop test2.txt acidhc4 compo smearsumrel diving-NDSdist 1e-1 100 1


TODO
----

Graficar resultados on-the-fly apretando tecla para avanzar **(Matias)**:
  - [x] Destacar caja seleccionada en cada iteracion (top)
  - [ ] Tener la opcion de mostrar UB como funcion escalonada/puntos
  - [x] Mostrar conjunto de cajas UB
  - [ ] Mostrar conjunto de cajas Sout (si es que hubiera)
  - [x] Graficar recta lb dentro de cajas: z1 + a*z2=w_lb

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

Criterio de termino:
  - [x] Definir criterio relativo: abs_prec = rel_prec * min(wid(z1), wid(z2))
  - [ ] Calcular hipervolumen relativo de la solucion

Biseccion:
  - [x] Adaptar LSmear (tecnica de biseccion)
  - Algunos metodos de biseccion no bisectan algunas variables

Discarding boxes:
  - [x] Lowerbounding usando restriccion auxiliar z1+a*z2=w  
  - [x] Usar esta distancia para heuristica de seleccion de caja y criterio de termino

**(Ignacio, Damir)** Upperbounding:
  - [x] Criterio dinamico para establecer cantidad de puntos que se generan
  - [x] Adaptar para trabajar con ecuaciones (Hansen-Sengupta)
  - [ ] Encontrar recta factible en x usando simplex,
  para luego obtener segmentos upperbound de la curva asociada en y
  - [ ] Implementar metodos para manejar el set de segmentos no dominados (nuevo UB)



(Ignacio)
Definir estructura de papers.
Paper 1. Nonlinear biobjective optimization. Improvements to interval Branch&Bounds algorithms  (contribuciones):
  - Propiedad de las soluciones (ub):
  any feasible solution x is eps-dominated by at least one solution x' in the ub set, i.e.,
  for all x feasible, there exists at least one x' in ub_set, such that: f1(x') <= f1(x) + eps  and f2(x') <= f2(x) + eps
  - Crierios de seleccion del siguiente nodo (ub_distance + diving)
  - Upperbounding usando simplex (min f1 + min f2 + midpoint)
  - Discarding boxes by using an auxiliary constraint: z1+a*z2=w
  
  
Experiments

  - Compare hypervolumes and times between using the line z1+a*z2=w for improving the lower bound and the traditional 
  method (using just the boxes lb)
  - Compare the new heuristic (max_distance) with other approaches. Compare also the anytime behaviour.
  - Compare the upperbounding using the midpoint, simplex with n points and the dynamic version of simplex

Paper 2. Nonlinear biobjective optimization. Improving the precision of the nondominated set by using edges. (contribuciones):
  - Definicion del ub_set usando segmentos factibles
  - Proponer algoritmo eficiente (busqueda binaria) para encontrar soluciones factibles asociadas a puntos en los segmentos



(Ignacio)
Crear métodos que permitan mantener y actualizar set (para paper 2)
de segmentos no dominados (en principio para el UB).
Segmentos se representan con conjunto de puntos.
- Agregar segmento (recibe dos puntos), debe actualizar los segmentos del set

Notar que x aumenta e y disminuye en el ub_set.


    insert_segment(p1, p2)
      p1+ <- (p1.x,inf)
      p2+ <- (p2.y,inf)
      v1 <- lb_x de p1 en ub_set
      v2 <- next(ub_set)
      in <- false
      if s <- intersect(v1-v2, p1+-p1)
      in <- true
      add_point(s), add_point(p1)
      delete_point(v2)  

      while(v1.y>p2.y)

        if s <- intersect(v1-v2, p1-p2)
          in=!in
          add_point(s)

        if(in) delete point(v2)

        if v2.y > p2.y && s <- intersect(v1-v2, p2-p2+)
          in<-false;
          add_point(s); add_point(p2)

        v1 <- v2
        v2 <- next(ub_set)  

[Algoritmo para encontrar interseccion entre 2 segmentos](https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect)

    point intersect(p, p2, q, q2)
      if p.x=-inf
      if q.x=-inf return q
        else return NULL
      if p2.y=-inf
      if q2.y=-inf return q2
        else return NULL

      r=p2-p
      s=q2-q
      //now we find a solution for the equation p+tr = q+us,
      t=(q-p) x s/(r x s)

      if r x s!=0 and t in [0,1]
        return p+tr
      else
        return NULL

Ver Vincent et al. (2013) Multiple objective branch and bound for mixed 0-1 linear programming:
Corrections and improvements for the biobjective case
