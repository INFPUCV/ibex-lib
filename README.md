[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org

./waf configure --with-optim  --with-ampl --with-affine --prefix=. --gaol-dir= --optim-lib=soplex

./waf install

export PKG_CONFIG_PATH=/home/directorio_ibex/ibex-2.3.4/share/pkgconfig   -> Esto se tiene que hacer cada vez que se conecta a la maquina

In plugins/optim/examples:
make optimizer-mop

Then, to run an example:
./optimizer-mop test2.txt acidhc4 compo smearsumrel diving 1e-1 100 1

Para compilar todo en una sola linea yo hago lo siguiente:
1. ingreso al directorio raiz de ibex (solo la primera vez)
2. cd plugins/optim/examples (solo la primera vez)
3. cd -; sudo ./waf install; cd -; rm optimizer-mop; make optimizer-mop (cada vez que quiero re-compilar todo)

Y luego resolver un problema de ejemplo:
./optimizer-mop test2.txt acidhc4 compo smearsumrel diving 1e-1 100 1


TODO

Graficar resultados paso a paso:
  - Bisección (remplazar caja por dos subcajas)
  - Reducción o eliminación de caja (filtrado)
  - Actualización UB (graficar como función escalera)
  - Graficar recta lb en caja
  
Estudiar y agregar técnicas de selección de nodo (para comparar con DivingMOP)
- Estrategia del paper Constraint propagation using dominance in interval (2016)
http://ben-martin.fr/files/publications/2016/EJOR_2016.pdf)
Escoger nodo que minimiza: (z1.lb-z1_init.lb)/wid(z1_init) +  (z2.lb-z2_init.lb)/wid(z2_init)
Puede que baste con agregar función de comparación a CellSetBuffer
- Estrategias de acá https://drive.google.com/open?id=0B9JSHx01XN1rSXNOWXZXbEswWVE&authuser=0
No encontré más estrategias...
Para que quede una comparación pulenta yo agregaría:
- Escoger caja no dominada random
- Escoger caja no dominada con area máxima
- Caja que maximiza distancia a UB

Criterio de parada: distancia de caja a UB < epsilon 
(criterios que solo ven el tamaño de las cajas no permiten 
realizar una comparación justa)

Crear métodos que permitan mantener y actualizar set 
de segmentos no dominados (en principio para el UB). 
Segmentos se representan con conjunto de puntos.
- Agregar segmento (recibe dos puntos), debe actualizar los segmentos del set
Ver Vincent et al. (2013) Multiple objective branch and bound for mixed 0-1 linear programming: 
Corrections and improvements for the biobjective case






