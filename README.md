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


TODO (13 de octubre)


TODO  (DA, MC): Generalizar CellFeasibleDiving con template (funcion de comparacion para nodos hermanos)
        Crear CellNS_SetMOP (CellBuffer) consistente en un set de cajas no dominadas 
        (ordenadas segun funcion de comparacion) y una liste de nodos dominados 
  
Graficar resultados (ub points y lb boxes)

