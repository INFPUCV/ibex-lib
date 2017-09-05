[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org

./waf configure --with-optim  --with-ampl --with-affine --prefix=. --gaol-dir=  --soplex-path=../ibex-dag/soplex-1.7.2

./waf install

export PKG_CONFIG_PATH=/home/directorio_ibex/ibex-2.3.4/share/pkgconfig   -> Esto se tiene que hacer cada vez que se conecta a la maquina

In plugins/optim/examples:
make defaultoptimizer

Then, to run an example:
./defaultoptimizer ../benchs/easy/ex14_1_2.bch


Para compilar todo en una sola linea yo hago lo siguiente:
1. ingreso al directorio raiz de ibex (solo la primera vez)
2. cd plugins/optim/examples (solo la primera vez)
3. cd -; sudo ./waf install; cd -; rm defaultoptimizer; make defaultoptimizer (cada vez que quiero re-compilar todo)


Optimizer:
plugins/optim/src/strategy/ibex_Optimizer

Ejemplos de uso:
plugins/optim/examples

Una vez hecho los cambios.

1.- En la Raiz ./waf install

2.- In plugins/optim/examples:
rm defaultoptimizer
make defaultoptimizer

3.- Then, to run an example:
./defaultoptimizer ../benchs/easy/ex14_1_2.bch
