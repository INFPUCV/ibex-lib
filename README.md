[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org

Instalation
-----------

    sudo apt install python2.7 flex bison gcc g++ make pkg-config
    cp /usr/bin/python2.7 /usr/bin/python

download soplex version 1.7.2 from http://soplex.zib.de/

    tar -zxvf filename.tgz

    sudo apt-get install libz-dev zlib1g-dev

requirements for plot3.py
-------------------------

´´´´
  sudo apt-get install python3 python3-tk
  sudo apt-get install python3-virtualenv virtualenv
  virtualenv env --python=python3
  source env/bin/activate
  pip install -r requirements.txt

  ./waf configure --with-optim --with-optim-mop  --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex

  ./waf install
´´´´

Luego resolver un problema de ejemplo:

       __build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/binh.txt --cy-contract-full --eps-contract --eps 0.1 -b largestfirst --nb_ub_sols 2 --w2 0.01 --SEGMENTS

Para graficar resultados:

       python3 plugins/optim-mop/main/plot3.py
