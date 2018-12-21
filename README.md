[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org

Installation
-----------

````
    sudo apt install python3 flex bison gcc g++ make pkg-config libz-dev zlib1g-dev python3-tk
    pip3 install -r requirements.txt
    
    ./waf configure --with-optim --with-optim-mop --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex
    ./waf install
````

Luego resolver un problema de ejemplo:

       ./__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/binh.txt --eps-contract  --HAMBURGER --eps=0.001 --cy-contract-full

Para graficar resultados:

       python3 plugins/optim-mop/main/plot3.py
       
Ibex-MOP
--------

*ibexMop* is an interval branch & bound solver for **Nonlinear BiObjective Optimization** problems.

[comment]: <> *ibexMop* returns a set of solutions X and its images Y
[comment]: <> guaranteeing a maximal distance *eps* between
[comment]: <> any *non-dominated* feasible vector and the returned set Y.

*ibexMop* constructs an **upper envelope** for the non-dominated set
by following a branch & bound strategy starting with an initial *box* (containing the variable domains)
and building a search tree. In each iteration of the algorithm,
a node is selected and treated by classical *filtering*, *upper-bounding*
and *splitting* techniques.

*ibexMop* includes some methods to take into account the upper bound of the
*envelope* for filtering dominated solutions (e.g., dominance peeler).

*ibexMop* also offers several improvements related to other NLBOO algorithms:

* Uses a termination criteria directly related with the
precision of the *envelope*.

* Includes an *additional dynamic constraint* **cy** for better defining the feasible
objective region related to each box. This constraint is used
by the filtering procedures improving the perfomance of the solver.

![Cy Comparison](https://i.imgur.com/yLIxyUV.png)
The *envelope* for the instance *kim* with *eps*=1.
In the left side the strategy without the additional constraint **cy**.
In the right side the strategy using **cy** which allows
to approximate better the *envelope* of non-dominated solutions.
Red points corresponds to the set of feasible vectors found by the strategies.
Remark that no feasible vector can be found under or left the lower envelope (blue segments).

* Includes *NDSdist*, a new strategy for selecting nodes. *NDSdist*
selects in each iteration the node/box maximizing its distance to the
upper envelope.
*NDSdist* has an *anytime behaviour*, i.e., it can return valid solutions
even if it is interrupted before it ends. See the figure below:

![Cy Comparison](https://i.imgur.com/uyZq6gB.png)
Comparison of the anytime behavior of the search strategies.
Figures show the envelope of the non-dominated set for the instances
[*osy*](https://github.com/INFPUCV/ibex-lib/blob/master/plugins/optim-mop/benchs/osy.txt)
after 100 iterations (top) and [*tan*](https://github.com/INFPUCV/ibex-lib/blob/master/plugins/optim-mop/benchs/tan.txt)
after 50 iterations (down),
using the [OC search strategy](http://www.sciencedirect.com/science/article/pii/S0377221716303824) (left)
and the *NDSdist* search strategy (right).

* Includes a [inner polytope algorithm](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.653.5777&rep=rep1&type=pdf)
for finding feasible solutions.
The algorithm constructs a feasible and convex polytope and then it finds
two feasible vectors inside this polytope by minimizing a linearization of each one of the
objective functions.
Then it finds a set of $n-2$ equidistant feasible solutions
between this two vectors.
