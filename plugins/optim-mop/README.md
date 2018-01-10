# ibexMop

This plugin implements *ibexMop*, an 
interval branch & bound solver for **Nonlinear BiObjective Optimization** problems 
in the [Ibex library](https://github.com/ibex-team/ibex-lib).

*ibexMop* returns a set of solutions X and its images Y
guaranteeing a maximal distance *eps* between
any *non-dominated* feasible vector and the returned set Y. 

*ibexMop* constructs an **envelope** for the non-dominated set
by following a branch & bound strategy starting with an initial *box* (containing the variable domains) 
and building a search tree. In each iteration of the algorithm, 
a node is selected and treated by classical *filtering*, *upper-bounding* 
and *splitting* techniques. 

*ibexMop* includes some methods to take into account the upper bound of the 
*envelope* for filtering dominated solutions (e.g., well-known discarding tests).

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

## Installation

### Requirements

*ibexMop* is a plugin of the Ibex Library, so if you already have the library, 
you only have to add this folder 
([optim-mop.zip](https://github.com/INFPUCV/ibex-lib/blob/master/plugins/optim-mop/optim-mop.zip)) in the plugin
folder of Ibex. 
Otherwise you can download the entire library (including ibexMop) from [here](https://github.com/INFPUCV/ibex-lib).

### Configure

Once you have the plugin in the plugins' folder, you should go to the root folder of the  Ibex library 
and run the following line in your terminal.

```
./waf configure --with-optim --with-optim-mop --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex
```
```
./waf install
```

This will configure Ibex to use the *IbexMop* plugin and then will build it in:
```
<ibex-root>/__build__/plugins/optim-mop/
```

## Uses
```
./__build__/plugins/optim-mop/ibexmop {OPTIONS} [filename]
```
  OPTIONS (the most important ones):

      -h, --help                        Display this help menu
      -f[string], --filt=[string]       the filtering method (default: acidhc4)
      --linear-relax=[string]           the linear relaxation method (default: compo)
      -b[string], --bis=[string]        the bisection method (default: largestfirst)
      -s[string], --search=[string]     the search strategy (default: NDSdist)
      --eps=[float]                     the precision (default: 0.01)
      -t[float], --time=[float]         timelimit (default: 100)
      --cy-contract-full                Contract using the additional constraint.
      --eps-contract                    Contract using eps.
      --nb_ub_sols=[int]                Max number of solutions added by the
                                        inner-polytope algorithm (default: 50)
      filename                          The name of the MINIBEX file.
      "--" can be used to terminate flag options and force all following
      arguments to be treated as positional options

### Filtering Method (-f):
 + hc4
 + acidhc4
### Linear Relaxation Method (--linear-relax):
 + no
 + compo (a method combining two linearization techniques: AF2 and XNewton)
### Search Strategy (-s):
 + weighted_sum (or the [OC search strategy](http://www.sciencedirect.com/science/article/pii/S0377221716303824))
 + NDSdist
### Bisection Method (-b):
 + largestfirst

## Run an example:

To run an example with the default parmeters just write this line in your terminal in the root directory of the Ibex library:
```
./__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/binh.txt
```

## Format of the instances (Minibex):

Instances can be written in the [Minibex language](http://www.ibex-lib.org/doc/minibex.html), 
considering that the objectives *must correspond* to the first two constraints with the following syntax:
```
Constraints
<expression of the first objective function> = z1;
<expression of the second objective function> = z2;
// other constraints...
```
You can see some examples in [benchs](https://github.com/INFPUCV/ibex-lib/tree/master/plugins/optim-mop/benchs).

## Authors:
 - Ignacio Araya - <ignacio.araya@pucv.cl>
 - Damir Aliquintui - <damir.aliquintui.p@mail.ucv.cl>
 - Jose Campusano - <jose.campusano.c@mail.ucv.cl>
