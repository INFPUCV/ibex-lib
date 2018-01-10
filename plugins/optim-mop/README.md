# Description

This plugin implements *ibexMop*, an 
interval branch & bound solver for **Nonlinear BiObjective Optimization** problems 
in the library [Ibex](https://github.com/ibex-team/ibex-lib).

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

* Includes an additional dynamic constraint for better defining the feasible 
objective region related to each box. This constraint is used 
by the filtering procedures improving the perfomance of the solver.

![Cy Comparison](https://i.imgur.com/yLIxyUV.png)
The *envelope* for the instance *kim* with *eps*=1. 
In the left side the strategy whitout the additional constraint. 
In the right side the strategy using the additional constraint which allows  
to approximate better the *envelope* of non-dominated solutions.

* Includes *NDSdist* a new strategy for selecting nodes. *NDSdist*
selects in each iteration the node/box maximizing its distance to the
upper envelope.
*NDSdist* has an anytime behaviour, i.e., it can return valid solutions
even if it is interrupted before it ends. See the figure below:

![Cy Comparison](https://i.imgur.com/uyZq6gB.png)
Comparison of the anytime behavior of the search strategies.
Figures show the envelope of the non-dominated set for the instances 
[*osy*](https://github.com/INFPUCV/ibex-lib/blob/master/plugins/optim-mop/benchs/MOP/osy.txt) 
after 100 iterations (top) and [*tan*](https://github.com/INFPUCV/ibex-lib/blob/master/plugins/optim-mop/benchs/MOP/tan.txt) 
after 50 iterations (down), 
using the [OC search strategy](http://www.sciencedirect.com/science/article/pii/S0377221716303824) (left) 
and the *NDSdist* search strategy (right).

* Includes a [inner polytope algorithm](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.653.5777&rep=rep1&type=pdf) 
for updating the upper envelope using feasible solutions. 
The algorithm constructs a feasible and convex polytope and then it finds 
two feasible vectors inside the polytope by minimizing a linearization of each one of the 
objective functions in the polytope. 
Then it finds a set of $n-2$ equidistant feasible solutions 
between this two vectors.

# Instalation

## Requirements

Optim-mop is a Plugin for the Ibex-Lib, so it's neccesary to have it before you could install it.
This repo includes it, but if you just want to add this pluggin for your own installation, you could download the ibex-mop.zip that is in this folder and add it to the plugins folder in your ibex root directory.

## Configure

Once you have the optim-mop in your plugins folder, you should go to the root folder of your ibex-lib and run the following line in your terminal.

```
./waf configure --with-optim --with-optim-mop --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex
```
```
./waf install
```

This will configure Ibex to use the Optimizer MOP plugin and then build it in:
```
<ibex-root>/__build__/plugins/optim-mop/
```

# Uses
```
./ibexmop {OPTIONS} [filename]
```
  OPTIONS:

      -h, --help                        Display this help menu
      -f[string], --filt=[string]       the filtering method
      --linear-relax=[string]           the linear relaxation method
      -b[string], --bis=[string]        the bisection method
      -s[string], --search=[string]     the search strategy
      --eps=[float]                     eps (the precision of the pareto front)
      --eps_x=[float]                   eps_x (the precision of the x boxes)
      -t[float], --time=[float]         timelimit
      --plot                            Save a python plot.
      --no-bisecty                      Do not bisect y variables.
      --cy-contract                     Contract using the box y+cy, w_ub=+inf.
      --cy-contract-full                Contract using the box y+cy.
      --eps-contract                    Contract using eps.
      --nb_ub_sols=[int]                Max number of solutions added by the
                                        inner-simplex
      --weight2=[float], --w2=[float]   Min distance between two non dominated
                                        points to be considered (default: 0.01)
      --min_ub_dist=[float]             Min distance between two non dominated
                                        points to be considered (default:
                                        eps/10)
      --hv                              Compute the hypervolume
      --trace                           Activate trace. Updates of loup/uplo are
                                        printed while minimizing.
      filename                          The name of the MINIBEX file.
      "--" can be used to terminate flag options and force all following
      arguments to be treated as positional options

### Filtering Method (-f):
 + hc4
 + acidhc4
### Linear Relaxation Method (--linear-relax):
 + compo
### Search Strategy (-s):
 + weighted_sum
 + NDSdist
 + diving-NDSdist
### Bisection Method (-b):
 + lsmear
 + largestfirst

### Benchs:
Inside benchs/MOP folder you could find Multi Objective Problems to solve.

## Examples:

To run an example just write this in your terminal inside the ibex root:
```
./__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/MOP/binh.txt --cy-contract --eps 1 -b largestfirst --nb_ub_sols 10 --w2 0.01
```





## Authors:
 - Ignacio Araya - <ignacio.araya@pucv.cl>
 - Damir Aliquintui - <damir.aliquintui.p@mail.ucv.cl>
 - Jose Campusano - <jose.campusano.c@mail.ucv.cl>
