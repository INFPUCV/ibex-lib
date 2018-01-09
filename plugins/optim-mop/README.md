# Description

This plugin implements a **Nonlinear BiObjective Optimization** (NLBOO) interval branch & bound 
Solver in the library [Ibex](https://github.com/ibex-team/ibex-lib).

The solver finds a thin envelope
containing the set of non-dominated vectors and a set of feasible 
solutions delimiting the upper bound of the envelope.

It follows a branch & bound strategy starting with an initial box 
and building a search tree. In each iteration a node is selected and treated by *filtering* 
and *upper-bounding* procedures. The plugin include some methods to take into account the upper 
envelope for filtering dominated solutions.
If the box has not been discarded by the filtering process, 
it is split into two sub-boxes by dividing the domain of one 
variable and generating two child nodes in the search tree.
The procedure iterates until a termination criteria is reached.

The plugin offers several improvements related to other algorithms:

* It uses a termination criteria directly related with the 
precision of the envelope containing the non-dominated vectors.

* We include an additional constraint for better defining the feasible 
objective region related to each box. This constraint is used 
by the filtering procedures improving the perfomance of the solver.

* In each iteration we select the node maximizing the distance to
the upper envelope of the non-dominated vectors.
Our criteria has an anytime behaviour, i.e., it can return valid solutions
even if it is interrupted before it ends. Also, the search strategy
allows to improve the precision of the envelope in a homogeneous and 
quite uniform way.

* For upperbounding we use a inner polytope algorithm used for monobjective optimization. 
The algorithm constructs a feasible and convex polytope and then it finds 
a feasible solution inside the polytope minimizing a linearization of the 
objective function. We adapted the algorithm for finding n feasible solutions: two solution vectors minimizing each one of the 
linearized objectives and a set of $n-2$ equidistant feasible solutions 
between this two vectors. As the first two vectors are inside the
feasible and convex polytope, we ensure that any solution between these
two vectors is also feasible.


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

![Cy Comparison](https://i.imgur.com/yLIxyUV.png)
A part of the set Y and the lower envelope for the instance kim with e= 1 and n= 50. In the left side the reference strategy, whitout the cy-contractor (79 boxes for reaching the precision). In the right side the strategy using the constraint cy which allows us to approximate better the lower envelope (only 19 boxes for reaching the precision).

![Cy Comparison](https://i.imgur.com/uyZq6gB.png)

Comparison of the anytime behavior of the search strategies.
Figures show the envelope of the non-dominated set for the instances osy after 100 iterations (top) and tan after 50 iterations (down), using the weighted sum search strategy (left) and the NDSdist search strategy (right).

## Authors:
 - Ignacio Araya - <ignacio.araya@pucv.cl>
 - Damir Aliquintui - <damir.aliquintui.p@mail.ucv.cl>
 - Jose Campusano - <jose.campusano.c@mail.ucv.cl>
