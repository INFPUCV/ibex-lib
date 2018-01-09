# Description

Interval based solvers are commonly used for solving 
single-objective non linear optimization problems. 
Their reliability and increasing performance make them useful
when proofs of infeasibility and/or certification of solutions
are a must.
On the other hand, there exist only a few approaches dealing with
non linear optimization problems, when they consider multiple objectives.

Interval branch & bound solvers start with an initial box 
and build a search tree. 
In each iteration a node is selected and treated by *filtering* 
and *upper-bounding* procedures.
Filtering (or contraction) consists in removing inconsistent 
values from the bounds of the box, while the upper-bounding 
procedures attempt to find good 
feasible solutions for improving the upper envelope.
Some filtering methods take into account the upper 
envelope for filtering dominated solutions (e.g., discarding 
tests~).
If the box has not been discarded by the filtering process, 
it is split into two sub-boxes by dividing the domain of one 
variable and generating two child nodes in the search tree.
The procedure iterates until a termination criteria is reached.
In current solvers, boxes are no longer treated when their
sizes reach a given precision.

# Instalation

## Requirements

Optim-mop is a Plugin for the Ibex-Lib, so it's neccesary to have it before you could install it.
This repo includes it, but if you just want to add this pluggin for your own installation, you should add this folder to the plugins folder in your ibex directory.

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

### Contact:

+ Author: Ignacio Araya @rilianx
+ mail: ignacio.araya@pucv.cl
