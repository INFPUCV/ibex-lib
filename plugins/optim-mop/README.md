# Description

Interval based solvers are commonly used for solving 
single-objective non linear optimization problems. 
Their reliability and increasing performance make them useful
when proofs of infeasibility and/or certification of solutions
are a must.
On the other hand, there exist only a few approaches dealing with
non linear optimization problems, when they consider multiple objectives. 

# Instalation

## Requirements

Optim-mop is a Plugin for the Ibex-Lib, so it's neccesary to have it before you could install it.
This repo includes it, but if you just want to add this pluggin for your own installation, you should add this folder to the plugins folder in your ibex directory.

Additionaly, you shoul have soplex installed in yoursystem.

## Configure

Once you have the optim-mop folder in your plugins folder, you should go to the root folder of your ibex-lib and run the following line in your terminal.

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

### Contact:

+ Author: Ignacio Araya
+ mail: ignacio.araya@pucv.cl
