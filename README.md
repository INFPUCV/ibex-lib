# Ibex-MOP

*ibexMop* is an interval branch & bound solver for **Nonlinear BiObjective Optimization** problems.

*ibexMop* constructs an **upper envelope** for the non-dominated set
by following a branch & bound strategy starting with an initial *box* (containing the variable domains)
and building a search tree. In each iteration of the algorithm,
a node is selected and treated by *filtering*, *upper-bounding*
and *splitting* techniques.

The solver returns a set of vectors Y' guaranteeing a maximal distance *eps* between
any *non-dominated* feasible vector and the returned set Y'. It includes 
some methods to take into account the *upper envelope* 
for filtering dominated solutions (e.g., dominance peeler).

*ibexMop* includes three upper bounding methods:

  * (inner_polytope) The upper envelope Y' is represented by using a dominace-free set of *vectors*. An
  [inner polytope algorithm](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.653.5777&rep=rep1&type=pdf)
  is used for finding feasible solutions. The algorithm constructs a feasible and convex polytope and then it finds
two feasible vectors inside this polytope by minimizing a linearization of each one of the objective functions.
Then it finds a set of N-2 (by default N=50) equidistant feasible solutions.

  * (ub1) The upper envelope Y' is represented by using a set of *upper line segments*. 
  It is warrantied that every vector in any segment is epsilon-dominated 
  by at least one feasible solution. For finding upper line segments, the 
  algorithm first finds a feasible line segment in the domain space using an inner polytope algorithm 
  and then it generates a line segment passing over the image of this feasible line in the objective space. 
  
  * (ub2) The upper envelope Y' is also represented by using a set of upper line segments. 
  The algorithm performs the same procedure than ub1 for finding upper line segments, however
  it goes further and it is able to find the whole upper envelope related to the feasible line segment.

![ub1](https://i.imgur.com/H6zAwpO.png)
Example of an upper line segment.
  
![Cy Comparison](https://i.imgur.com/Wphf10d.png)
Example of using ub2 for finding an upper envelope for the blue feasible curve in the objective space.

## Download
````
git clone https://github.com/INFPUCV/ibex-lib.git
git checkout -t origin/master-mop
````

## Installation

````
    sudo apt install python3 flex bison gcc g++ make pkg-config libz-dev zlib1g-dev python3-tk
    pip3 install -r requirements.txt
    
    ./waf configure --with-optim --with-optim-mop --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex
    ./waf install
````      
       
## Main options
```
./__build__/plugins/optim-mop/ibexmop {OPTIONS} [filename]
```
  OPTIONS:

      -h, --help                        Display this help menu
      -f[string], --filt=[string]       the filtering method (default: acidhc4)
      --linear-relax=[string]           the linear relaxation method (default: compo)
      -b[string], --bis=[string]        the bisection method (default: largestfirst)
      -s[string], --search=[string]     the search strategy (default: NDSdist)
      --eps=[float]                     the desired precision (default: 0.01)
      -t[float], --timelimit=[float]    timelimit (default: 100)
      --cy-contract-full                Contract using the additional constraint cy.
      --eps-contract                    Contract using eps.
      --N=[int]                         Max number of solutions added by the inner-polytope algorithm (default: 50)
      --cy-contract                     Contract using the box y+cy, w_ub=+inf.
      --ub=[string]                     Upper bounding strategy (default: ub2).
      filename                          The name of the MINIBEX file.

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
 + roundrobin
### Upper bounding method (-ub):
 + inner_polytope
 + ub1
 + ub2

## Run an example:

     ./__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/binh.txt --eps-contract --ub=ub1 --eps=0.001 --cy-contract-full

For plotting the non-dominated vectors returned by the solver

     python3 plugins/optim-mop/main/plot3.py


## Format of the instances (Minibex):

Instances can be written in the [Minibex language](http://www.ibex-lib.org/doc/minibex.html),
considering that the objectives *must correspond* to the first two constraints with the following syntax:
```
Constraints
<expression of the first objective function> = z1;
<expression of the second objective function> = z2;
// other constraints...
```
You can see the set of instances in [plugins/optim-mop/benchs](https://github.com/INFPUCV/ibex-lib/tree/master-mop/plugins/optim-mop/benchs).


## Authors:
 - Ignacio Araya - <ignacio.araya@pucv.cl>
 - Damir Aliquintui - <damir.aliquintui.p@mail.ucv.cl>
 - Jose Campusano - <jose.campusano.c@mail.ucv.cl>
