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

       ./__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/binh.txt --eps-contract --HAMBURGER --eps=0.001 --cy-contract-full

Para graficar resultados:

       python3 plugins/optim-mop/main/plot3.py
       
Ibex-MOP
--------

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

*ibexMop* includes three variants for representing the upper envelope Y' and performing the upper-bounding:

  * The upper envelope is represented by using a dominace-free set of *vectors*. An
  [inner polytope algorithm](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.653.5777&rep=rep1&type=pdf)
  is used for finding feasible solutions. The algorithm constructs a feasible and convex polytope and then it finds
two feasible vectors inside this polytope by minimizing a linearization of each one of the objective functions.
Then it finds a set of N-2 (by default N=50) equidistant feasible solutions.

  * (ub1) The upper envelope is represented by using a set of *upper line segments*. 
  It is warrantied that every vector in any segment is epsilon-dominated 
  by at least one feasible solution. For finding upper line segments, the 
  algorithm first finds a feasible line segment in the domain space using an inner polytope algorithm 
  and then it generates a line segment passing over the image of this feasible line in the objective space. 
  
  * (ub2) The upper envelope is also represented by using a set of upper line segments. 
  The algorithm performs the same procedure than ub1 for finding upper line segments, however
  it goes further and is able to find the whole upper envelope related to the feasible line segment.

![ub1](https://imgur.com/H6zAwpO)
Example of an upper line segment.
  
![Cy Comparison](https://imgur.com/Wphf10d)
Example of using ub2 for finding an upper envelope for the blue feasible curve in the objective space.

