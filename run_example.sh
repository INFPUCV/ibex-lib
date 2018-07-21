# THIS SCRIPT EXECUTE AN IBEX-MOP EXAMPLE

# configure
# ./waf configure --with-optim --with-optim-mop --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex

# install
./waf install

# run example
__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/binh.txt --cy-contract-full --eps-contract --eps 0.1 -b largestfirst --nb_ub_sols 2 --w2 0.01 --SEGMENTS
