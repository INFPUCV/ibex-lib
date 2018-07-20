# THIS SCRIPT EXECUTE AN IBEX-MOP EXAMPLE

# configure
# ./waf configure --with-optim --with-optim-mop --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex

# install
./waf install

# run example
__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/binh.txt --cy-contract --eps 0.05 -b largestfirst --nb_ub_sols 10 --plot --w2 0.01 --hamb | grep -e "aux\|PROCESS_NODE\|error\|dist\|bis\|point lb"

__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/binh.txt --cy-contract-full --eps-contract --eps 0.05 -b largestfirst --nb_ub_sols 10 --plot --w2 0.1 --hamb | grep -e "aux\|PROCESS_NODE\|error\|dist\|bis\|point lb"
