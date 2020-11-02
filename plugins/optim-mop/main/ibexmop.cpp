//============================================================================
//                                  I B E X
// File        : ibexmop.cpp
// Author      : Ignacio Araya, Matias Campusano, Damir Aliquintui
// Copyright   : PUCV (Chile)
// License     : See the LICENSE file
// Created     : Jan 01, 2017
// Last Update : Jan 01, 2017
//============================================================================


#include "ibex.h"
#include "args.hxx"

// Server side C/C++ program to demonstrate Socket programming
#include <unistd.h>
#include <stdio.h>
#include <sys/socket.h>
#include <stdlib.h>
#include <netinet/in.h>
#include <string.h>
#include <cstring>
#include <set> 
#include <iterator>
#include <vector>


#ifndef _IBEX_WITH_OPTIM_MOP_
#error "You need the plugin Optim MOP to run this example."
#endif

const double default_relax_ratio = 0.2;

using namespace std;
using namespace ibex;

string Print_Vector(vector<double> Vec) 
{ 
	string vector = "";
    for (int i = 0; i < Vec.size(); i++) { 
		vector = vector + std::to_string(Vec[i]);
		if( i == 0){
			vector = vector + ",";
		}else{
			vector = vector + "/";
		}
    } 
    return vector; 
} 


struct graph_compare {
    bool operator() (const vector<double>& a, const vector<double>& b) const {
        if( a[0] != b[0]){
            return a[0] < b[0];
        }else{
            return a[1] > b[1];
        }
        return true;
    }
};

int main(int argc, char** argv){


	try {

	if (argc<1) {
		cerr << "usage: optimizer04 filename filtering linear_relaxation bisection strategy prec timelimit "  << endl;
		exit(1);
	}

	char fns[4] = "fns";

	args::ArgumentParser parser("********* IbexMop *********.", "Solve a NLBOOP (Minibex file).");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<std::string> _filtering(parser, "string", "the filtering method (default: acidhc4)", {'f', "filt"});
	args::ValueFlag<std::string> _linear_relax(parser, "string", "the linear relaxation method (default: compo)", {"linear-relax"});
	args::ValueFlag<std::string> _bisector(parser, "string", "the bisection method (default: largestfirst)", {'b', "bis"});
	args::ValueFlag<std::string> _strategy(parser, "string", "the search strategy (default: NDSdist)", {'s', "search"});
	args::ValueFlag<double> _eps(parser, "float", "the desired precision (default: 0.01)", {"eps"});
	args::ValueFlag<double> _eps_rpm(parser, "float", "precision of the relative point method (default: 1e-5)", {"eps_rpm"});
	args::ValueFlag<double> _eps_r(parser, "float", "the desired relative precision (default: 0.01)", {"eps_r"});
	args::ValueFlag<double> _timelimit(parser, "float", "timelimit (default: 100)", {'t',"timelimit"});
	args::Flag _cy_contract(parser, "cy-contract", "Contract using the box y+cy, w_ub=+inf.", {"cy-contract"});
	args::Flag _cy_contract_full(parser, "cy-contract", "Contract using the additional constraint cy.", {"cy-contract-full"});
	args::Flag _eps_contract(parser, "eps-contract", "Contract using eps.", {"eps-contract"});
	//args::ValueFlag<int> _nb_ub_sols(parser, "int", "Max number of solutions added by the inner-polytope algorithm (default: 50)", {"N"});
	//args::ValueFlag<double> _weight2(parser, "float", "Weight of the second objective (default: 0.01)", {"w2","weight2"});
	//args::ValueFlag<double> _min_ub_dist(parser, "float", "Min distance between two non dominated points to be considered (default: eps/10)", {"min_ub_dist"});
	args::Flag _nobisecty(parser, "nobisecty", "Do not bisect y variables.", {"no-bisecty"});
	args::ValueFlag<std::string> _upperbounding(parser, "string", "Upper bounding strategy (no|ub1|ub2) (default: no).", {"ub"});
	args::ValueFlag<double> _rh(parser, "float", "Termination criteria for the ub2 algorithm (dist < rh*ini_dist)", {"rh"});
	args::Flag _server_mode(parser, "server", "Server Mode (some options are discativated).",{"server_mode"});
	args::ValueFlag<std::string> _output_file(parser, "string", "Server Output File ", {"server_out"});
	args::ValueFlag<std::string> _instructions_file(parser, "string", "Server Instructions File", {"server_in"});
	args::ValueFlag<std::string> _input_file(parser, "string", "Loading file", {"input_file"});
	args::ValueFlag<std::string> _demo(parser, "string", "Demo file", {"demo"});

	args::ValueFlag<int> _port(parser, "int", "Port for connection with the API", {"port"});

	args::Flag verbose(parser, "verbose", "Verbose output. Shows the dominance-free set of solutions obtained by the solver.",{'v',"verbose"});
	args::Flag _trace(parser, "trace", "Activate trace. Updates of loup/uplo are printed while minimizing.", {"trace"});
	args::Flag _plot(parser, "plot", "Save a file to be plotted by plot.py.", {"plot"});
	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");

	try
	{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

  	string file = filename.Get();



	int port= 0;

	if(_port){
			port = _port.Get();
			printf("port: %i\n", port);
	}


	// aqui tira el error

  //original system: 2 objective functions and constraints
	System ext_sys(file.c_str());



  //for accing the cy envelope constraint
	SystemFactory fac2;

	Variable w;
	Variable a;

	fac2.add_var(w);
	fac2.add_var(a);

	fac2.add_var(ext_sys.args[ext_sys.nb_var-2]);
	fac2.add_var(ext_sys.args[ext_sys.nb_var-1]);
	fac2.add_ctr(ext_sys.args[ext_sys.nb_var-2] + a * ext_sys.args[ext_sys.nb_var-1] - w = 0);

	OptimizerMOP::cy_contract_var= _cy_contract || _cy_contract_full;
	OptimizerMOP::_cy_upper= _cy_contract_full;

	System *_ext_sys;
	if(OptimizerMOP::cy_contract_var)	_ext_sys=new System(ext_sys, System(fac2));
	else _ext_sys =new System(ext_sys);




	//cout << *_ext_sys << endl;

	string filtering = (_filtering)? _filtering.Get() : "acidhc4";
	string linearrelaxation= (_linear_relax)? _linear_relax.Get() : "compo";
	string bisection= (_bisector)? _bisector.Get() : "largestfirst";
	string strategy= (_strategy)? _strategy.Get() : "NDSdist";
	double eps= (_eps)? _eps.Get() : 0.01 ;
	double eps_rpm= (_eps_rpm)? _eps_rpm.Get() : 0.01 ;
	double rel_eps= (_server_mode && !_eps_r)? 0.01: ((_eps_r)? _eps_r.Get() : 0.0 );
	double eps_x= 1e-8 ;
	double timelimit = (_timelimit)? _timelimit.Get() : 100 ;
	double eqeps= 1.e-8;
	double rh=(_rh)? _rh.Get():0.1;
	bool _segments=(_upperbounding && _upperbounding.Get()=="ub1");
	bool _hamburger=(_upperbounding && _upperbounding.Get()=="ub2");

	OptimizerMOP::_plot = _plot;

	int nb_ub_sols = 10 ;
	OptimizerMOP::_min_ub_dist = 0.1;
	LoupFinderMOP::_weight2 = 0.01 ;
	bool no_bisect_y  = _nobisecty;
	OptimizerMOP::_eps_contract = _eps_contract;
	if(_server_mode) OptimizerMOP::_eps_contract = false;

	if(bisection=="largestfirst_noy"){
		bisection="largestfirst";
		no_bisect_y=true;
	}

	RNG::srand(0);

	cout << "Instance: " << argv[1] << endl;
	cout << "Filtering: " << filtering << endl;
	cout << "Linear Relax: " << linearrelaxation << endl;
	cout << "Bisector: " << bisection << endl;
	cout << "Strategy: " << strategy << endl;
	cout << "eps: " << eps << endl;
	cout << "rel_eps: " << rel_eps << endl;
	cout << "eps_x: " << eps_x << endl;
	cout << "nb_ub_sols: " << nb_ub_sols << endl;
	cout << "min_ub_dist: " << OptimizerMOP::_min_ub_dist << endl;
	cout << "plot: " <<  ((OptimizerMOP::_plot)? "yes":"no") << endl;
	cout << "weight f2: " << LoupFinderMOP::_weight2 << endl;
	cout << "bisect y?: " << ((no_bisect_y)? "no":"yes") << endl;
	cout << "cy_contract?: " << ((OptimizerMOP::cy_contract_var)? "yes":"no") << endl;
	cout << "eps_contract?: " << ((OptimizerMOP::_eps_contract)? "yes":"no") << endl;
	cout << "segments?: " << ((_segments)? "yes":"no") << endl;
	cout << "hamburger?: " << ((_hamburger)? "yes":"no") << endl;


	SystemFactory fac;

	Array<const ExprNode> symbs;
	Array<const ExprSymbol> _x1x2;
	for(int i=0; i<ext_sys.args.size()-2; i++ ){
		const ExprSymbol& x=ExprSymbol::new_(ext_sys.args[i].name);
		fac.add_var(x);
		symbs.add(x);

	}

	for(int j=2; j<ext_sys.nb_ctr; j++ ){
		const ExprNode& e = ext_sys.ctrs[j].f.expr();
		Array<const ExprSymbol> _x1x2;
		for(int i=0;i<ext_sys.args.size()-2;i++) _x1x2.add(ext_sys.ctrs[j].f.args()[i]);

		const ExprNode& new_e = ExprCopy().copy(_x1x2, symbs, e);
		ExprCtr cc(new_e, ext_sys.ctrs[j].op);
		fac.add_ctr(cc);
		//delete &new_e;

	}


	System sys(fac);
	for(int i=0; i<sys.nb_var; i++ )
		sys.box[i] = ext_sys.box[i];


	//cout << sys << endl;

	IntervalVector box = ext_sys.box.mid();
	box[sys.nb_var]=0;
	box[sys.nb_var+1]=0;

	LoupFinderMOP finder(sys, ext_sys.ctrs[0].f, ext_sys.ctrs[1].f, 1e-8, nb_ub_sols);

	CellBufferOptim* buffer;
	if(strategy=="OC3")
	  buffer = new CellSet<OC3>;
	else if(strategy=="OC4")
	  buffer = new CellSet<OC4>;
	else if(strategy=="weighted_sum")
	  buffer = new CellSet<weighted_sum>;
	else if(strategy=="NDSdist")
	  buffer = new DistanceSortedCellBufferMOP;




	// Build the bisection heuristic
	// --------------------------

	Bsc * bs;

	Vector p(ext_sys.nb_var, eps_x);
	if(no_bisect_y){
		p[ext_sys.nb_var-2]=POS_INFINITY;
		p[ext_sys.nb_var-1]=POS_INFINITY;
	}

	if (bisection=="roundrobin")
	  bs = new RoundRobin (0);
	else if (bisection== "largestfirst"){
          bs= new LargestFirst(p);
	}else if (bisection=="smearsum")
	  bs = new SmearSum(ext_sys,p);
	else if (bisection=="smearmax")
	  bs = new SmearMax(ext_sys,p);
	else if (bisection=="smearsumrel")
	  bs = new SmearSumRelative(ext_sys,p);
	else if (bisection=="smearmaxrel")
	  bs = new SmearMaxRelative(ext_sys,p);
	//else if (bisection=="lsmear")
	//  bs = new LSmear(ext_sys,p);
	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}





	// The contractors

	// the first contractor called
	CtcHC4 hc4(_ext_sys->ctrs,0.01,true);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(_ext_sys->ctrs,0.1,true);
	// hc4 inside xnewton loop
	CtcHC4 hc44xn (_ext_sys->ctrs,0.01,false);

	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4")
	Ctc3BCid c3bcidhc4(hc44cid);
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4"
	CtcCompo hc43bcidhc4 (hc4, c3bcidhc4);

	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	CtcAcid acidhc4(*_ext_sys,hc44cid,true);
	// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4"
	CtcCompo hc4acidhc4 (hc4, acidhc4);



	Ctc* ctc;
	if (filtering == "hc4")
	  ctc= &hc4;
	else if
	  (filtering =="acidhc4")
	  ctc= &hc4acidhc4;
	else if
	  (filtering =="3bcidhc4")
	  ctc= &hc43bcidhc4;
	else {cout << filtering <<  " is not an implemented  contraction  mode "  << endl; return -1;}

	Linearizer* lr;
	if (linearrelaxation=="art")
	  lr= new LinearizerCombo(*_ext_sys,LinearizerCombo::ART);
	else if  (linearrelaxation=="compo")
	  lr= new LinearizerCombo(*_ext_sys,LinearizerCombo::COMPO);
	else if (linearrelaxation=="xn")
	  lr= new LinearizerXTaylor (*_ext_sys, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM_OPP);

	CtcFixPoint* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          {
		cxn_poly = new CtcPolytopeHull(*lr);

		cxn_compo =new CtcCompo(*cxn_poly, hc44xn);
		cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
	  }
	//  the actual contractor  ctc + linear relaxation
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          ctcxn= new CtcCompo  (*ctc, *cxn);
	else
	  ctcxn = ctc;



	// the optimizer : the same precision goalprec is used as relative and absolute precision
	OptimizerMOP* o;
	if(!_server_mode){
		o = new OptimizerMOP(sys.nb_var,ext_sys.ctrs[0].f,ext_sys.ctrs[1].f, *ctcxn,*bs,*buffer,finder,
					(_hamburger)?  OptimizerMOP::HAMBURGER: (_segments)? OptimizerMOP::SEGMENTS:OptimizerMOP::POINTS,
					OptimizerMOP::MIDPOINT,	eps, rel_eps);
		if(strategy=="NDSdist")
			dynamic_cast<DistanceSortedCellBufferMOP*>(buffer)->set(o->ndsH);
			// the trace
			o->trace=(_trace)? _trace.Get() : false;
			// the allowed time for search
			o->timeout=timelimit;
			o->optimize(ext_sys.box);
			o->report();


	}else{
		OptimizerMOP_I* o;
		o = new OptimizerMOP_I(sys.nb_var,ext_sys.ctrs[0].f,ext_sys.ctrs[1].f, *ctcxn,*bs,*buffer,finder,
							(_hamburger)?  OptimizerMOP::HAMBURGER: (_segments)? OptimizerMOP::SEGMENTS:OptimizerMOP::POINTS,
							OptimizerMOP::MIDPOINT);
		if(strategy=="NDSdist")
			dynamic_cast<DistanceSortedCellBufferMOP*>(buffer)->set(o->ndsH);


		

		//Definicion de variables para el servidor
		int server_fd, new_socket, valread;
	    struct sockaddr_in address;
	    int opt = 1;
	    int addrlen = sizeof(address);
	    char bufferserver[1024] = {0};
	    string instruction = "run" ;

		// Creating socket file descriptor
	    if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) == 0)
	    {
	        perror("socket failed");
	        exit(EXIT_FAILURE);
	    }

	    // Forcefully attaching socket to the port 8080
	    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT,
	                                                  &opt, sizeof(opt)))
	    {
	        perror("setsockopt");
	        exit(EXIT_FAILURE);
	    }

	    address.sin_family = AF_INET;
	    address.sin_addr.s_addr = INADDR_ANY;
	    address.sin_port = htons( port );

	    // Forcefully attaching socket to the port 8080
	    if (bind(server_fd, (struct sockaddr *)&address,
	                                 sizeof(address))<0)
	    {
	        perror("bind failed");
	        exit(EXIT_FAILURE);
	    }
	    if (listen(server_fd, 3) < 0)
	    {
	        perror("listen");
	        exit(EXIT_FAILURE);
	    }
			string line;
			string mensaje;
			std::string respuesta;

			int data;

			IntervalVector focus = o->load(ext_sys.box);
			double ini_eps = focus.size() / 100.0;
			double iters=10;

			double eps; Vector ref(2);
			//o->run(iters, ini_eps);
			while(instruction != "fns"){
				close(new_socket);

				if ((new_socket = accept(server_fd, (struct sockaddr *)&address,
								(socklen_t*)&addrlen))<0)
				{
					perror("accept");
					exit(EXIT_FAILURE);
				}

				// Traduccion de instruccion
				valread = read( new_socket , bufferserver, 1024);
				std::string resp (bufferserver);
				cout << "message:" << resp << endl;
				stringstream message;
				message << resp;

				int max_nb_changes = 10;

        		string instruction;
				message >> instruction;

				// En el caso de que se solicitan los datos
				if( instruction == "glw"){
					list<vector<double>> lowerList;
					double eps;
					message >> eps;
					string newlower = "";
					o->lower_envelope_tostring(lowerList, eps);
					
					for (auto it = lowerList.begin(); it != lowerList.end(); it++){ 
						newlower = "/" + to_string((*it)[0]) + "," +to_string((*it)[1]); 
						char response [newlower.size()];
						strcpy(response, newlower.c_str());
						send(new_socket , response , strlen(response) , 0 );
					}
					char response [1024];
					strcpy(response, "/fs");
					send(new_socket, response, strlen(response), 0);
	
				}
				else if( instruction == "gup"){
					list<vector<double>> upperList; 
					double eps;
					message >> eps;
					cout << "eps: " << eps << endl;
					string newUpper = "";
					o->upper_envelope_tostring(upperList, eps);
					
					for (auto it = upperList.begin(); it != upperList.end(); it++){ 
						newUpper = "/" + to_string((*it)[0]) + "," +to_string((*it)[1]); 
						char response [newUpper.size()];
						strcpy(response, newUpper.c_str());
						send(new_socket , response , strlen(response) , 0 );
					}
					char response [1024];
					strcpy(response, "/fs");
					send(new_socket, response, strlen(response), 0);
				}

				// En el caso de que quieren el espacio de búsqueda
				else if( instruction == "zoo"){
					cout << "entro al zoo" << endl;
					char response [1024];
					string t = mensaje;
					istringstream iss(t);
					string word;
					int x; int y;
					message >> x;
					message >> y;
					message >> eps;
					cout << "x: " << x << endl;
					cout << "y: " << y << endl;
					cout << "eps: " << eps << endl;

					ref[0] = x;
					ref[1] = y;

					cout << "ref: " << ref << endl;
					o->update_refpoint(ref, eps);
					strcpy(response, "Puntos de referencia actualizados");
					send(new_socket , response , strlen(response) , 0 );

				}

				else if( instruction == "run"){
					char response [1024];


					message >> iters;
					message >> eps;

					cout << "iters: " << iters << endl;
					cout << "eps: " << eps << endl;

					OptimizerMOP_I::IStatus status = o->run(iters, eps);

					float precision = o->current_precision;
					cout << "precision: " << precision << endl;
					cout << "status: " << status << endl;


					strcpy(response, "Run executed\n");
					send(new_socket , response , strlen(response) , 0 );


				}else if(instruction == "plot"){
					o->plot();
				}
			}
		}
	}catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}


	return 0;


}