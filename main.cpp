#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "grammar.hpp"
#include "codegen.hpp"

using namespace std;

start_node * grammar::start = NULL;

int main (int argc, char **argv) {
	string outfile ("--out-file");
	string out_name ("out.cu");
	string datatype ("--datatype");
	string data_type ("double");
	string unroll ("--unroll");
	map<string, int> unroll_decls;
	vector<string> iters;
	string dist_rhs ("--distribute-rhs");
	bool distribute_rhs = true;
	string split_accs ("--split-size");
	int split_size = 1; 
	string top_sort ("--topo-sort");
	bool topo_sort = false;
	string gen_fma ("--gen-fma");
	bool fma = true;
	string print_intrinsics ("--print-intrinsics");
	bool intrinsics = true;
	string vect_size ("--vect-size");
	int vsize = 256;
	string sort_func ("--sort-function");
	string sort_metric ("b6");

	if (DEBUG) printf ("filename : %s\n", argv[1]);
	for (int i=2; i<argc; i++) {
		if (outfile.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing output file name, using default out.cu\n");
			else 
				out_name = argv[++i];
		}
		if (datatype.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing datatype, using default double\n");
			else 
				data_type = argv[++i];
		}
		if (dist_rhs.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing distribute rhs value, using default true\n");
			else {
				string tmp = argv[++i];
				distribute_rhs = (tmp.compare ("true") == 0) ? true : false;
			} 
		}
		if (split_accs.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing split size value, using default 1\n");
			else {
				split_size = atoi (argv[++i]);
			} 
		}
		if (top_sort.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing topological sort value, using default false\n");
			else {
				string tmp = argv[++i];
				topo_sort = (tmp.compare ("true") == 0) ? true : false;
			}
		}
		if (gen_fma.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing generate fma value, using default true\n");
			else {
				string tmp = argv[++i];
				fma = (tmp.compare ("true") == 0) ? true : false;
			}
		}
		if (print_intrinsics.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing printing instrincis value, using default true\n");
			else {
				string tmp = argv[++i];
				intrinsics = (tmp.compare ("true") == 0) ? true : false;
			}
		}
		if (vect_size.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing vector size, using default 256\n");
			else 
				vsize = atoi (argv[++i]);
		}
		if (sort_func.compare (argv[i]) == 0) {
			if (i == argc-1) 
				printf ("Missing metric sort value, using default b6(Jacobi-like)\n");
			else {
				sort_metric = argv[++i];
			}
		}
		if (unroll.compare (argv[i]) == 0) {
			if (i != argc-1) {
				string tmp = argv[++i];
				while (tmp.find (",") != string::npos) {
					size_t pos = tmp.find (",");
					string uf = tmp.substr(0,pos);
					size_t vpos = uf.find ("=");
					string dimension = uf.substr(0,vpos);
					int val = atoi(uf.substr(vpos+1).c_str());
					unroll_decls[dimension] = val;
					iters.push_back (dimension);
					tmp = tmp.substr(pos+1);		
				}
                                size_t vpos = tmp.find ("=");
                                string dimension = tmp.substr(0,vpos);
                                int val = atoi(tmp.substr(vpos+1).c_str());
                                unroll_decls[dimension] = val;
				iters.push_back (dimension);
			}
		}
	}
	if (DEBUG) printf ("output file : %s\n", out_name.c_str());
	FILE *in = fopen (argv[1], "r");
	string slc_name = "slc-" + out_name;
	string slc_acc_name = "slc-acc-" + out_name;
	ofstream slc_out (slc_name.c_str(), ofstream::out);
	ofstream slc_acc_out (slc_acc_name.c_str(), ofstream::out);
	ofstream reorder_out (out_name.c_str(), ofstream::out);

	grammar::set_input (in);
	grammar::parse ();
	codegen *sp_gen = new codegen (grammar::start);
	stringstream reorder, slc, slc_acc;
	DATA_TYPE gdata_type = DOUBLE; 
	if (data_type.compare ("float") == 0) 
		gdata_type = FLOAT;
	reverse (iters.begin(), iters.end());
	sp_gen->generate_code (reorder, slc, slc_acc, unroll_decls, iters, gdata_type, vsize, split_size, distribute_rhs, topo_sort, fma, intrinsics, sort_metric);
	slc_out << slc.rdbuf ();
	slc_acc_out << slc_acc.rdbuf ();
	reorder_out << reorder.rdbuf ();
	slc_out.close ();
	slc_acc_out.close ();
	reorder_out.close ();
	fclose (in);
	return 0;
}
