#ifndef __EXPRNODE_HPP__
#define __EXPRNODE_HPP__
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <cmath>
#include <iostream>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include "utils.hpp" 

class expr_node {
	protected:
		ETYPE expr_type;
		DATA_TYPE type=BOOL;
		std::string name="";
		bool nested=false, gen_fma=true, print_intrinsics=true;
		int vsize = 256;
	public:
		ETYPE get_expr_type (void);
		void set_expr_type (ETYPE);
		virtual void set_type (DATA_TYPE) {}
		void set_name (std::string);
		void set_name (char *);
		std::string get_name (void);	
		DATA_TYPE get_type (void);
		void set_nested (void);
		bool is_nested (void);
		virtual bool is_data_type (DATA_TYPE gdata_type);
		virtual bool is_data_type (void);
		virtual bool is_shiftvec_type (DATA_TYPE gdata_type);
		virtual bool is_shiftvec_type (void);
		virtual bool simple_nondecomposable_expr (void);
		virtual bool is_id_type (DATA_TYPE gdata_type);
		virtual bool is_id_type (void); 
		virtual void print_node (std::stringstream &) {}
		virtual void print_node (std::stringstream &, std::vector<std::string> &, std::vector<std::string> &, bool, bool) {}
		virtual void print_node (std::map<std::string, std::string> &, std::stringstream &) {} 
		virtual void decompose_node (std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<expr_node*> &, expr_node *, STMT_OP, int &, DATA_TYPE, bool &, bool &, bool) {}
		virtual void create_labels (std::map<std::string, expr_node*> &) {}
		virtual void print_initializations (std::stringstream &, std::vector<std::string> &, std::vector<std::string>, bool, bool) {}
		virtual void create_labels (std::map<std::string, int> &, std::map<std::string, expr_node*> &, bool) {}
		virtual void stringify_accesses (std::vector<std::string> &) {}
		virtual void stringify_accesses (std::vector<std::string> &, std::string &) {}
		virtual void array_access_info (std::string &, int &) {}
		virtual void gather_participating_labels (std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string>) {}
		virtual expr_node *unroll_expr (std::string, int, std::vector<std::string>, std::map<std::string,int> &, bool) {}
		virtual expr_node *deep_copy (void) {}
		virtual void set_vect_size (int);
		virtual void set_gen_fma (bool);	
		virtual void set_print_intrinsics (bool);
		virtual void set_codegen_parameters (int, bool, bool); 
};

inline ETYPE expr_node::get_expr_type (void) {
	return expr_type;
}

inline void expr_node::set_expr_type (ETYPE type) {
	expr_type = type;
}

inline void expr_node::set_nested (void) {
	nested = true;
}

inline bool expr_node::is_nested (void) {
	return nested;
}

inline void expr_node::set_vect_size (int v) {
	vsize = v;
}

inline void expr_node::set_gen_fma (bool g) {
	gen_fma = g;
}

inline void expr_node::set_print_intrinsics (bool p) {
	print_intrinsics = p;
}

inline void expr_node::set_codegen_parameters (int v, bool g, bool p) {
	vsize = v;
	gen_fma = g;
	print_intrinsics = p;
}

template <typename T>
class datatype_node : public expr_node {
	private:
		T value;
	public:
		datatype_node (T, DATA_TYPE);
		datatype_node (T, DATA_TYPE, int, bool, bool);
		T get_value (void);
		void set_type (DATA_TYPE);
		void print_node (std::stringstream &);
		void print_node (std::stringstream &, std::vector<std::string> &, std::vector<std::string> &, bool, bool);
		void print_node (std::map<std::string, std::string> &, std::stringstream &);
		void decompose_node (std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<expr_node*> &, expr_node *, STMT_OP, int &, DATA_TYPE, bool &, bool &, bool);
		void array_access_info (std::string &, int &);
		void stringify_accesses (std::vector<std::string> &, std::string &);
		expr_node *unroll_expr (std::string, int, std::vector<std::string>, std::map<std::string,int> &, bool) {return this;}
		expr_node *deep_copy (void) {return this;}
		void set_codegen_parameters (int, bool, bool);
};

template <typename T>
inline datatype_node<T>::datatype_node (T val, DATA_TYPE dtype) {
	value = val;
	type = dtype;
	expr_type = T_DATATYPE;
}

template <typename T>
inline datatype_node<T>::datatype_node (T val, DATA_TYPE dtype, int v, bool g, bool p) {
	value = val;
	type = dtype;
	expr_type = T_DATATYPE;
	vsize = v;
	gen_fma = g;
	print_intrinsics = p;
}

template <typename T>
inline T datatype_node<T>::get_value (void) {
	return value;
}			

template <typename T>
inline void datatype_node<T>::set_type (DATA_TYPE dtype) {
	type = dtype;
}

template<typename T>
inline void datatype_node<T>::set_codegen_parameters (int v, bool g, bool p) {
        vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

template<typename T>
inline void datatype_node<T>::print_node (std::stringstream &out) {
	std::stringstream val;
	val << value;
	if (type == FLOAT) {
		if (val.str().find (".") == std::string::npos && 
                        val.str().find("e") == std::string::npos && val.str().find("E") == std::string::npos)
			val << ".0";
		val << "f";
	}
	if (type == DOUBLE) {
		if (val.str().find (".") == std::string::npos &&
                        val.str().find("e") == std::string::npos && val.str().find("E") == std::string::npos)
			val << ".0";
	}
	if (print_intrinsics) {
		if (type == FLOAT) 
			out << "_mm" << vsize << "_set1_ps (" << val.str() << ")";
		else if (type == DOUBLE)
			out << "_mm" << vsize << "_set1_pd (" << val.str() << ")";
		else
			out << "_mm" << vsize << "_set1_epi32 (" << val.str() << ")";
	}
	else out << val.str ();
}

template<typename T>
inline void datatype_node<T>::print_node (std::stringstream &out, std::vector<std::string> &initialized_labels, std::vector<std::string> & iters, bool perform_load, bool is_lhs) {
	std::stringstream val;
	val << value;
	if (type == FLOAT) {
		if (val.str().find (".") == std::string::npos &&
			val.str().find("e") == std::string::npos && val.str().find("E") == std::string::npos)
			val << ".0";
		val << "f";
	}
	if (type == DOUBLE) {
		if (val.str().find (".") == std::string::npos && 
			val.str().find("e") == std::string::npos && val.str().find("E") == std::string::npos)
			val << ".0";
	}
	if (print_intrinsics) {
		if (type == FLOAT)
			out << "_mm" << vsize << "_set1_ps (" << val.str() << ")";
		else if (type == DOUBLE)
			out << "_mm" << vsize << "_set1_pd (" << val.str() << ")";
		else
			out << "_mm" << vsize << "_set1_epi32 (" << val.str() << ")";
	}
	else out << val.str ();
}

template<typename T>
inline void datatype_node<T>::print_node (std::map<std::string, std::string> &reg_map, std::stringstream &out) {
	std::stringstream val;
	val << value;
	if (type == FLOAT) {
		if (val.str().find (".") == std::string::npos &&
                        val.str().find("e") == std::string::npos && val.str().find("E") == std::string::npos)
			val << ".0";
		val << "f";
	}
	if (type == DOUBLE) {
		if (val.str().find (".") == std::string::npos &&
                        val.str().find("e") == std::string::npos && val.str().find("E") == std::string::npos)
			val << ".0";
	}
	if (print_intrinsics) {
		if (type == FLOAT)
			out << "_mm" << vsize << "_set1_ps (" << val.str() << ")";
		else if (type == DOUBLE)
			out << "_mm" << vsize << "_set1_pd (" << val.str() << ")";
		else
			out << "_mm" << vsize << "_set1_epi32 (" << val.str() << ")";
	}
	else out << val.str ();
}

template<typename T>
inline void datatype_node<T>::decompose_node (std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &tstmt, std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &init, std::vector<expr_node*> &temp_vars, expr_node *alhs, STMT_OP cur_op, int &id, DATA_TYPE gdata_type, bool &local_assigned, bool &global_assigned, bool flip) {
	cur_op = get_cur_op (cur_op, flip); 
	// Create a node as alhs += 4.2f;
	if (!global_assigned) {
		tstmt.push_back (std::make_tuple (alhs, this, ST_EQ));
		local_assigned = true;
		global_assigned = true;
	}
	else tstmt.push_back (std::make_tuple (alhs, this, cur_op));
	// Infer types
	alhs->set_type (gdata_type);
}

template<typename T>
inline void datatype_node<T>::array_access_info (std::string &id, int &offset) {
	if (DEBUG) assert (type == INT && "Array access has non-int offset (array_access_info)");
	offset += value;
}

template<typename T>
inline void datatype_node<T>::stringify_accesses (std::vector<std::string> &labels, std::string &expr_label) {
	std::stringstream val;
	print_node (val);
	expr_label = expr_label + val.str ();
}

class id_node : public expr_node {
	private:
		std::string label;
	public:
		id_node (char *);
		id_node (std::string);
		id_node (char *, DATA_TYPE, bool, std::string);
		id_node (std::string, DATA_TYPE, bool, std::string);
		id_node (char *, int, bool, bool);
		id_node (std::string, int, bool, bool);
		id_node (char *, DATA_TYPE, int, bool, bool);
		id_node (std::string, DATA_TYPE, int, bool, bool);
		id_node (char *, DATA_TYPE, bool, std::string, int, bool, bool);
		id_node (std::string, DATA_TYPE, bool, std::string, int, bool, bool);
		void set_type (DATA_TYPE);
		void print_node (std::stringstream &);
		void print_node (std::stringstream &, std::vector<std::string> &, std::vector<std::string> &, bool, bool);
		void print_initializations (std::stringstream &, std::vector<std::string> &, std::vector<std::string>, bool, bool);
		void print_node (std::map<std::string, std::string> &, std::stringstream &);
		void decompose_node (std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<expr_node*> &, expr_node *, STMT_OP, int &, DATA_TYPE, bool &, bool &, bool);
		void create_labels (std::map<std::string, expr_node*> &);
		void create_labels (std::map<std::string, int> &, std::map<std::string, expr_node*> &, bool);
		void stringify_accesses (std::vector<std::string> &); 
		void stringify_accesses (std::vector<std::string> &, std::string &);
		void array_access_info (std::string &, int &);
		void gather_participating_labels (std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string>);
		expr_node *unroll_expr (std::string, int, std::vector<std::string>, std::map<std::string,int> &, bool);
		expr_node *deep_copy (void);
		void set_label (std::string s) {label = s;}
		void set_codegen_parameters (int, bool, bool);	
};

inline id_node::id_node (char *s) {
	name = std::string (s);
	label = std::string (s);
	expr_type = T_ID;
}

inline id_node::id_node (char *s, DATA_TYPE t, bool nest, std::string l) {
	name = std::string (s);
	expr_type = T_ID;
	type = t;
	nested = nest;
	label = l;
}

inline id_node::id_node (std::string s) {
	name = s;
	label = s;
	expr_type = T_ID;
}

inline id_node::id_node (std::string s, DATA_TYPE t, bool nest, std::string l) {
	name = s;
	expr_type = T_ID;
	type = t;
	nested = nest;
	label = l;
}

inline id_node::id_node (char *s, int v, bool g, bool p) {
        name = std::string (s);
        label = std::string (s);
        expr_type = T_ID;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline id_node::id_node (char *s, DATA_TYPE t, int v, bool g, bool p) {
        name = std::string (s);
        label = std::string (s);
        expr_type = T_ID;
	type = t;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline id_node::id_node (char *s, DATA_TYPE t, bool nest, std::string l, int v, bool g, bool p) {
        name = std::string (s);
        expr_type = T_ID;
        type = t;
        nested = nest;
        label = l;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline id_node::id_node (std::string s, int v, bool g, bool p) {
        name = s;
        label = s;
        expr_type = T_ID;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline id_node::id_node (std::string s, DATA_TYPE t, int v, bool g, bool p) {
        name = s;
        label = s;
        expr_type = T_ID;
	type = t;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline id_node::id_node (std::string s, DATA_TYPE t, bool nest, std::string l, int v, bool g, bool p) {
        name = s;
        expr_type = T_ID;
        type = t;
        nested = nest;
        label = l;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline void id_node::set_type (DATA_TYPE dtype) {
	type = dtype;
}

inline void id_node::set_codegen_parameters (int v, bool g, bool p) {
        vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

class uminus_node : public expr_node {
	private:
		expr_node *base_expr;
	public:
		uminus_node (expr_node *);
		uminus_node (expr_node *, char *, DATA_TYPE, bool);
		uminus_node (expr_node *, std::string, DATA_TYPE, bool);
		uminus_node (expr_node *, int, bool, bool);
		uminus_node (expr_node *, char *, DATA_TYPE, bool, int, bool, bool);
		uminus_node (expr_node *, std::string, DATA_TYPE, bool, int, bool, bool);
		expr_node *get_base_expr (void);
		void set_type (DATA_TYPE);
		bool is_data_type (DATA_TYPE gdata_type);
		bool is_data_type (void);
		bool simple_nondecomposable_expr (void);
		bool is_shiftvec_type (DATA_TYPE gdata_type);
		bool is_shiftvec_type (void);
		bool is_id_type (DATA_TYPE gdata_type);
		bool is_id_type (void);
		void print_node (std::stringstream &);
		void print_node (std::stringstream &, std::vector<std::string> &, std::vector<std::string> &, bool, bool);
		void print_initializations (std::stringstream &, std::vector<std::string> &, std::vector<std::string>, bool, bool);
		void print_node (std::map<std::string, std::string> &, std::stringstream &);
		void decompose_node (std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<expr_node*> &, expr_node *, STMT_OP, int &, DATA_TYPE, bool &, bool &, bool);
		void create_labels (std::map<std::string, expr_node*> &);
		void create_labels (std::map<std::string, int> &, std::map<std::string, expr_node*> &, bool);
		void stringify_accesses (std::vector<std::string> &, std::string &);
		void array_access_info (std::string &, int &);
		void gather_participating_labels (std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string>);
		expr_node *unroll_expr (std::string, int, std::vector<std::string>, std::map<std::string,int> &, bool);
		expr_node *deep_copy (void);
		void set_codegen_parameters (int, bool, bool);
};

inline uminus_node::uminus_node (expr_node *expr) {
	base_expr = expr;
	expr_type = T_UMINUS;
}

inline uminus_node::uminus_node (expr_node *expr, char *s, DATA_TYPE t, bool nest) {
	base_expr = expr;
	expr_type = T_UMINUS;
	name = std::string (s);
	type = t;
	nested = nest;
}

inline uminus_node::uminus_node (expr_node *expr, std::string s, DATA_TYPE t, bool nest) {
	base_expr = expr;
	expr_type = T_UMINUS;
	name = s;
	type = t;
	nested = nest;
}

inline uminus_node::uminus_node (expr_node *expr, int v, bool g, bool p) {
	base_expr = expr;
	expr_type = T_UMINUS;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline uminus_node::uminus_node (expr_node *expr, char *s, DATA_TYPE t, bool nest, int v, bool g, bool p) {
	base_expr = expr;
	expr_type = T_UMINUS;
	name = std::string (s);
	type = t;
	nested = nest;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline uminus_node::uminus_node (expr_node *expr, std::string s, DATA_TYPE t, bool nest, int v, bool g, bool p) {
	base_expr = expr;
	expr_type = T_UMINUS;
	name = s;
	type = t;
	nested = nest;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline void uminus_node::set_type (DATA_TYPE dtype) {
	base_expr->set_type (dtype);
	type = dtype;
}

inline void uminus_node::set_codegen_parameters (int v, bool g, bool p) {
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
	base_expr->set_codegen_parameters (v, g, p);
}

class binary_node : public expr_node {
	private:
		OP_TYPE op;
		expr_node *lhs, *rhs;
	public:
		binary_node (OP_TYPE, expr_node *, expr_node *);
		binary_node (OP_TYPE, expr_node *, expr_node *, char *, DATA_TYPE, bool);
		binary_node (OP_TYPE, expr_node *, expr_node *, std::string, DATA_TYPE, bool);
		binary_node (OP_TYPE, expr_node *, expr_node *, int, bool, bool);
		binary_node (OP_TYPE, expr_node *, expr_node *, char *, DATA_TYPE, bool, int, bool, bool);
		binary_node (OP_TYPE, expr_node *, expr_node *, std::string, DATA_TYPE, bool, int, bool, bool);
		OP_TYPE get_operator (void);
		expr_node *get_rhs (void);
		expr_node *get_lhs (void);
		void set_type (DATA_TYPE);
		bool is_data_type (DATA_TYPE gdata_type);
		bool is_data_type (void);
		bool is_shiftvec_type (DATA_TYPE gdata_type);
		bool is_shiftvec_type (void);
		bool is_id_type (DATA_TYPE gdata_type);
		bool is_id_type (void);
		bool simple_nondecomposable_expr (void);
		void print_node (std::stringstream &);
		void print_node (std::stringstream &, std::vector<std::string> &, std::vector<std::string> &, bool, bool);
		void print_initializations (std::stringstream &, std::vector<std::string> &, std::vector<std::string>, bool, bool);
		void print_node (std::map<std::string, std::string> &, std::stringstream &);
		void decompose_node (std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<expr_node*> &, expr_node *, STMT_OP, int &, DATA_TYPE, bool &, bool &, bool);
		void create_labels (std::map<std::string, expr_node*> &);
		void create_labels (std::map<std::string, int> &, std::map<std::string, expr_node*> &, bool);
		void stringify_accesses (std::vector<std::string> &, std::string &);
		void array_access_info (std::string &, int &);
		void gather_participating_labels (std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string>);
		expr_node *unroll_expr (std::string, int, std::vector<std::string>, std::map<std::string,int> &, bool);
		expr_node *deep_copy (void);
		void set_codegen_parameters (int, bool, bool);
};

inline binary_node::binary_node (OP_TYPE a_type, expr_node *a_lhs, expr_node *a_rhs) {
	op = a_type;
	lhs = a_lhs;
	rhs = a_rhs;
	expr_type = T_BINARY;
}

inline binary_node::binary_node (OP_TYPE a_type, expr_node *a_lhs, expr_node *a_rhs, char *s, DATA_TYPE t, bool nest) {
	op = a_type;
	lhs = a_lhs;
	rhs = a_rhs;
	expr_type = T_BINARY;
	name = std::string (s);
	type = t;
	nested = nest;
}

inline binary_node::binary_node (OP_TYPE a_type, expr_node *a_lhs, expr_node *a_rhs, std::string s, DATA_TYPE t, bool nest) {
	op = a_type;
	lhs = a_lhs;
	rhs = a_rhs;
	expr_type = T_BINARY;
	name = s;
	type = t;
	nested = nest;
}

inline binary_node::binary_node (OP_TYPE a_type, expr_node *a_lhs, expr_node *a_rhs, int v, bool g, bool p) {
	op = a_type;
	lhs = a_lhs;
	rhs = a_rhs;
	expr_type = T_BINARY;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline binary_node::binary_node (OP_TYPE a_type, expr_node *a_lhs, expr_node *a_rhs, char *s, DATA_TYPE t, bool nest, int v, bool g, bool p) {
	op = a_type;
	lhs = a_lhs;
	rhs = a_rhs;
	expr_type = T_BINARY;
	name = std::string (s);
	type = t;
	nested = nest;
	vsize =v;
        gen_fma = g;
        print_intrinsics = p;
}

inline binary_node::binary_node (OP_TYPE a_type, expr_node *a_lhs, expr_node *a_rhs, std::string s, DATA_TYPE t, bool nest, int v, bool g, bool p) {
	op = a_type;
	lhs = a_lhs;
	rhs = a_rhs;
	expr_type = T_BINARY;
	name = s;
	type = t;
	nested = nest;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline OP_TYPE binary_node::get_operator (void) {
	return op;
}

inline expr_node *binary_node::get_lhs (void) {
	return lhs;
}

inline expr_node *binary_node::get_rhs (void) {
	return rhs;
}

inline void binary_node::set_type (DATA_TYPE dtype) {
	lhs->set_type (dtype);
	rhs->set_type (dtype);
	type = dtype;
}

class shiftvec_node : public expr_node {
	private:
		std::vector<expr_node*> indices;
		std::string label;
	public:
		shiftvec_node ();
		shiftvec_node (char *);
		shiftvec_node (char *, DATA_TYPE, bool, std::string); 
		shiftvec_node (std::string);
		shiftvec_node (std::string, DATA_TYPE, bool, std::string);
		shiftvec_node (int, bool, bool);
		shiftvec_node (char *, int, bool, bool);
		shiftvec_node (char *, DATA_TYPE, bool, std::string, int, bool, bool); 
		shiftvec_node (std::string, int, bool, bool);
		shiftvec_node (std::string, DATA_TYPE, bool, std::string, int, bool, bool);
		void set_type (DATA_TYPE); 
		std::vector<expr_node*>& get_indices (void);
		void push_index (expr_node *);
		void print_node (std::stringstream &);
		void print_node (std::stringstream &, std::vector<std::string> &, std::vector<std::string> &, bool, bool);
		void print_initializations (std::stringstream &, std::vector<std::string> &, std::vector<std::string>, bool, bool);
		void print_node (std::map<std::string, std::string> &, std::stringstream &);
		int get_index_size (void);
		void decompose_node (std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<expr_node*> &, expr_node *, STMT_OP, int &, DATA_TYPE, bool &, bool &, bool);
		void create_labels (std::map<std::string, expr_node*> &);
		void create_labels (std::map<std::string, int> &, std::map<std::string, expr_node*> &, bool);
		void stringify_accesses (std::vector<std::string> &);
		void stringify_accesses (std::vector<std::string> &, std::string &);
		void gather_participating_labels (std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string>);
		void lexical_index_offsets (std::vector<int> &);
		expr_node *unroll_expr (std::string, int, std::vector<std::string>, std::map<std::string,int> &, bool);
		expr_node *deep_copy (void);
		void set_label (std::string s) {label = s;}
		std::string print_array (void);
};

inline shiftvec_node::shiftvec_node () {
	expr_type = T_SHIFTVEC;
}

inline shiftvec_node::shiftvec_node (char *s) {
	name = std::string (s);
	expr_type = T_SHIFTVEC;
}

inline shiftvec_node::shiftvec_node (std::string s) {
	name = s;
	expr_type = T_SHIFTVEC;
}

inline shiftvec_node::shiftvec_node (std::string s, DATA_TYPE t, bool nest, std::string l) {
	name = s;
	expr_type = T_SHIFTVEC;
	type = t;
	nested = nest;
	label = l;
}


inline shiftvec_node::shiftvec_node (char *s, DATA_TYPE t, bool nest, std::string l) {
	name = std::string (s);
	expr_type = T_SHIFTVEC;
	type = t;
	nested = nest;
	label = l;
}

inline shiftvec_node::shiftvec_node (int v, bool g, bool p) {
	expr_type = T_SHIFTVEC;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline shiftvec_node::shiftvec_node (char *s, int v, bool g, bool p) {
	name = std::string (s);
	expr_type = T_SHIFTVEC;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline shiftvec_node::shiftvec_node (char *s, DATA_TYPE t, bool nest, std::string l, int v, bool g, bool p) {
	name = std::string (s);
	expr_type = T_SHIFTVEC;
	type = t;
	nested = nest;
	label = l;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline shiftvec_node::shiftvec_node (std::string s, int v, bool g, bool p) {
	name = s;
	expr_type = T_SHIFTVEC;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline shiftvec_node::shiftvec_node (std::string s, DATA_TYPE t, bool nest, std::string l, int v, bool g, bool p) {
	name = s;
	expr_type = T_SHIFTVEC;
	type = t;
	nested = nest;
	label = l;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}


inline void shiftvec_node::set_type (DATA_TYPE dtype) {
	type = dtype;
}

inline void shiftvec_node::push_index (expr_node *node) {
	indices.push_back (node);
}

class function_node : public expr_node {
	private:
		expr_node *arg;
	public:
		function_node (char *, expr_node *);
		function_node (char *, expr_node *, DATA_TYPE, bool);
		function_node (std::string, expr_node *);
		function_node (std::string, expr_node *, DATA_TYPE, bool);
		function_node (char *, expr_node *, int, bool, bool);
		function_node (char *, expr_node *, DATA_TYPE, bool, int, bool, bool);
		function_node (std::string, expr_node *, int, bool, bool);
		function_node (std::string, expr_node *, DATA_TYPE, bool, int, bool, bool);
		bool is_data_type (DATA_TYPE gdata_type);
		bool is_data_type (void);
		void set_type (DATA_TYPE);
		bool is_shiftvec_type (DATA_TYPE gdata_type);
		bool is_id_type (DATA_TYPE gdata_type);
		void print_node (std::stringstream &);
		void print_node (std::stringstream &, std::vector<std::string> &, std::vector<std::string> &, bool, bool);
		void print_initializations (std::stringstream &, std::vector<std::string> &, std::vector<std::string>, bool, bool);
		void print_node (std::map<std::string, std::string> &, std::stringstream &);
		void decompose_node (std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<std::tuple<expr_node*, expr_node*, STMT_OP>> &, std::vector<expr_node*> &, expr_node *, STMT_OP, int &, DATA_TYPE, bool &, bool &, bool);
		void create_labels (std::map<std::string, expr_node*> &);
		void create_labels (std::map<std::string, int> &, std::map<std::string, expr_node*> &, bool);
		void stringify_accesses (std::vector<std::string> &);
		void stringify_accesses (std::vector<std::string> &, std::string &);
		void gather_participating_labels (std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string>);
		expr_node *unroll_expr (std::string, int, std::vector<std::string>, std::map<std::string,int> &, bool);
		expr_node *deep_copy (void);
		void set_codegen_parameters (int, bool, bool);
};

inline function_node::function_node (char *s, expr_node *t) {
	name = std::string (s);
	arg = t;
	expr_type = T_FUNCTION; 
}

inline function_node::function_node (char *s, expr_node *ag, DATA_TYPE t, bool nest) {
	name = std::string (s);
	arg = ag;
	expr_type = T_FUNCTION; 
	type = t;
	nested = nest;
}

inline function_node::function_node (std::string s, expr_node *t) {
	name = s;
	arg = t;
	expr_type = T_FUNCTION; 
}

inline function_node::function_node (std::string s, expr_node *ag, DATA_TYPE t, bool nest) { 
	name = s;
	arg = ag; 
	expr_type = T_FUNCTION;
	type = t;
	nested = nest;
}

inline function_node::function_node (char *s, expr_node *t, int v, bool g, bool p) {
	name = std::string (s);
	arg = t;
	expr_type = T_FUNCTION; 
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline function_node::function_node (char *s, expr_node *ag, DATA_TYPE t, bool nest, int v, bool g, bool p) {
	name = std::string (s);
	arg = ag;
	expr_type = T_FUNCTION; 
	type = t;
	nested = nest;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline function_node::function_node (std::string s, expr_node *t, int v, bool g, bool p) {
	name = s;
	arg = t;
	expr_type = T_FUNCTION; 
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline function_node::function_node (std::string s, expr_node *ag, DATA_TYPE t, bool nest, int v, bool g, bool p) { 
	name = s;
	arg = ag; 
	expr_type = T_FUNCTION;
	type = t;
	nested = nest;
	vsize = v;
        gen_fma = g;
        print_intrinsics = p;
}

inline void function_node::set_type (DATA_TYPE dtype) {
	arg->set_type (dtype);
	type = dtype;
}

#endif
