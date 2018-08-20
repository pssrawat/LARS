LEX=lex
YACC=yacc
CC=g++

CFLAGS=-g -O3 -std=c++11
OPTFLAGS=-DDEBUG=true -DDROP_COEFS=true -DASSOC_MULT=false -DAVAIL_EXPR_OPT=false -DRETAIN_SIMPLE_OPS=false -DSPLICE_EQUALITY=false -DPAR_LOADS=1 -DEXPLICIT_LOADS=false
INTEROPTFLAGS=-DRESTRICT_INTER_OPT=true -DSPLICE_TEMP_LABELS=false -DINTRA_TYPE_INTER_OPT=false 
REGALLOCFLAGS=-DFIRST_LEVEL=false -DSECOND_LEVEL=true -DOPERATION_VIEW=false

default : test

exprnode.o : exprnode.cpp exprnode.hpp utils.hpp sort.hpp 
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o exprnode.o -c exprnode.cpp
funcdefn.o : funcdefn.cpp funcdefn.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o funcdefn.o -c funcdefn.cpp
vardecl.o : vardecl.cpp vardecl.hpp  utils.hpp sort.hpp 
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o vardecl.o -c vardecl.cpp
codegen.o : codegen.cpp codegen.hpp utils.hpp sort.hpp 
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o codegen.o -c codegen.cpp
test : lex.yy.c y.tab.c main.cpp utils.hpp sort.hpp exprnode.o funcdefn.o vardecl.o codegen.o
	$(CC) $(CFLAGS) $(OPTFLAGS) $(REGALLOCFLAGS) -o test main.cpp exprnode.o codegen.o funcdefn.o vardecl.o lex.yy.c y.tab.c

all : lex.yy.c y.tab.c main.cpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o exprnode.o -c exprnode.cpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o funcdefn.o -c funcdefn.cpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o vardecl.o -c vardecl.cpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o codegen.o -c codegen.cpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(INTEROPTFLAGS) $(REGALLOCFLAGS) -o test main.cpp exprnode.o codegen.o funcdefn.o vardecl.o lex.yy.c y.tab.c

lex.yy.c : scanner.l
	$(LEX) scanner.l

y.tab.c : grammar.y
	$(YACC) -d grammar.y

clean:
	-@rm *.o lex.yy.* y.tab.* out.cu orig_out.cu test 2>/dev/null || true
