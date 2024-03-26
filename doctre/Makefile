all: ptre

ptre.tab.c ptre.tab.h:	ptre.y
	bison -t -v -d ptre.y

lex.yy.c: ptre.l ptre.tab.h
	flex ptre.l

ptre: lex.yy.c ptre.tab.c ptre.tab.h ptre.cpp helperFunctions.cpp A_Polyhedron.cpp pzone.cpp I_Polyhedron.cpp itre.cpp
	g++ -o ptre ptre.tab.c lex.yy.c helperFunctions.cpp A_Polyhedron.cpp ptre.cpp pzone.cpp I_Polyhedron.cpp itre.cpp -lppl -lgmp

clean:
	rm ptre ptre.tab.c lex.yy.c ptre.tab.h ptre.output