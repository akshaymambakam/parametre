#include <ppl.hh>
#include "A_Polyhedron.hpp"

#include <ctime>

using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;
using namespace std;

vector<C_Polyhedron> ptre_union(vector<C_Polyhedron> ptre1, vector<C_Polyhedron> ptre2);
vector<C_Polyhedron> ptre_intersection(vector<C_Polyhedron> &ptre1, vector<C_Polyhedron> &ptre2);
vector<C_Polyhedron> ptre_concatenation(vector<C_Polyhedron> &ptre1, vector<C_Polyhedron> &ptre2);
vector<C_Polyhedron> ptre_duration_restriction(vector<C_Polyhedron> &ptre, Constraint A, Constraint B);
vector<C_Polyhedron> ptre_KleenePlus(vector<C_Polyhedron> &plyIn);
vector<C_Polyhedron> ptre_KleenePlusSquaring(vector<C_Polyhedron> &plyR);
vector<C_Polyhedron> ptre_postcondition(vector<C_Polyhedron> ptre1, vector<C_Polyhedron> ptre2);
vector<C_Polyhedron> ptre_precondition(vector<C_Polyhedron> ptre1, vector<C_Polyhedron> ptre2);
vector<C_Polyhedron> ptre_filter(vector<C_Polyhedron> ptre);
void ptre_print(vector<C_Polyhedron> ptre);
void ptre_print(vector<A_Polyhedron> ptre);

// Examples
void example_intersection();