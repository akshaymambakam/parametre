#include <ppl.hh>
#include "I_Polyhedron.hpp"

using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;
using namespace std;

vector<I_Polyhedron> itre_union(vector<I_Polyhedron> itre1, vector<I_Polyhedron> itre2);
vector<I_Polyhedron> itre_intersection(vector<I_Polyhedron> &itre1, vector<I_Polyhedron> &itre2);
vector<I_Polyhedron> itre_concatenation(vector<I_Polyhedron> &itre1, vector<I_Polyhedron> &itre2);
vector<I_Polyhedron>& itre_duration_restriction(vector<I_Polyhedron> &itre, Linear_Expression &A, Linear_Expression &B);
// vector<I_Polyhedron> itre_KleenePlus(vector<I_Polyhedron> &plyIn);
vector<I_Polyhedron> itre_KleenePlusSquaring(vector<I_Polyhedron> &plyR);
void itre_print(vector<I_Polyhedron> &itre);