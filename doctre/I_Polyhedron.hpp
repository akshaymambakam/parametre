#include <ppl.hh>

using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;
using namespace std;

#ifndef IPOLY_H_INCLUDED
#define IPOLY_H_INCLUDED

class I_Polyhedron{
  public:
    int begin, end;
    C_Polyhedron phed;
    I_Polyhedron(C_Polyhedron phed, int beg, int end);
    int contains(I_Polyhedron &opoly);
};

#endif