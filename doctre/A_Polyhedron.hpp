#include <ppl.hh>

using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;
using namespace std;

class A_Polyhedron{
  Coefficient c1_n,c1_d;
  Coefficient c2_n,c2_d;
  Coefficient c5_n,c5_d;
  Coefficient c6_n,c6_d;
  public:
    C_Polyhedron porig;
    A_Polyhedron(C_Polyhedron porig);
    void getEarlyEnd(Coefficient &coe_n, Coefficient &coe_d);
    void getEarlyBegin(Coefficient &coe_n, Coefficient &coe_d);
    void getLateBegin(Coefficient &coe_n, Coefficient &coe_d);
    void getLateEnd(Coefficient &coe_n, Coefficient &coe_d);
};