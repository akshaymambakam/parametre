#include <ppl.hh>

#define NUM_BOUNDS (6)

using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;
using namespace std;

class pzone{
  Coefficient c1_n,c1_d;
  Coefficient c2_n,c2_d;
  Coefficient c3_n,c3_d;
  Coefficient c4_n,c4_d;
  Coefficient c5_n,c5_d;
  Coefficient c6_n,c6_d;

  public:
  	vector<vector<Linear_Expression>> cbounds;
    C_Polyhedron cons;
    pzone(vector<Linear_Expression> a,vector<Linear_Expression> b,vector<Linear_Expression> c,
		vector<Linear_Expression> d,vector<Linear_Expression> e,vector<Linear_Expression> f,
		C_Polyhedron pcons);
    pzone(vector<vector<Linear_Expression>> cbounds,C_Polyhedron pcons);
    void print();
    void refineB(int iOrig, int i1, char sn, int i2, char ineq);
    void refineC(int iOrig, int i1, char sn, int i2, char ineq);
    void refineTauto(int iOrig, int i1, char sn, int i2, char ineq);
    void refine();
    int is_empty();
    int is_emptyBackup();
    pzone intersection(pzone &p);
    pzone sequential_composition(pzone &p);
    pzone duration_restriction(Linear_Expression A, Linear_Expression B);
    void getEarlyEnd(Coefficient &coe_n, Coefficient &coe_d);
    void getEarlyBegin(Coefficient &coe_n, Coefficient &coe_d);
    void getLateEnd(Coefficient &coe_n, Coefficient &coe_d);
    void getLateBegin(Coefficient &coe_n, Coefficient &coe_d);
    C_Polyhedron to_polyhedron();
};

pzone ptre_poly2pzone(Constraint_System cs);