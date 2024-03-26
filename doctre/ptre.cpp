#include "ptre.hpp"
#include "pzone.hpp"
#include <chrono> 
using namespace std::chrono;

extern double seqcomp_time;
extern double intersect_time;
extern double pconcat_time;
extern double pintersect_time;
extern double isempty_time;
extern double contained_time;

#define FILTER (1)
#define PLANE_SWEEP (1)

bool lessCompareBounds(Coefficient c1_n, Coefficient c1_d, Coefficient c2_n, Coefficient c2_d){
	return c1_n*c2_d<c2_n*c1_d;
}

bool geqCompareBounds(Coefficient c1_n, Coefficient c1_d, Coefficient c2_n, Coefficient c2_d){
	return c1_n*c2_d>=c2_n*c1_d;
}

struct {
    bool operator()(A_Polyhedron &a, A_Polyhedron &b) const
    {   
    	Coefficient eb1_n, eb1_d;
    	Coefficient eb2_n, eb2_d;
    	a.getEarlyBegin(eb1_n, eb1_d);
    	b.getEarlyBegin(eb2_n, eb2_d);
        return eb1_n*eb2_d-eb2_n*eb1_d < 0;
    }  
}compareEarlyBeginAug;

struct {
    bool operator()(A_Polyhedron &a, A_Polyhedron &b) const
    {   
    	Coefficient eb1_n, eb1_d;
    	Coefficient eb2_n, eb2_d;
    	a.getEarlyEnd(eb1_n, eb1_d);
    	b.getEarlyEnd(eb2_n, eb2_d);
        return eb1_n*eb2_d-eb2_n*eb1_d < 0;
    }  
}compareEarlyEndAug;

struct {
    bool operator()(A_Polyhedron &a, A_Polyhedron &b) const
    {   
    	Coefficient eb1_n, eb1_d;
    	Coefficient eb2_n, eb2_d;
    	a.getLateBegin(eb1_n, eb1_d);
    	b.getLateBegin(eb2_n, eb2_d);
        return eb1_n*eb2_d-eb2_n*eb1_d < 0;
    }  
}compareLateBeginAug;

struct {
    bool operator()(A_Polyhedron &a, A_Polyhedron &b) const
    {   
    	Coefficient eb1_n, eb1_d;
    	Coefficient eb2_n, eb2_d;
    	a.getLateEnd(eb1_n, eb1_d);
    	b.getLateEnd(eb2_n, eb2_d);
        return eb1_n*eb2_d-eb2_n*eb1_d < 0;
    }  
}compareLateEndAug;

void ptre_print(vector<C_Polyhedron> ptre){
	cout<<"ptre size:"<<ptre.size()<<endl;
	for(int i=0;i<ptre.size();i++){
		cout<<ptre[i]<<endl;
	}
	cout<<"----ptre----"<<endl;
}

void ptre_print(vector<A_Polyhedron> ptre){
    cout<<"ptre size:"<<ptre.size()<<endl;
    for(int i=0;i<ptre.size();i++){
        cout<<ptre[i].porig<<endl;
    }
    cout<<"----ptre----"<<endl;
}

void poly_duration_restriction(C_Polyhedron &ply, Constraint A, Constraint B){
    ply.add_constraint(A);
    ply.add_constraint(B);
}

C_Polyhedron poly_intersection(C_Polyhedron ply1, C_Polyhedron &ply2){
    auto start = high_resolution_clock::now();

	ply1.intersection_assign(ply2);

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    intersect_time += duration.count();

    return ply1;
}

C_Polyhedron poly_sequential_composition_fm(C_Polyhedron &ply1, C_Polyhedron &ply2){
    static int trace = 0;
    pzone pz1 = ptre_poly2pzone(ply1.minimized_constraints());
    pzone pz2 = ptre_poly2pzone(ply2.minimized_constraints());
    pzone seqcomp = pz1.sequential_composition(pz2);
    cout<<"here:"<<trace<<endl;
    trace++;
    return seqcomp.to_polyhedron();
}

C_Polyhedron poly_sequential_composition_cyl(C_Polyhedron &ply1, C_Polyhedron &ply2){
    auto start = high_resolution_clock::now();

    Constraint_System cs1 = ply1.minimized_constraints();
    Constraint_System cs2 = ply2.minimized_constraints();

    Constraint_System cssc;
    Constraint_System csnew;

    Variable t0(0),t1(1),t2(2);
    // Swap t' with t''
    for(Constraint_System_const_iterator csit=cs1.begin();csit!=cs1.end();csit++){
        Constraint ctemp=*csit;
        ctemp.swap_space_dimensions(t1,t2);
        csnew.insert(ctemp);
    }
    // Swap t with t''
    for(Constraint_System_const_iterator csit=cs2.begin();csit!=cs2.end();csit++){
        Constraint ctemp=*csit;
        ctemp.swap_space_dimensions(t0,t2);
        csnew.insert(ctemp);
    }

    C_Polyhedron plyRes(csnew);
    plyRes.unconstrain(t2);

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    seqcomp_time += duration.count();

    return plyRes;
}

C_Polyhedron poly_precondition(C_Polyhedron ply1, C_Polyhedron ply2){
    Constraint_System cs1 = ply1.minimized_constraints();

    Constraint_System cssc;
    Constraint_System csnew;

    Variable t0(0),t1(1),t2(2);
    // Swap t with t' and then t' and t''
    for(Constraint_System_const_iterator csit=cs1.begin();csit!=cs1.end();csit++){
        Constraint ctemp=*csit;
        ctemp.swap_space_dimensions(t0,t1);
        ctemp.swap_space_dimensions(t1,t2);
        csnew.insert(ctemp);
    }
    csnew.insert(t2<=t0);

    C_Polyhedron ptemp(csnew); 
    ptemp.unconstrain(t2);

    ply2.intersection_assign(ptemp);

    return ply2;
}

C_Polyhedron poly_postcondition(C_Polyhedron ply1, C_Polyhedron ply2){
    Constraint_System cs2 = ply2.minimized_constraints();

    Constraint_System cssc;
    Constraint_System csnew;

    Variable t0(0),t1(1),t2(2);

    // Swap t with t' and then t with t''
    for(Constraint_System_const_iterator csit=cs2.begin();csit!=cs2.end();csit++){
        Constraint ctemp=*csit;
        ctemp.swap_space_dimensions(t0,t1);
        ctemp.swap_space_dimensions(t0,t2);
        csnew.insert(ctemp);
    }
    csnew.insert(t2>=t1);

    C_Polyhedron ptemp(csnew);
    ptemp.unconstrain(t2);

    ply1.intersection_assign(ptemp);
    return ply1;
}

vector<C_Polyhedron> ptre_postcondition(vector<C_Polyhedron> ptre1, vector<C_Polyhedron> ptre2){
    vector<C_Polyhedron> ptre_post;
    for(int i=0;i<ptre1.size();i++){
        for(int j=0;j<ptre2.size();j++){
            C_Polyhedron polyRes = poly_postcondition(ptre1[i],ptre2[j]);
            if(!polyRes.is_empty())
                ptre_post.push_back(polyRes);
        }
    }
    return ptre_post;
}

vector<C_Polyhedron> ptre_precondition(vector<C_Polyhedron> ptre1, vector<C_Polyhedron> ptre2){
    vector<C_Polyhedron> ptre_pre;
    for(int i=0;i<ptre1.size();i++){
        for(int j=0;j<ptre2.size();j++){
            C_Polyhedron polyRes = poly_precondition(ptre1[i],ptre2[j]);
            if(!polyRes.is_empty())
                ptre_pre.push_back(polyRes);
        }
    }
    return ptre_pre;
}

#if PLANE_SWEEP

vector<C_Polyhedron> ptre_filter(vector<C_Polyhedron> plyIn){
    vector<A_Polyhedron> zs;

    for(int i=0;i<plyIn.size();i++){
        zs.push_back(A_Polyhedron(plyIn[i]));
    }

    sort(zs.begin(), zs.end(), compareEarlyBeginAug);

    vector<A_Polyhedron> active, active_temp;

    vector<C_Polyhedron> _retvalue;

    for(auto z1 : zs){

        bool already_included = std::any_of(active.begin(), active.end(), 
            [&z1](A_Polyhedron &z2){return z2.porig.contains(z1.porig);});

        if(!already_included){

            active.erase( std::remove_if(active.begin(), active.end(), 
                [&z1](A_Polyhedron &z2){return z1.porig.contains(z2.porig);}), active.end());
            active.push_back(z1);

            active_temp.clear();
            for(auto z2 : active){
                Coefficient z2_n,z2_d;
                z2.getLateBegin(z2_n,z2_d);
                Coefficient z1_n,z1_d;
                z1.getEarlyBegin(z1_n,z1_d);
                if(lessCompareBounds(z2_n,z2_d,z1_n,z1_d)){
                    if(!z2.porig.is_empty())
                        _retvalue.push_back(z2.porig);
                }
                else {
                    active_temp.push_back(z2);
                }
            }
            active = active_temp;
        }
    }
    for(auto z2 : active){
        if(!z2.porig.is_empty())
            _retvalue.push_back(z2.porig);
    }

    return _retvalue;
}

vector<C_Polyhedron> ptre_concatenation(vector<C_Polyhedron> &ptre1, vector<C_Polyhedron> &ptre2){
    auto start = high_resolution_clock::now();

    vector<A_Polyhedron> zs1,zs2;
    
    for(int i=0;i<ptre1.size();i++){
    	zs1.push_back(A_Polyhedron(ptre1[i]));
    }

    for(int i=0;i<ptre2.size();i++){
    	zs2.push_back(A_Polyhedron(ptre2[i]));
    }

    std::vector <A_Polyhedron> act_1, act_2;
    vector<C_Polyhedron> _retvalue;
  
    std::sort(zs1.begin(), zs1.end(), compareEarlyEndAug);
    std::sort(zs2.begin(), zs2.end(), compareEarlyBeginAug);
  
  
    size_t i = 0, j = 0;
    while (i < zs1.size() && j < zs2.size()) {
      Coefficient c_n, c_d;
      Coefficient cp_n, cp_d;
      zs1[i].getEarlyEnd(c_n,c_d);
      zs2[j].getEarlyBegin(cp_n,cp_d);
  
  
      if (lessCompareBounds(c_n,c_d,cp_n,cp_d)) { //  z1.emin < z2.bmin
        act_1.push_back(zs1[i]);
        act_2.erase(std::remove_if(act_2.begin(), act_2.end(), [ & ](A_Polyhedron &z2) {
          Coefficient yp_lb_n, yp_lb_d;
          z2.getLateBegin(yp_lb_n, yp_lb_d);
          return lessCompareBounds(yp_lb_n,yp_lb_d,c_n,c_d);
        }), act_2.end()); // remove if z2.bmax < z1.emin
        for (auto z2: act_2) {
          C_Polyhedron seqcomp = poly_sequential_composition_cyl(zs1[i].porig,z2.porig);
          
          auto start = high_resolution_clock::now();
          
          bool seq_empty = seqcomp.is_empty();
          if(!seq_empty){
              _retvalue.push_back(seqcomp);
          }
              
          auto stop = high_resolution_clock::now();
  
          auto duration = duration_cast<microseconds>(stop - start);
          isempty_time += duration.count();

        }
        ++i; // Dont forget to proceed
      } else {
        act_2.push_back(zs2[j]);
        act_1.erase(std::remove_if(act_1.begin(), act_1.end(), [ & ](A_Polyhedron &z1) {
          Coefficient y_n,y_d;
          z1.getLateEnd(y_n,y_d);
          return lessCompareBounds(y_n,y_d,cp_n,cp_d);
        }), act_1.end()); // remove if z1.emax < z2.bmin
        for (auto z1: act_1) {
          C_Polyhedron seqcomp = poly_sequential_composition_cyl(z1.porig,zs2[j].porig);
        
          auto start = high_resolution_clock::now();
          
          bool seq_empty = seqcomp.is_empty();
          if(!seq_empty){
              _retvalue.push_back(seqcomp);
          }
              
          auto stop = high_resolution_clock::now();
  
          auto duration = duration_cast<microseconds>(stop - start);
          isempty_time += duration.count();

        }
        ++j; // Dont forget to proceed
      }
    }
    while(i<zs1.size()){
  	  Coefficient c_n, c_d;
  	  zs1[i].getEarlyEnd(c_n,c_d);
  
  	  act_1.push_back(zs1[i]);
  	  act_2.erase(std::remove_if(act_2.begin(), act_2.end(), [ & ](A_Polyhedron &z2) {
          Coefficient yp_lb_n, yp_lb_d;
          z2.getLateBegin(yp_lb_n, yp_lb_d);
          return lessCompareBounds(yp_lb_n,yp_lb_d,c_n,c_d);
        }), act_2.end()); // remove if z2.bmax < z1.emin
	  for (auto z2: act_2) {
	    C_Polyhedron seqcomp = poly_sequential_composition_cyl(zs1[i].porig,z2.porig);

	    auto start = high_resolution_clock::now();
        
        bool seq_empty = seqcomp.is_empty();
        if(!seq_empty){
            _retvalue.push_back(seqcomp);
        }
            
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(stop - start);
        isempty_time += duration.count();

      }
  	  ++i;
    }
    while(j<zs2.size()){
  	  Coefficient cp_n, cp_d;
  	  zs2[j].getEarlyBegin(cp_n,cp_d);
  
  	  act_2.push_back(zs2[j]);
        act_1.erase(std::remove_if(act_1.begin(), act_1.end(), [ & ](A_Polyhedron &z1) {
          Coefficient y_n,y_d;
          z1.getLateEnd(y_n,y_d);
          return lessCompareBounds(y_n,y_d,cp_n,cp_d);
        }), act_1.end()); // remove if z1.emax < z2.bmin
      for (auto z1: act_1) {
        C_Polyhedron seqcomp = poly_sequential_composition_cyl(z1.porig,zs2[j].porig);

        auto start = high_resolution_clock::now();
        
        bool seq_empty = seqcomp.is_empty();
        if(!seq_empty){
            _retvalue.push_back(seqcomp);
        }
            
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(stop - start);
        isempty_time += duration.count();

      }
      ++j;
    }

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    pconcat_time += duration.count();

    #if FILTER
    return ptre_filter(_retvalue);
    #else
    return _retvalue;
    #endif
}

vector<C_Polyhedron> ptre_intersection(vector<C_Polyhedron> &ptre1, vector<C_Polyhedron> &ptre2){
    auto start = high_resolution_clock::now();

    vector<A_Polyhedron> zs1,zs2;
    
    for(int i=0;i<ptre1.size();i++){
    	zs1.push_back(A_Polyhedron(ptre1[i]));
    }

    for(int i=0;i<ptre2.size();i++){
    	zs2.push_back(A_Polyhedron(ptre2[i]));
    }

    std::vector <A_Polyhedron> act_1, act_2;
    vector<C_Polyhedron> _retvalue;
  
    std::sort(zs1.begin(), zs1.end(), compareEarlyBeginAug);
    std::sort(zs2.begin(), zs2.end(), compareEarlyBeginAug);
  
  
    size_t i = 0, j = 0;
    while (i < zs1.size() && j < zs2.size()) {
      Coefficient c_n, c_d;
      Coefficient cp_n, cp_d;
      zs1[i].getEarlyBegin(c_n,c_d);
      zs2[j].getEarlyBegin(cp_n,cp_d);
  
  
      if (lessCompareBounds(c_n,c_d,cp_n,cp_d)) { //  z1.emin < z2.bmin
        act_1.push_back(zs1[i]);
        act_2.erase(std::remove_if(act_2.begin(), act_2.end(), [ & ](A_Polyhedron &z2) {
          Coefficient yp_lb_n, yp_lb_d;
          z2.getLateBegin(yp_lb_n, yp_lb_d);
          return lessCompareBounds(yp_lb_n,yp_lb_d,c_n,c_d);
        }), act_2.end()); // remove if z2.bmax < z1.emin
        for (auto z2: act_2) {
          C_Polyhedron intersec = poly_intersection(zs1[i].porig,z2.porig);

          auto start = high_resolution_clock::now();
          
          bool inter_empty = intersec.is_empty();
          if(!inter_empty){
              _retvalue.push_back(intersec);
          }
              
          auto stop = high_resolution_clock::now();
  
          auto duration = duration_cast<microseconds>(stop - start);
          isempty_time += duration.count();

        }
        ++i; // Dont forget to proceed
      } else {
        act_2.push_back(zs2[j]);
        act_1.erase(std::remove_if(act_1.begin(), act_1.end(), [ & ](A_Polyhedron &z1) {
          Coefficient y_n,y_d;
          z1.getLateBegin(y_n,y_d);
          return lessCompareBounds(y_n,y_d,cp_n,cp_d);
        }), act_1.end()); // remove if z1.emax < z2.bmin
        for (auto z1: act_1) {
          C_Polyhedron intersec = poly_intersection(z1.porig,zs2[j].porig);

          auto start = high_resolution_clock::now();
          
          bool inter_empty = intersec.is_empty();
          if(!inter_empty){
              _retvalue.push_back(intersec);
          }
              
          auto stop = high_resolution_clock::now();
  
          auto duration = duration_cast<microseconds>(stop - start);
          isempty_time += duration.count();

        }
        ++j; // Dont forget to proceed
      }
    }
  
    while(i<zs1.size()){
  	  Coefficient c_n, c_d;
  	  zs1[i].getEarlyBegin(c_n,c_d);
  
  	  act_1.push_back(zs1[i]);
  	  act_2.erase(std::remove_if(act_2.begin(), act_2.end(), [ & ](A_Polyhedron &z2) {
          Coefficient yp_lb_n, yp_lb_d;
          z2.getLateBegin(yp_lb_n, yp_lb_d);
          return lessCompareBounds(yp_lb_n,yp_lb_d,c_n,c_d);
        }), act_2.end()); // remove if z2.bmax < z1.emin
	  for (auto z2: act_2) {
	    C_Polyhedron intersec = poly_intersection(zs1[i].porig,z2.porig);

        auto start = high_resolution_clock::now();
        
        bool inter_empty = intersec.is_empty();
        if(!inter_empty){
            _retvalue.push_back(intersec);
        }
            
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(stop - start);
        isempty_time += duration.count();

	  }
  	  ++i;
    }
  
    while(j<zs2.size()){
  	  Coefficient cp_n, cp_d;
  	  zs2[j].getEarlyBegin(cp_n,cp_d);
  
  	  act_2.push_back(zs2[j]);
        act_1.erase(std::remove_if(act_1.begin(), act_1.end(), [ & ](A_Polyhedron &z1) {
          Coefficient y_n,y_d;
          z1.getLateBegin(y_n,y_d);
          return lessCompareBounds(y_n,y_d,cp_n,cp_d);
        }), act_1.end()); // remove if z1.emax < z2.bmin
      for (auto z1: act_1) {
        C_Polyhedron intersec = poly_intersection(z1.porig,zs2[j].porig);

        auto start = high_resolution_clock::now();
        
        bool inter_empty = intersec.is_empty();
        if(!inter_empty){
            _retvalue.push_back(intersec);
        }
            
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(stop - start);
        isempty_time += duration.count();

      }
      ++j;
    }

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    pintersect_time += duration.count();

    #if FILTER
    return ptre_filter(_retvalue);
    #else
    return _retvalue;
    #endif
}

#else

vector<C_Polyhedron> ptre_filter(vector<C_Polyhedron> plyIn){
    vector<C_Polyhedron> plyFil;
    for(int i=0;i<plyIn.size();i++){
        int toAdd = 1;
        for(int j=0;j<plyIn.size();j++){
            if(i!=j && plyIn[j].strictly_contains(plyIn[i])){
                toAdd = 0;
            }else if(j<i && plyIn[j].contains(plyIn[i])){
                toAdd = 0;
            }
        }
        if(toAdd){
            plyFil.push_back(plyIn[i]);
        }
    }
    return plyFil;
}

vector<C_Polyhedron> ptre_concatenation(vector<C_Polyhedron> &ptre1, vector<C_Polyhedron> &ptre2){
    vector<C_Polyhedron> ptre_cat;
    for(int i=0;i<ptre1.size();i++){
        for(int j=0;j<ptre2.size();j++){
            C_Polyhedron ptemp =poly_sequential_composition_cyl(ptre1[i],ptre2[j]);

            auto start = high_resolution_clock::now();
        
            bool seqcomp_empty = ptemp.is_empty();
            if(!seqcomp_empty){
                ptre_cat.push_back(ptemp);
            }
            
            auto stop = high_resolution_clock::now();

            auto duration = duration_cast<microseconds>(stop - start);
            isempty_time += duration.count();

        }
    }
    #if FILTER
    return ptre_filter(ptre_cat);
    #else
    return ptre_cat;
    #endif
}

vector<C_Polyhedron> ptre_intersection(vector<C_Polyhedron> &ptre1, vector<C_Polyhedron> &ptre2){
    vector<C_Polyhedron> ptre_inter;
    for(int i=0;i<ptre1.size();i++){
        for(int j=0;j<ptre2.size();j++){
            C_Polyhedron pinter = poly_intersection(ptre1[i],ptre2[j]);

            auto start = high_resolution_clock::now();
        
            bool inter_empty = pinter.is_empty();
            if(!inter_empty){
                ptre_inter.push_back(pinter);
            }
            
            auto stop = high_resolution_clock::now();

            auto duration = duration_cast<microseconds>(stop - start);
            isempty_time += duration.count();

        }
    }
    #if FILTER
    return ptre_filter(ptre_inter);
    #else
    return ptre_inter;
    #endif
}

#endif

vector<C_Polyhedron> ptre_union(vector<C_Polyhedron> ptre1, vector<C_Polyhedron> ptre2){
    ptre1.insert(ptre1.end(), ptre2.begin(), ptre2.end());
    #if FILTER
    return ptre_filter(ptre1);
    #else
    return ptre1;
    #endif
}

vector<C_Polyhedron> ptre_duration_restriction(vector<C_Polyhedron> &ptre, Constraint A, Constraint B){
    vector<C_Polyhedron> ret_dr; 
    for(int i=0;i<ptre.size();i++){
        poly_duration_restriction(ptre[i],A,B);
        if(!ptre[i].is_empty()){
            ret_dr.push_back(ptre[i]);
        }
    }
    #if FILTER
    return ptre_filter(ret_dr);
    #else
    return ret_dr;
    #endif
}

int ptre_contained(vector<C_Polyhedron> plyNew, vector<C_Polyhedron> plyOld){
    auto start = high_resolution_clock::now();

	int all_contained = 1;
	for(int i=0;i<plyNew.size();i++){
		int contained = 0;	
		for(int j=0;j<plyOld.size();j++){
			contained = (contained || plyOld[j].contains(plyNew[i]));
		}
		all_contained = (all_contained && contained);
	}

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    contained_time += duration.count();

	return all_contained;
}

vector<C_Polyhedron> ptre_KleenePlus(vector<C_Polyhedron> &plyIn){
	vector<C_Polyhedron> plyNew=plyIn, plyOld=plyIn;
	int count =0;
	do{
		plyOld = plyNew;
		plyNew = ptre_union(ptre_concatenation(plyOld, plyIn), plyIn);
		cout<<"1,2:"<<plyOld.size()<<","<<plyNew.size()<<endl;
		count++;
	}while(!ptre_contained(plyNew, plyOld));

	cout<<"Kleene count:"<<count<<endl;
	return plyNew;
}

vector<C_Polyhedron> ptre_KleenePlusSquaring(vector<C_Polyhedron> &plyR){
	vector<C_Polyhedron> plyQ,plyP;
	plyQ = plyR;
	plyP = ptre_concatenation(plyR, plyR);
	int count =0;
	while(!ptre_contained(plyP, plyQ)){
		plyQ = ptre_union(ptre_union(plyP,plyQ), ptre_concatenation(plyP,plyQ));
		plyP = ptre_concatenation(plyP,plyP);
		cout<<"1,2:"<<plyP.size()<<","<<plyQ.size()<<endl;
		count++;
	}
	cout<<"Kleene count:"<<count<<endl;
	return plyQ;
}