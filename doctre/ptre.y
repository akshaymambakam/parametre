%{
#include "helperFunctions.hpp"
#include "ptre.hpp"
#include "itre.hpp"

extern int yylex();
extern int yyparse();
extern FILE* yyin;
int max_dot=0;
int inputType = 0;
double pi_comp_time=0;
double seqcomp_time=0;
double intersect_time=0;
double pconcat_time=0;
double pintersect_time=0;
double isempty_time=0;
double contained_time=0;

vector<int> timeVec;
vector<int> eventVec;
vector<vector<string>> valueVec;
vector<vector<C_Polyhedron>> bZoneL,bNegZoneL;
vector<pair<int,int>> labelIntervals;

vector<Constraint> paramS;

vector<Linear_Expression> lestack;

vector<vector<C_Polyhedron>> ptreStack;

vector<vector<I_Polyhedron>> itreStack;

vector<C_Polyhedron> labelZones;

void yyerror(const char* s);
%}

%union {
	int ival;
	float fval;
}

%token<ival> T_INT
%token<ival> P_STRING X_STRING Y_STRING Z_STRING
%token<fval> T_FLOAT
%token<ival> GEQ GE LEQ LE T_MIN T_MAX T_DIFF
%token T_PLUS K_PLUS T_MINUS T_NOT T_MULTIPLY T_DIVIDE T_LEFT T_RIGHT
%token T_SLEFT T_SRIGHT
%token T_PRE T_POST T_EPS
%token<ival> T_RISE T_FALL
%token T_AND T_OR T_SEQCOMP
%left T_PLUS T_MINUS
%left T_MULTIPLY T_DIVIDE

%type<ival> expression compop edgeop

%start full_form;

%%

full_form: constraints ptre;

constraints: constraints constraint
      | T_LEFT constraint T_RIGHT
      | constraint
      ;

constraint: expression GEQ expression {
        Linear_Expression le2=lestack.back();
        lestack.pop_back();
        Linear_Expression le1=lestack.back();
        lestack.pop_back();
        paramS.push_back(le1>=le2);
      } 
      | expression LEQ expression{
        Linear_Expression le2=lestack.back();
        lestack.pop_back();
        Linear_Expression le1=lestack.back();
        lestack.pop_back();
        paramS.push_back(le1<=le2);
      }
      ;

ptre: T_EPS {
            if(inputType == 0 || inputType == 2){
                Variable t0(0),t1(1);
                vector<C_Polyhedron> pzRes;
                Constraint_System zoneTemp;
                for(int i=0;i<paramS.size();i++){
                    zoneTemp.insert(paramS[i]);
                }
                zoneTemp.insert(t1-t0==0);
                pzRes.push_back(C_Polyhedron(zoneTemp));
                ptreStack.push_back(pzRes);
            }else{
                vector<I_Polyhedron> izRes = epsToIntervalPoly(paramS, timeVec, eventVec);
                itreStack.push_back(izRes);
            }
            
      }
      | ptre T_PRE ptre {
            vector<C_Polyhedron> pz2 = ptreStack.back();
            ptreStack.pop_back();
            vector<C_Polyhedron> pz1 = ptreStack.back();
            ptreStack.pop_back();
            vector<C_Polyhedron> pzRes = ptre_precondition(pz1,pz2);
            ptreStack.push_back(pzRes);
      }
      | ptre T_POST ptre {
            vector<C_Polyhedron> pz2 = ptreStack.back();
            ptreStack.pop_back();
            vector<C_Polyhedron> pz1 = ptreStack.back();
            ptreStack.pop_back();
            vector<C_Polyhedron> pzRes = ptre_postcondition(pz1,pz2);
            ptreStack.push_back(pzRes);     
      }
      | T_RISE X_STRING{
            vector<C_Polyhedron> ptreRet;
            ptreRet = boolToEdgeZone(paramS, timeVec, valueVec[$2], 1);
            ptreStack.push_back(ptreRet);
      }
      | T_FALL X_STRING{
            vector<C_Polyhedron> ptreRet;
            ptreRet = boolToEdgeZone(paramS, timeVec, valueVec[$2], 0);
            ptreStack.push_back(ptreRet);
      }
      | ptre K_PLUS {
            vector<C_Polyhedron> pz;
            vector<I_Polyhedron> iz;
            if(inputType == 0 || inputType == 2){
                pz = ptreStack.back();
                ptreStack.pop_back();
            }else{
                iz = itreStack.back();
                itreStack.pop_back();
            }
            
            if(inputType == 0 || inputType == 2){
                vector<C_Polyhedron> pzRes = ptre_KleenePlusSquaring(pz);
                ptreStack.push_back(pzRes);
            }else{
                vector<I_Polyhedron> izRes = itre_KleenePlusSquaring(iz);
                itreStack.push_back(izRes);
            }
            
      }
      | T_LEFT ptre T_RIGHT {

      }
      | ptre T_SEQCOMP ptre{
            vector<C_Polyhedron> pz1, pz2;
            vector<I_Polyhedron> iz1, iz2;
            if(inputType == 0 || inputType == 2){
                pz2 = ptreStack.back();
                ptreStack.pop_back();
                pz1 = ptreStack.back();
                ptreStack.pop_back();
            }else{
                iz2 = itreStack.back();
                itreStack.pop_back();
                iz1 = itreStack.back();
                itreStack.pop_back();
            }

            if(inputType == 0 || inputType == 2){
                vector<C_Polyhedron> pzRes = ptre_concatenation(pz1,pz2);
                ptreStack.push_back(pzRes);
            }else{
                vector<I_Polyhedron> izRes = itre_concatenation(iz1, iz2);
                itreStack.push_back(izRes);
            }
      }
      | ptre T_AND ptre{
            vector<C_Polyhedron> pz1, pz2;
            vector<I_Polyhedron> iz1, iz2;
            if(inputType == 0 || inputType == 2){
                pz2 = ptreStack.back();
                ptreStack.pop_back();
                pz1 = ptreStack.back();
                ptreStack.pop_back();
            }else{
                iz2 = itreStack.back();
                itreStack.pop_back();
                iz1 = itreStack.back();
                itreStack.pop_back();
            }

            if(inputType == 0 || inputType == 2){
                vector<C_Polyhedron> pzRes = ptre_intersection(pz1,pz2);
                ptreStack.push_back(pzRes);
            }else{
                vector<I_Polyhedron> izRes = itre_intersection(iz1, iz2);
                itreStack.push_back(izRes);
            }
      }
      | ptre T_OR ptre{
            vector<C_Polyhedron> pz1, pz2;
            vector<I_Polyhedron> iz1, iz2;
            if(inputType == 0 || inputType == 2){
                pz2 = ptreStack.back();
                ptreStack.pop_back();
                pz1 = ptreStack.back();
                ptreStack.pop_back();
            }else{
                iz2 = itreStack.back();
                itreStack.pop_back();
                iz1 = itreStack.back();
                itreStack.pop_back();
            }

            if(inputType == 0 || inputType == 2){
                vector<C_Polyhedron> pzRes = ptre_union(pz1,pz2);
                ptreStack.push_back(pzRes);
            }else{
                vector<I_Polyhedron> izRes = itre_union(iz1, iz2);
                itreStack.push_back(izRes);
            }
      }
      | porv
      | ptre T_SLEFT expression expression T_SRIGHT{
            Variable t1(0);
            Variable t2(1);
            vector<C_Polyhedron> pz;
            vector<I_Polyhedron> iz;
            if(inputType == 0 || inputType == 2){
                pz = ptreStack.back();
                ptreStack.pop_back();
            }else{
                iz = itreStack.back();
                itreStack.pop_back();
            }
            Linear_Expression le2 = lestack.back();
            lestack.pop_back();
            Linear_Expression le1 = lestack.back();
            lestack.pop_back();

            if(inputType == 0 || inputType == 2){
                Constraint A(le1<=t2-t1);
                Constraint B(t2-t1<=le2);
                vector<C_Polyhedron> pzRes = ptre_duration_restriction(pz, A, B);
                ptreStack.push_back(pzRes);
            }else{
                vector<I_Polyhedron> izRes = itre_duration_restriction(iz, le1, le2);
                itreStack.push_back(izRes);
            }
            
      }
      | X_STRING {
            vector<C_Polyhedron> ptreRet;
            ptreRet = boolToZone(paramS, timeVec, valueVec[$1], 1);
            ptreStack.push_back(ptreRet);
      }
      | T_NOT X_STRING {
            vector<C_Polyhedron> ptreRet;
            ptreRet = boolToZone(paramS, timeVec, valueVec[$2], 0);
            ptreStack.push_back(ptreRet);
      }
      | Y_STRING {
            /*
            vector<C_Polyhedron> ptreRet;
            ptreRet = eventToZone(paramS, timeVec, eventVec, $1);
            ptreStack.push_back(ptreRet);
            */
            vector<I_Polyhedron> itreRet;
            itreRet = eventToIntervalPoly(paramS, timeVec, eventVec, $1);
            itreStack.push_back(itreRet);
      }
      | Z_STRING {
            vector<C_Polyhedron> ptreIn, ptreRet;
            ptreIn = bZoneL[$1];

            Constraint_System zoneTemp;
            for(int i=0;i<paramS.size();i++){
                zoneTemp.insert(paramS[i]);
            }

            for(int i=0;i<ptreIn.size();i++){
                C_Polyhedron pt(zoneTemp);
                pt.add_constraints(ptreIn[i].constraints());
                ptreRet.push_back(pt);
            }
            ptreStack.push_back(ptreRet);
      }
      | T_NOT Z_STRING {
            vector<C_Polyhedron> ptreIn, ptreRet;
            ptreIn = bNegZoneL[$2];

            Constraint_System zoneTemp;
            for(int i=0;i<paramS.size();i++){
                zoneTemp.insert(paramS[i]);
            }

            for(int i=0;i<ptreIn.size();i++){
                C_Polyhedron pt(zoneTemp);
                pt.add_constraints(ptreIn[i].constraints());
                ptreRet.push_back(pt);
            }
            ptreStack.push_back(ptreRet);     
      }
      ;
compop: GEQ| LEQ;
edgeop: T_RISE | T_FALL;

porv:   expression LEQ T_INT X_STRING LEQ expression{
            vector<C_Polyhedron> ptreRet;
            if($1 == 0 && $6 == 0){
                vector<string> fvec = valueVec[$4];
                vector<string> fvecBool;
                Linear_Expression le2 = lestack.back();
                lestack.pop_back();
                Linear_Expression le1 = lestack.back();
                lestack.pop_back(); 
                for(int i=0;i<fvec.size();i++){
                    Constraint_System zt;
                    int exp10 = exponentOf10(placesAfterDot(fvec[i]));
                    int fvectemp = removeDot(fvec[i]);
                    zt.insert(exp10*le1 <= $3*fvectemp);
                    zt.insert($3*fvectemp <= le2*exp10);
                    C_Polyhedron cpt(zt);
                    if(cpt.is_universe()){
                        fvecBool.push_back("1");
                    }else if(cpt.is_empty()){
                        fvecBool.push_back("0");
                    }else{
                        cerr<<"Unforseen error in porv"<<endl;
                    }
                }
                ptreRet = boolToZone(paramS, timeVec, fvecBool, 1);
            }else{
                // Incomplete
            }
            ptreStack.push_back(ptreRet);
        }
        | T_INT X_STRING compop expression 
        {
            vector<string> fvec = valueVec[$2];
            vector<C_Polyhedron> ptreRet;
            if($4==0){
                vector<string> fvecBool;
                for(int i=0;i<fvec.size();i++){
                    Constraint_System zt;
                    int exp10 = exponentOf10(placesAfterDot(fvec[i]));
                    int fvectemp = removeDot(fvec[i]);
                    if($3==1){
                        zt.insert($1*fvectemp <= exp10*lestack.back());
                        C_Polyhedron cpt(zt);
                        if(cpt.is_universe()){
                            fvecBool.push_back("1");
                        }else if(cpt.is_empty()){
                            fvecBool.push_back("0");
                        }else{
                            cerr<<"Unforseen error in porv"<<endl;
                        }
                    }else if($3==3){
                        zt.insert($1*fvectemp >= exp10*lestack.back());
                        C_Polyhedron cpt(zt);
                        if(cpt.is_universe()){
                            fvecBool.push_back("1");
                        }else if(cpt.is_empty()){
                            fvecBool.push_back("0");
                        }else{
                            cerr<<"Unforseen error in porv"<<endl;
                        }
                    }
                }
                ptreRet = boolToZone(paramS, timeVec, fvecBool, 1);
            }else{
                vector<C_Polyhedron> ptreTemp;
                // Take the subsignal
                if(labelIntervals.size()==1){
                    vector<int> ntv;
                    vector<string> nfv;
                    compute_subsignal(labelIntervals[0], timeVec, fvec, ntv, nfv);
                    ptreTemp = porvZones(ntv, nfv, paramS, $1, lestack.back(), 0, ntv.size()-2, $3);
                }else{
                    ptreTemp = porvZones(timeVec, fvec, paramS, $1, lestack.back(), 0, timeVec.size()-2, $3);    
                }
                
                // TODO: Intersection with triangles corresponding to labelling
                if(labelIntervals.size()!=0){
                    if(labelZones.size()==0){
                        labelZones = labellingToZones(paramS, labelIntervals);
                    }
                    ptreRet = ptre_intersection(labelZones, ptreTemp);
                }else{
                    ptreRet = ptreTemp;
                }
                // TODO: End
            }
            ptreStack.push_back(ptreRet);
            lestack.pop_back();
        }
        | T_INT T_MAX X_STRING GEQ expression
        {
            vector<string> fvec = valueVec[$3];
            vector<C_Polyhedron> ptreRet;
            if($5 == 0){
                vector<string> fvecBool;
                for(int i=0;i<fvec.size();i++){
                    Constraint_System zt;
                    int exp10 = exponentOf10(placesAfterDot(fvec[i]));
                    int fvectemp = removeDot(fvec[i]);

                    zt.insert($1*fvectemp <= exp10*lestack.back());
                    C_Polyhedron cpt(zt);
                    if(cpt.is_universe()){
                        fvecBool.push_back("1");
                    }else if(cpt.is_empty()){
                        fvecBool.push_back("0");
                    }else{
                        cerr<<"Unforseen error in porv"<<endl;
                    }
                }
                ptreRet = aggrToZone(paramS, timeVec, fvecBool);
            }else{
                // incomplete
            }
            ptreStack.push_back(ptreRet);
            lestack.pop_back();

            /*
            vector<string> fvec = valueVec[$3];
            vector<C_Polyhedron> ptreRet;

            ptreRet = aggrZones(timeVec, fvec, paramS, $1, lestack.back(), 0, timeVec.size()-2, $2);

            ptreStack.push_back(ptreRet);
            lestack.pop_back();
            */
        }
        | T_INT T_MIN X_STRING LEQ expression
        {
            vector<string> fvec = valueVec[$3];
            vector<C_Polyhedron> ptreRet;
            if($5 == 0){
                vector<string> fvecBool;
                for(int i=0;i<fvec.size();i++){
                    Constraint_System zt;
                    int exp10 = exponentOf10(placesAfterDot(fvec[i]));
                    int fvectemp = removeDot(fvec[i]);

                    zt.insert($1*fvectemp >= exp10*lestack.back());
                    C_Polyhedron cpt(zt);
                    if(cpt.is_universe()){
                        fvecBool.push_back("1");
                    }else if(cpt.is_empty()){
                        fvecBool.push_back("0");
                    }else{
                        cerr<<"Unforseen error in porv"<<endl;
                    }
                }
                ptreRet = aggrToZone(paramS, timeVec, fvecBool);
            }else{
                // incomplete
            }
            ptreStack.push_back(ptreRet);
            lestack.pop_back();

            /*
            vector<string> fvec = valueVec[$3];
            vector<C_Polyhedron> ptreRet;

            ptreRet = aggrZones(timeVec, fvec, paramS, $1, lestack.back(), 0, timeVec.size()-2, $2);

            ptreStack.push_back(ptreRet);
            lestack.pop_back();
            */
        }
        | T_INT T_MAX X_STRING LEQ expression
        {
            vector<string> fvec = valueVec[$3];
            vector<C_Polyhedron> ptreRet;
            if($5 == 0){
                vector<string> fvecBool;
                for(int i=0;i<fvec.size();i++){
                    Constraint_System zt;
                    int exp10 = exponentOf10(placesAfterDot(fvec[i]));
                    int fvectemp = removeDot(fvec[i]);

                    zt.insert($1*fvectemp <= exp10*lestack.back());
                    C_Polyhedron cpt(zt);
                    if(cpt.is_universe()){
                        fvecBool.push_back("1");
                    }else if(cpt.is_empty()){
                        fvecBool.push_back("0");
                    }else{
                        cerr<<"Unforseen error in porv"<<endl;
                    }
                }
                ptreRet = boolToZone(paramS, timeVec, fvecBool, 1);
            }else{
                // incomplete
            }
            ptreStack.push_back(ptreRet);
            lestack.pop_back();
        }
        | T_INT T_MIN X_STRING GEQ expression
        {
            vector<string> fvec = valueVec[$3];
            vector<C_Polyhedron> ptreRet;
            if($5 == 0){
                vector<string> fvecBool;
                for(int i=0;i<fvec.size();i++){
                    Constraint_System zt;
                    int exp10 = exponentOf10(placesAfterDot(fvec[i]));
                    int fvectemp = removeDot(fvec[i]);

                    zt.insert($1*fvectemp >= exp10*lestack.back());
                    C_Polyhedron cpt(zt);
                    if(cpt.is_universe()){
                        fvecBool.push_back("1");
                    }else if(cpt.is_empty()){
                        fvecBool.push_back("0");
                    }else{
                        cerr<<"Unforseen error in porv"<<endl;
                    }
                }
                ptreRet = boolToZone(paramS, timeVec, fvecBool, 1);
            }else{
                // incomplete
            }
            ptreStack.push_back(ptreRet);
            lestack.pop_back();
        }
        | edgeop T_INT X_STRING compop expression 
        {
            vector<string> fvec = valueVec[$3]; //TODO: Changed from valueVec[$2] to valueVec[$3]?
            vector<C_Polyhedron> ptreRet;
            if($5==0){
                vector<string> fvecBool;
                for(int i=0;i<fvec.size();i++){
                    Constraint_System zt;
                    int exp10 = exponentOf10(placesAfterDot(fvec[i]));
                    int fvectemp = removeDot(fvec[i]);
                    if($4==1){
                        zt.insert($2*fvectemp <= exp10*lestack.back());
                        C_Polyhedron cpt(zt);
                        if(cpt.is_universe()){
                            fvecBool.push_back("1");
                        }else if(cpt.is_empty()){
                            fvecBool.push_back("0");
                        }else{
                            cerr<<"Unforseen error in porv"<<endl;
                        }
                    }else if($4==3){
                        zt.insert($2*fvectemp >= exp10*lestack.back());
                        C_Polyhedron cpt(zt);
                        if(cpt.is_universe()){
                            fvecBool.push_back("1");
                        }else if(cpt.is_empty()){
                            fvecBool.push_back("0");
                        }else{
                            cerr<<"Unforseen error in porv"<<endl;
                        }
                    }
                }
                ptreRet = boolToEdgeZone(paramS, timeVec, fvecBool, $1);
            }else{
                cout<<"Error in edge porv."<<endl;
                exit(0);
            }
            ptreStack.push_back(ptreRet);
            lestack.pop_back();
        }
        | T_INT T_DIFF X_STRING LEQ expression
        {
            vector<string> fvec = valueVec[$3];
            vector<C_Polyhedron> ptreRet;

            ptreRet = diffLeqZones(timeVec, fvec, paramS, $1, lestack.back());

            ptreStack.push_back(ptreRet);
            lestack.pop_back();
        }
        | T_INT T_DIFF X_STRING GEQ expression
        {
            vector<string> fvec = valueVec[$3];
            vector<C_Polyhedron> ptreRet;

            ptreRet = diffGeqZones(timeVec, fvec, paramS, $1, lestack.back());

            ptreStack.push_back(ptreRet);
            lestack.pop_back();
        }
        | T_LEFT porv T_RIGHT
;

expression: 
      T_INT				
      {
            lestack.push_back(Linear_Expression($1));
            $$=0;
      }
	  | expression T_PLUS expression
      {
            $$ = $1 || $3;
            Linear_Expression le1=lestack.back();
            lestack.pop_back();
            Linear_Expression le2=lestack.back();
            lestack.pop_back();
            lestack.push_back(le1+le2);
      }
	  | expression T_MINUS expression	
      {
            $$ = $1 || $3;
            Linear_Expression le1=lestack.back(); 
            lestack.pop_back();
            Linear_Expression le2=lestack.back();
            lestack.pop_back();
            lestack.push_back(le2-le1);
      }
	  | T_INT T_MULTIPLY expression 
      {     
            $$ = 1;
            Linear_Expression le1=lestack.back(); 
            lestack.pop_back(); 
            lestack.push_back($1*le1);
      }
	  | T_LEFT expression T_RIGHT
      {
            $$ = $2;
      }
	  | P_STRING 
      {
            $$ = 1;
            lestack.push_back(Linear_Expression(Variable($1+3)));
      }
;
%%

void print_polyhedron_mat(C_Polyhedron ply_in, int dim1, int dim2, int dim3){
    C_Polyhedron ply = ply_in;
    // project on dim1,dim2,dim3 space
    int num_dim = ply.space_dimension();
    for(int i=0;i<num_dim;i++){
        if(i!=dim1 && i!=dim2 && i!=dim3){
            ply.unconstrain(Variable(i));
        }
    }
    // cout<<ply<<endl;
    Constraint_System cs = ply.minimized_constraints();
    for(Constraint_System_const_iterator csit=cs.begin();csit!=cs.end();csit++){
        cout<<Linear_Expression::zero()-csit->coefficient(Variable(dim1))<<",";
        cout<<Linear_Expression::zero()-csit->coefficient(Variable(dim2))<<",";
        cout<<Linear_Expression::zero()-csit->coefficient(Variable(dim3))<<",";
        cout<<Linear_Expression::zero()+csit->inhomogeneous_term();
        cout<<":";
    }
    cout<<endl;
}

vector<C_Polyhedron> compute_params_separately(char *forfile, vector<pair<int,int>> allLabels){
    vector<C_Polyhedron> ret;
    for(int i=0;i<allLabels.size();i++){
        labelIntervals={};
        labelZones={};
        paramS={};
        lestack={};
        ptreStack={};
        labelIntervals.push_back(allLabels[i]);

        cout<<"i:"<<i<<endl;
        
        yyin = fopen(forfile,"r");
        do {
            yyparse();
        } while(!feof(yyin));
        fclose(yyin);
        
        // Compute the projection on labelling and intersect
        vector<C_Polyhedron> pres = ptreStack[0];
        vector<C_Polyhedron> projected = ptre_label_project(pres, paramS, labelIntervals);
        if(i==0){
            ret = projected;
        }else{
            ret = poly_list_intersection(ret, projected);
        }
    }
    return ret;
}

int main(int argc, char** argv) {
    // TODO: Clean here
    vector<pair<int,int>> allLabels;

    if(argc<5){
        cout<<"Improper arguments"<<endl;
        exit(0);
    }
    if(argc==6){
        // TODO: Clean here
        // readLabelFile(argv[5], labelIntervals);
        readLabelFile(argv[5], allLabels);
    }

    inputType = atoi(argv[4]);

    if(inputType == 0){
        readDataFile(argv[2], timeVec, valueVec);
    }else if(inputType == 1){
        readEventFile(argv[2], timeVec, eventVec);
    }else if(inputType == 2){
        readBoolZones(argv[2], bZoneL, bNegZoneL);
    }

    int numParams = atoi(argv[3]);

    if(argc==6){
        vector<C_Polyhedron> param_sep = compute_params_separately(argv[1], allLabels);
        ptre_print(param_sep);
        exit(0);
    }
    
    yyin = fopen(argv[1],"r");
	do {
		yyparse();
	} while(!feof(yyin));

    if(inputType == 0 || inputType == 2){
        vector<C_Polyhedron> pres = ptreStack[0];
        ptre_print(pres);
        if(argc==6){
            vector<C_Polyhedron> pspace = ptre_label_project(pres, paramS, labelIntervals);
            ptre_print(pspace);
        }
    }else{
        itre_print(itreStack[0]);
    }
	return 0;
}

void yyerror(const char* s) {
	fprintf(stderr, "Parse error: %s\n", s);
	exit(1);
}