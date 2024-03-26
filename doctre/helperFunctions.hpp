#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <ppl.hh>
#include <fstream>
#include "I_Polyhedron.hpp"

using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;
using namespace std;

int placesAfterDot(string sf);
int exponentOf10(int expo);
int removeDot(string fnum);
void readLabelFile(char *fName, vector<pair<int,int>> &labellingList);
void readEventFile(char *fName, vector<int> &timeVec, vector<int> &eventNumVec);
void readDataFile(char* fName, vector<int> &timeVec, vector<vector<string>> &valueVec);
vector<C_Polyhedron> labellingToZones(vector<Constraint> paramS, vector<pair<int,int>> &labellingList);
void readBoolZones(char *fName, vector<vector<C_Polyhedron>> &bZonelist, 
        vector<vector<C_Polyhedron>> &bNegZonelist);
void ptreUnconstrain(vector<C_Polyhedron> &ptre, int numParams);
vector<C_Polyhedron> eventToZone(vector<Constraint> paramS, vector<int> timeVec,
    vector<int> eventVec, int eventNum);
vector<I_Polyhedron> eventToIntervalPoly(vector<Constraint> paramS, vector<int> timeVec,
    vector<int> eventVec, int eventNum);
vector<I_Polyhedron> epsToIntervalPoly(vector<Constraint> paramS, vector<int> timeVec,
    vector<int> eventVec);
vector<C_Polyhedron> boolToZone(vector<Constraint> paramS, vector<int> timeVec, 
	vector<string> fvec, int prevOneIn);
vector<C_Polyhedron> aggrToZone(vector<Constraint> paramS, vector<int> timeVec, 
    vector<string> fvec);
vector<C_Polyhedron> boolToEdgeZone(vector<Constraint> paramS, vector<int> timeVec,
    vector<string> fvec, int prevOneIn);
vector<C_Polyhedron> projectParams(vector<C_Polyhedron> ptre);
vector<C_Polyhedron> porvZones(vector<int> timeVec, vector<string> fvec, vector<Constraint> paramS,
		int coeff, Linear_Expression le, int begin, int end, int mxFlag);
vector<C_Polyhedron> aggrZones(vector<int> timeVec, vector<string> fvec, vector<Constraint> paramS,
        int coeff, Linear_Expression le, int begin, int end, int mxFlag);
vector<C_Polyhedron> diffLeqZones(vector<int> timeVec, vector<string> fvec, vector<Constraint> paramS,
        int coeff, Linear_Expression le);
vector<C_Polyhedron> diffGeqZones(vector<int> timeVec, vector<string> fvec, vector<Constraint> paramS,
        int coeff, Linear_Expression le);
C_Polyhedron labelZoneCombination(C_Polyhedron plabel, C_Polyhedron matchZone, int numParams);
vector<C_Polyhedron> labelZoneList(vector<C_Polyhedron> vLabel, vector<C_Polyhedron> vZone, int numParams);

vector<C_Polyhedron> poly_list_filter(vector<C_Polyhedron> plyIn);
vector<C_Polyhedron> poly_list_intersection(vector<C_Polyhedron> &ptre1, vector<C_Polyhedron> &ptre2);
vector<C_Polyhedron> ptre_label_project(vector<C_Polyhedron> ptre, vector<Constraint> paramS,
    vector<pair<int,int>> labelIntervals);

void compute_subsignal(pair<int,int> interval, vector<int> orig_tvec, vector<string> orig_fvec, vector<int> &new_tvec,
    vector<string> &new_fvec);