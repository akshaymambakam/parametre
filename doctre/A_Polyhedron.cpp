#include "A_Polyhedron.hpp"
#include <chrono> 
using namespace std::chrono;

extern double pi_comp_time;

A_Polyhedron::A_Polyhedron(C_Polyhedron porig){
	Variable t0(0),t1(1);
	this->porig = porig;
	bool minim,maxim;

	auto start = high_resolution_clock::now();

	porig.minimize(t0,c6_n,c6_d,minim);
	porig.maximize(t0,c1_n,c1_d,maxim);
	porig.minimize(t1,c5_n,c5_d,minim);
	porig.maximize(t1,c2_n,c2_d,maxim);

	auto stop = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(stop - start);
	pi_comp_time += duration.count();
}

void A_Polyhedron::getEarlyEnd(Coefficient &coe_n, Coefficient &coe_d){
	coe_n = c5_n;
	coe_d = c5_d;
}

void A_Polyhedron::getEarlyBegin(Coefficient &coe_n, Coefficient &coe_d){
	coe_n = c6_n;
	coe_d = c6_d;
}

void A_Polyhedron::getLateEnd(Coefficient &coe_n, Coefficient &coe_d){
	coe_n = c2_n;
	coe_d = c2_d;
}

void A_Polyhedron::getLateBegin(Coefficient &coe_n, Coefficient &coe_d){
	coe_n = c1_n;
	coe_d = c1_d;
}