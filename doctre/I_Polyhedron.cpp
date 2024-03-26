#include "I_Polyhedron.hpp"

I_Polyhedron::I_Polyhedron(C_Polyhedron phed, int begin, int end){
	this->phed = phed;
	this->begin   = begin;
	this->end   = end;
}

int I_Polyhedron::contains(I_Polyhedron &opoly){
	return (this->begin == opoly.begin) && (this->end == opoly.end) && (this->phed).contains(opoly.phed);
}