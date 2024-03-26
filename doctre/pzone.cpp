#include "pzone.hpp"

pzone ptre_poly2pzone(Constraint_System cs){
	vector<Linear_Expression> c1,c2,c3,c4,c5,c6;

	Variable t0(0),t1(1),t2(2);
	Constraint_System cspoly;

	for(Constraint_System_const_iterator csit=cs.begin();csit!=cs.end();csit++){
		Linear_Expression le1(csit->expression());
		Linear_Expression lePlus = t0-t1;
		if(csit->coefficient(t0)==1 && csit->coefficient(t1)==0){
			le1 = le1 + (-t0);
			c6.push_back(-le1);
		}else if(csit->coefficient(t0)==-1 && csit->coefficient(t1)==0){
			le1 = le1 + t0;
			c1.push_back(le1);
		}else if(csit->coefficient(t0)==0 && csit->coefficient(t1)==1){
			le1 = le1 + (-t1);
			c5.push_back(-le1);
		}else if(csit->coefficient(t0)==0 && csit->coefficient(t1)==-1){
			le1 = le1 + t1;
			c2.push_back(le1);
		}else if(csit->coefficient(t0)==-1 && csit->coefficient(t1)==1){
			le1 = le1 + (t0-t1);
			c4.push_back(-le1);
		}else if(csit->coefficient(t0)==1 && csit->coefficient(t1)==-1){
			le1 = le1 + (-t0+t1);
			c3.push_back(le1);
		}else if(csit->coefficient(t0)==0 && csit->coefficient(t1)==0){
			cspoly.insert(*csit);
		}else{
			cout<<"Error with constraint:"<<(*csit)<<endl;
			exit(0);
		}
	}
	return pzone(c1,c2,c3,c4,c5,c6,C_Polyhedron(cspoly));
}

bool minimize_le(C_Polyhedron cp, Linear_Expression &le, Coefficient &coe_n, Coefficient &coe_d){
	bool minim;
	return cp.minimize(le, coe_n, coe_d, minim);
}

bool maximize_le(C_Polyhedron cp, Linear_Expression &le, Coefficient &coe_n, Coefficient &coe_d){
	bool maxim;
	return cp.maximize(le, coe_n, coe_d, maxim);
}

void coeff_minmax(int min, Coefficient coe1_n, Coefficient coe1_d, Coefficient coe2_n, Coefficient coe2_d, 
	Coefficient &cm_n, Coefficient &cm_d){
	if(coe1_n*coe2_d < coe2_n*coe1_d){
		if(min){
			cm_n = coe1_n;
			cm_d = coe1_d;
		}else{
			cm_n = coe2_n;
			cm_d = coe2_d;
		}
	}else{
		if(min){
			cm_n = coe2_n;
			cm_d = coe2_d;	
		}else{
			cm_n = coe1_n;
			cm_d = coe1_d;
		}
	}
}
void linexl_minmax(int min, C_Polyhedron cp, Coefficient &cm_n, Coefficient &cm_d, vector<Linear_Expression> &lel){
	for(int i=0;i<lel.size();i++){
		if(i==0){
			if(min)
				minimize_le(cp, lel[i], cm_n, cm_d);
			else
				maximize_le(cp, lel[i], cm_n, cm_d);
		}else{
			Coefficient ct_n, ct_d;
			if(min)
				minimize_le(cp, lel[i], ct_n, ct_d);
			else
				maximize_le(cp, lel[i], ct_n, ct_d);
			coeff_minmax(min, cm_n, cm_d, ct_n, ct_d, cm_n, cm_d);
		}
	}
}

pzone::pzone(vector<Linear_Expression> a,vector<Linear_Expression> b,vector<Linear_Expression> c,
		vector<Linear_Expression> d,vector<Linear_Expression> e,vector<Linear_Expression> f,
		C_Polyhedron pcons){
	this->cbounds.push_back(a);
	this->cbounds.push_back(b);
	this->cbounds.push_back(c);
	this->cbounds.push_back(d);
	this->cbounds.push_back(e);
	this->cbounds.push_back(f);
	this->cons = pcons;

	linexl_minmax(0, pcons, c1_n, c1_d, a);
	linexl_minmax(0, pcons, c2_n, c2_d, b);
	linexl_minmax(0, pcons, c3_n, c3_d, c);
	linexl_minmax(1, pcons, c4_n, c4_d, d);
	linexl_minmax(1, pcons, c5_n, c5_d, e);
	linexl_minmax(1, pcons, c6_n, c6_d, f);
}

pzone::pzone(vector<vector<Linear_Expression>> cbounds,C_Polyhedron pcons){
	this->cbounds=cbounds;
	this->cons = pcons;

	linexl_minmax(0, pcons, c1_n, c1_d, cbounds[0]);
	linexl_minmax(0, pcons, c2_n, c2_d, cbounds[1]);
	linexl_minmax(0, pcons, c3_n, c3_d, cbounds[2]);
	linexl_minmax(1, pcons, c4_n, c4_d, cbounds[3]);
	linexl_minmax(1, pcons, c5_n, c5_d, cbounds[4]);
	linexl_minmax(1, pcons, c6_n, c6_d, cbounds[5]);
}

void pzone::getEarlyEnd(Coefficient &coe_n, Coefficient &coe_d){
	coe_n = c5_n;
	coe_d = c5_d;
}

void pzone::getEarlyBegin(Coefficient &coe_n, Coefficient &coe_d){
	coe_n = c6_n;
	coe_d = c6_d;
}
 
void pzone::getLateEnd(Coefficient &coe_n, Coefficient &coe_d){
	coe_n = c2_n;
	coe_d = c2_d;
}

void pzone::getLateBegin(Coefficient &coe_n, Coefficient &coe_d){
	coe_n = c1_n;
	coe_d = c1_d;
}   

void pzone::print()
{
    // write obj to stream
    for(int i=0;i<NUM_BOUNDS;i++){
    	for(int j=0;j<this->cbounds[i].size();j++){
    		cout<<this->cbounds[i][j]<<";";
    	}
    	cout<<endl;
    }
    cout<<this->cons.minimized_constraints()<<endl;
    cout<<"^^^^^^^^"<<endl;
    cout<<"c1:"<<Linear_Expression(c1_n)<<"/"<<Linear_Expression(c1_d)<<endl;
    cout<<"c2:"<<Linear_Expression(c2_n)<<"/"<<Linear_Expression(c2_d)<<endl;
    cout<<"c3:"<<Linear_Expression(c3_n)<<"/"<<Linear_Expression(c3_d)<<endl;
    cout<<"c4:"<<Linear_Expression(c4_n)<<"/"<<Linear_Expression(c4_d)<<endl;
    cout<<"c5:"<<Linear_Expression(c5_n)<<"/"<<Linear_Expression(c5_d)<<endl;
    cout<<"c6:"<<Linear_Expression(c6_n)<<"/"<<Linear_Expression(c6_d)<<endl;
    cout<<"<><><><>\n"<<endl;
}
void pzone::refineTauto(int iOrig, int i1, char sn, int i2, char ineq){
	// Get back to zero indexing.
	iOrig--;
	i1--;
	i2--;

	vector<Linear_Expression> cOrig = this->cbounds[iOrig];
	vector<Linear_Expression> cFirst = this->cbounds[i1];
	vector<Linear_Expression> cSecond = this->cbounds[i2];

	vector<Linear_Expression> rc;

	for(int i=0;i<cOrig.size();i++){
		int tautoAdd = 1;
		for(int j=0;j<cOrig.size();j++){
			if(ineq == '>'){
				if((cOrig[j]<cOrig[i]).is_tautological()){
					tautoAdd = 0;
				}
			}else if(ineq == '<'){
				if((cOrig[j]>cOrig[i]).is_tautological()){
					tautoAdd = 0;
				}
			}
			if((cOrig[j]==cOrig[i]).is_tautological() && j<i){
				tautoAdd = 0;
			}
		}
		for(int k=0;k<cFirst.size();k++){
			for(int l=0;l<cSecond.size();l++){
				Linear_Expression texpr;
				if(sn == '+')
					texpr = cFirst[k]+cSecond[l];
				else
					texpr = cFirst[k]-cSecond[l];

				if(ineq == '>'){
					if((texpr<cOrig[i]).is_tautological()){
						tautoAdd = 0;
					}
				}else if(ineq == '<'){
					if((texpr>cOrig[i]).is_tautological()){
						tautoAdd = 0;
					}
				}
				if((texpr==cOrig[i]).is_tautological()){
					tautoAdd = 0;
				}
			}
		}
		if(tautoAdd){
			rc.push_back(cOrig[i]);
		}
	}
	this->cbounds[iOrig] = rc;
}

void pzone::refineB(int iOrig, int i1, char sn, int i2, char ineq){
	// Get back to zero indexing.
	iOrig--;
	i1--;
	i2--;

	vector<Linear_Expression> cOrig = this->cbounds[iOrig];
	vector<Linear_Expression> cFirst = this->cbounds[i1];
	vector<Linear_Expression> cSecond = this->cbounds[i2];

	vector<Linear_Expression> rc;
	for(int i=0;i<cOrig.size();i++){
		int toAdd = 1;
		for(int j=0;j<cOrig.size();j++){
			C_Polyhedron tcp = this->cons;
			if(j!=i){
				if(ineq == '>')
					tcp.add_constraint(cOrig[j]>=cOrig[i]);
				else
					tcp.add_constraint(cOrig[j]<=cOrig[i]);
			}

			if(tcp.is_empty()){
				toAdd = 0;
			}
		}
		for(int k=0;k<cFirst.size();k++){
			for(int l=0;l<cSecond.size();l++){
				C_Polyhedron tcp = this->cons;
				Linear_Expression texpr;
				if(sn == '+')
					texpr = cFirst[k]+cSecond[l];
				else
					texpr = cFirst[k]-cSecond[l];
				
				if(ineq == '>')
					tcp.add_constraint(texpr>=cOrig[i]);
				else
					tcp.add_constraint(texpr<=cOrig[i]);

				if(tcp.is_empty()){
					toAdd = 0;
				}
			}
		}
		if(toAdd){
			rc.push_back(cOrig[i]);
		}
	}
	this->cbounds[iOrig] = rc;
}

void pzone::refineC(int iOrig, int i1, char sn, int i2, char ineq){
	// Get back to zero indexing.
	iOrig--;
	i1--;
	i2--;

	vector<Linear_Expression> cOrig = this->cbounds[iOrig];
	vector<Linear_Expression> cFirst = this->cbounds[i1];
	vector<Linear_Expression> cSecond = this->cbounds[i2];

	vector<Linear_Expression> rc;

	for(int i=0;i<cOrig.size();i++){
		Constraint_System cs = (this->cons).minimized_constraints();
		int tautoAdd = 1;
		for(int j=0;j<cOrig.size();j++){
			if(j!=i){
				if(ineq == '>')
					cs.insert(cOrig[j]>=cOrig[i]);
				else
					cs.insert(cOrig[j]<=cOrig[i]);
			}
			if(j<i && (cOrig[j]==cOrig[i]).is_tautological()){
				tautoAdd = 0;
			}
		}
		for(int k=0;k<cFirst.size();k++){
			for(int l=0;l<cSecond.size();l++){
				Linear_Expression texpr;
				if(sn == '+')
					texpr = cFirst[k]+cSecond[l];
				else
					texpr = cFirst[k]-cSecond[l];
				
				if(ineq == '>')
					cs.insert(texpr>=cOrig[i]);
				else
					cs.insert(texpr<=cOrig[i]);
			}
		}
		if(!C_Polyhedron(cs).is_empty() && (tautoAdd)){
			rc.push_back(cOrig[i]);
		}
	}
	this->cbounds[iOrig] = rc;
}

void pzone::refine(){
	/*
	refineTauto(1,2,'-',4,'>');
	refineTauto(2,3,'+',1,'>');
	refineTauto(3,2,'-',6,'>');
	refineTauto(4,5,'-',1,'<');
	refineTauto(5,4,'+',6,'<');
	refineTauto(6,5,'-',3,'<');
	refineC(1,2,'-',4,'>');
	refineC(2,3,'+',1,'>');
	refineC(3,2,'-',6,'>');
	refineC(4,5,'-',1,'<');
	refineC(5,4,'+',6,'<');
	refineC(6,5,'-',3,'<');
	*/
}

int pzone::is_empty(){
	C_Polyhedron cptemp = this->to_polyhedron();
	if(!cptemp.is_empty()){
		*this = ptre_poly2pzone(cptemp.minimized_constraints());	
	}
	return cptemp.is_empty();
}

int pzone::is_emptyBackup(){

	Constraint_System cs(this->cons.minimized_constraints());
	for(int i=0;i<this->cbounds[1].size();i++){
		for(int j=0;j<this->cbounds[4].size();j++){
			cs.insert(this->cbounds[4][j]<=this->cbounds[1][i]);	
		}
	}
	for(int i=0;i<this->cbounds[0].size();i++){
		for(int j=0;j<this->cbounds[5].size();j++){
			cs.insert(this->cbounds[5][j]<=this->cbounds[0][i]);	
		}
	}
	for(int i=0;i<this->cbounds[2].size();i++){
		for(int j=0;j<this->cbounds[3].size();j++){
			cs.insert(this->cbounds[3][j]<=this->cbounds[2][i]);	
		}
	}

	for(int i=0;i<this->cbounds[4].size();i++){
		for(int j=0;j<this->cbounds[0].size();j++){
			for(int k=0;k<this->cbounds[2].size();k++){
				cs.insert(this->cbounds[4][i]<=this->cbounds[0][j]+this->cbounds[2][k]);
			}
		}
	}

	for(int i=0;i<this->cbounds[1].size();i++){
		for(int j=0;j<this->cbounds[3].size();j++){
			for(int k=0;k<this->cbounds[5].size();k++){
				cs.insert(this->cbounds[3][j]+this->cbounds[5][k]<=this->cbounds[1][i]);
			}
		}
	}

	C_Polyhedron temp_poly(cs);
	if(!temp_poly.is_empty()){
		C_Polyhedron cptemp = this->to_polyhedron();
		*this = ptre_poly2pzone(cptemp.minimized_constraints());
	}
	
	// *this = ptre_poly2pzone(temp_poly.minimized_constraints());
	return temp_poly.is_empty();
}

/*

int pzone::is_empty(){
	Variable t0(0);
	Variable t1(1);
	Variable t2(2);
	Variable x(3);
	Variable y(3);

	C_Polyhedron temp_poly=this->cons;
	Constraint_System cs(this->cons.minimized_constraints());
	for(int i=0;i<this->cbounds[0].size();i++){
		temp_poly.add_constraint(t0<=this->cbounds[0][i]);
		cs.insert(t0<=this->cbounds[0][i]);
	}
	for(int i=0;i<this->cbounds[1].size();i++){
		temp_poly.add_constraint(t1<=this->cbounds[1][i]);
		cs.insert(t1<=this->cbounds[1][i]);
	}
	for(int i=0;i<this->cbounds[2].size();i++){
		temp_poly.add_constraint(t1-t0<=this->cbounds[2][i]);
		cs.insert(t1-t0<=this->cbounds[2][i]);
	}
	for(int i=0;i<this->cbounds[3].size();i++){
		temp_poly.add_constraint(this->cbounds[3][i]<=t1-t0);
		cs.insert(this->cbounds[3][i]<=t1-t0);
	}
	for(int i=0;i<this->cbounds[4].size();i++){
		temp_poly.add_constraint(this->cbounds[4][i]<=t1);
		cs.insert(this->cbounds[4][i]<=t1);
	}
	for(int i=0;i<this->cbounds[5].size();i++){
		temp_poly.add_constraint(this->cbounds[5][i]<=t0);
		cs.insert(this->cbounds[5][i]<=t0);
	}
	return temp_poly.is_empty();
}
*/

C_Polyhedron pzone::to_polyhedron(){
	Variable t0(0);
	Variable t1(1);
	Variable t2(2);
	Variable x(3);
	Variable y(3);

	Constraint_System cs=this->cons.minimized_constraints();
	for(int i=0;i<this->cbounds[0].size();i++){
		cs.insert(t0<=this->cbounds[0][i]);
	}
	for(int i=0;i<this->cbounds[1].size();i++){
		cs.insert(t1<=this->cbounds[1][i]);
	}
	for(int i=0;i<this->cbounds[2].size();i++){
		cs.insert(t1-t0<=this->cbounds[2][i]);
	}
	for(int i=0;i<this->cbounds[3].size();i++){
		cs.insert(this->cbounds[3][i]<=t1-t0);
	}
	for(int i=0;i<this->cbounds[4].size();i++){
		cs.insert(this->cbounds[4][i]<=t1);
	}
	for(int i=0;i<this->cbounds[5].size();i++){
		cs.insert(this->cbounds[5][i]<=t0);
	}
	return C_Polyhedron(cs);
}


// Assuming <Z>_[A,B]
pzone pzone::duration_restriction(Linear_Expression A, Linear_Expression B){
	pzone pdurest(this->cbounds,this->cons);
	vector<Linear_Expression> newc3 = pdurest.cbounds[2];
	newc3.push_back(B);
	pdurest.cbounds[2] = newc3;
	
	vector<Linear_Expression> newc4 = pdurest.cbounds[3];
	newc4.push_back(A);
	pdurest.cbounds[3] = newc4;

	return pdurest;
}

pzone pzone::intersection(pzone &p) {
	C_Polyhedron cons_inter = this->cons;
	cons_inter.add_constraints(p.cons.minimized_constraints());
	vector<vector<Linear_Expression>> new_cbounds;
	
	if(cons_inter.is_empty()){
		return pzone(new_cbounds, cons_inter);
	}

	for(int nc=0;nc<NUM_BOUNDS;nc++){
		vector<Linear_Expression> cboundsTemp = this->cbounds[nc];
		cboundsTemp.insert(cboundsTemp.end(), p.cbounds[nc].begin(), p.cbounds[nc].end());
		new_cbounds.push_back(cboundsTemp);
	}
	return pzone(new_cbounds,cons_inter);
}

pzone pzone::sequential_composition(pzone &p){
	C_Polyhedron cons_poly1 = this->cons;
	C_Polyhedron cons_poly2 = p.cons;

	C_Polyhedron cons_seqcomp = cons_poly1;
	cons_seqcomp.add_constraints(cons_poly2.minimized_constraints());

	vector<vector<Linear_Expression>> cbtemps;
	if(cons_seqcomp.is_empty()){
		return pzone(cbtemps, cons_seqcomp);
	}

	// Constraint_System cstemp;

	// Start with adding new constraints. c5<=c'1, c5<=c2, c'6<=c2, c'6<=c'1,  c4<=c3, c'4<=c'3,
	for(int i=0;i<this->cbounds[4].size();i++){
		for(int j=0;j<p.cbounds[0].size();j++){
			cons_seqcomp.add_constraint(Constraint(this->cbounds[4][i]<=p.cbounds[0][j]));
			// cstemp.insert(Constraint(this->cbounds[4][i]<=p.cbounds[0][j]));
		}
		for(int j=0;j<this->cbounds[1].size();j++){
			cons_seqcomp.add_constraint(Constraint(this->cbounds[4][i]<=this->cbounds[1][j]));
			// cstemp.insert(Constraint(this->cbounds[4][i]<=this->cbounds[1][j]));
		}
	}

	for(int i=0;i<p.cbounds[5].size();i++){
		for(int j=0;j<this->cbounds[1].size();j++){
			cons_seqcomp.add_constraint(Constraint(p.cbounds[5][i]<=this->cbounds[1][j]));
			// cstemp.insert(Constraint(p.cbounds[5][i]<=this->cbounds[1][j]));
		}
		for(int j=0;j<p.cbounds[0].size();j++){
			cons_seqcomp.add_constraint(Constraint(p.cbounds[5][i]<=p.cbounds[0][j]));
			// cstemp.insert(Constraint(p.cbounds[5][i]<=p.cbounds[0][j]));
		}
	}

	for(int i=0;i<this->cbounds[3].size();i++){
		for(int j=0;j<this->cbounds[2].size();j++){
			cons_seqcomp.add_constraint(Constraint(this->cbounds[3][i]<=this->cbounds[2][j]));
			// cstemp.insert(Constraint(this->cbounds[3][i]<=this->cbounds[2][j]));
		}
	}
	for(int i=0;i<p.cbounds[3].size();i++){
		for(int j=0;j<p.cbounds[2].size();j++){
			cons_seqcomp.add_constraint(Constraint(p.cbounds[3][i]<=p.cbounds[2][j]));
			// cstemp.insert(Constraint(p.cbounds[3][i]<=p.cbounds[2][j]));
		}
	}

	// Now create the new cbounds lists 

	// c1new = min(c1,c'1-c4,c2-c4)
	vector<Linear_Expression> newc1 = this->cbounds[0];
	for(int i=0;i<this->cbounds[3].size();i++){
		for(int j=0;j<p.cbounds[0].size();j++){
			Linear_Expression c1PrimeMinusC4 = p.cbounds[0][j]-this->cbounds[3][i];
			newc1.push_back(c1PrimeMinusC4);
		}
		for(int j=0;j<this->cbounds[1].size();j++){
			Linear_Expression c2MinusC4 = this->cbounds[1][j]-this->cbounds[3][i];
			newc1.push_back(c2MinusC4);
		}
	}

	// c2new = min(c'2,c2+c'3,c'1+c'3)
	vector<Linear_Expression> newc2 = p.cbounds[1];
	for(int i=0;i<p.cbounds[2].size();i++){
		for(int j=0;j<this->cbounds[1].size();j++){
			Linear_Expression c2PlusC3Prime = this->cbounds[1][j]+p.cbounds[2][i];
			newc2.push_back(c2PlusC3Prime);
		}
		for(int j=0;j<p.cbounds[0].size();j++){
			Linear_Expression c1PrimePlusC3Prime = p.cbounds[0][j]+p.cbounds[2][i];
			newc2.push_back(c1PrimePlusC3Prime);
		}
	}

	// c3new = c3+c'3
	vector<Linear_Expression> newc3;
	for(int i=0;i<this->cbounds[2].size();i++){
		for(int j=0;j<p.cbounds[2].size();j++){
			Linear_Expression c3Plusc3Prime = this->cbounds[2][i]+p.cbounds[2][j];
			newc3.push_back(c3Plusc3Prime);
		}
	}

	// c4new = c4+c'4
	vector<Linear_Expression> newc4;
	for(int i=0;i<this->cbounds[3].size();i++){
		for(int j=0;j<p.cbounds[3].size();j++){
			Linear_Expression c4Plusc4Prime = this->cbounds[3][i]+p.cbounds[3][j];
			newc4.push_back(c4Plusc4Prime);
		}
	}

	// c5new = max(c'5,c5+c'4,c'6+c'4)
	vector<Linear_Expression> newc5 = p.cbounds[4];
	for(int i=0;i<p.cbounds[3].size();i++){
		for(int j=0;j<this->cbounds[4].size();j++){
			Linear_Expression c5PlusC4Prime = this->cbounds[4][j]+p.cbounds[3][i];
			newc5.push_back(c5PlusC4Prime);
		}
		for(int j=0;j<p.cbounds[5].size();j++){
			Linear_Expression c6PrimePlusC4Prime = p.cbounds[5][j]+p.cbounds[3][i];
			newc5.push_back(c6PrimePlusC4Prime);
		}
	}

	// c6new = max(c6,c'6-c3,c5-c3)
	vector<Linear_Expression> newc6 = this->cbounds[5];
	for(int i=0;i<this->cbounds[2].size();i++){
		for(int j=0;j<p.cbounds[5].size();j++){
			Linear_Expression c6PrimeMinusC3 = p.cbounds[5][j]-this->cbounds[2][i];
			newc6.push_back(c6PrimeMinusC3);
		}
		for(int j=0;j<this->cbounds[4].size();j++){
			Linear_Expression c5MinusC3 = this->cbounds[4][j]-this->cbounds[2][i];
			newc6.push_back(c5MinusC3);
		}
	}

	return pzone(newc1,newc2,newc3,newc4,newc5,newc6,cons_seqcomp);
}
