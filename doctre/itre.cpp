#include "itre.hpp"

struct {
    bool operator()(I_Polyhedron &a, I_Polyhedron &b) const
    {   
    	if(a.begin < b.begin){
    		return true;
    	}else if(a.begin == b.begin){
    		return a.end < b.end;
    	}else{
    		return false;
    	}
    }  
}compareInterval;

struct {
    bool operator()(I_Polyhedron &a, I_Polyhedron &b) const
    {   
    	return a.begin < b.begin;
    }  
}compareBegin;

int findInterval(vector<I_Polyhedron> &iply, int start, int finish, int begin, int end){
	while(start<=finish){
		int mid = (start+finish)/2;
		if((iply[mid].begin == begin) && (iply[mid].end == end)){
			return mid;
		}else if(iply[mid].begin < begin){
			start = mid+1;
		}else if(iply[mid].begin == begin){
			if(iply[mid].end < end){
				start = mid+1;
			}else{
				finish = mid-1;
			}
		}else{
			finish = mid-1;
		}
	}
	return -1;
}

int findBegin(vector<I_Polyhedron> &iply, int start, int finish, int begin){
	while(start<=finish){
		int mid = (start+finish)/2;
		if(iply[mid].begin == begin){
			return mid;
		}else if(iply[mid].begin < begin){
			start = mid+1;
		}else{
			finish = mid-1;
		}
	}
	return -1;
}

I_Polyhedron ipoly_intersection(I_Polyhedron &iply1, I_Polyhedron &iply2){
	if((iply1.begin == iply2.begin) && (iply1.end == iply2.end)){
		C_Polyhedron pres = iply1.phed;
		pres.intersection_assign(iply2.phed);
		return I_Polyhedron(pres, iply1.begin, iply1.end);
	}else{
		cout<<"ipoly intersection anomaly"<<endl;
		exit(0);
		return iply1;
	}
}

I_Polyhedron ipoly_sequential_composition(I_Polyhedron &iply1, I_Polyhedron &iply2){
	if(iply1.end == iply2.begin){
		C_Polyhedron pres = iply1.phed;
		pres.intersection_assign(iply2.phed);
		return I_Polyhedron(pres, iply1.begin, iply2.end);
	}else{
		cout<<"ipoly sequential composition anomaly"<<endl;
		exit(0);
		return iply1;
	}
}

void ipoly_duration_restriction(I_Polyhedron &iply, Linear_Expression &l1, Linear_Expression &l2){
	iply.phed.add_constraint(l1 <= iply.end-iply.begin);
	iply.phed.add_constraint(iply.end-iply.begin <= l2);
}

vector<I_Polyhedron> itre_union(vector<I_Polyhedron> itre1, vector<I_Polyhedron> itre2){
	itre1.insert(itre1.end(), itre2.begin(), itre2.end());
	return itre1;
}

vector<I_Polyhedron> itre_intersection(vector<I_Polyhedron> &itre1, vector<I_Polyhedron> &itre2){
	std::sort(itre2.begin(), itre2.end(), compareInterval);
	vector<I_Polyhedron> ires;
	for(int i=0;i<itre1.size();i++){
		int begin1, end1;
		begin1 = itre1[i].begin;
		end1   = itre1[i].end;
		int i2 = findInterval(itre2, 0, itre2.size()-1, begin1, end1);
		if(i2!=-1){
			int j=i2;
			while( (j < itre2.size()) && (itre2[j].begin == begin1) && 
				(itre2[j].end == end1)){
				I_Polyhedron itemp = ipoly_intersection(itre1[i], itre2[j]);
				if(!itemp.phed.is_empty())
					ires.push_back(itemp);
				j++;
			}
			j = i2-1;
			while(j>=0 && (itre2[j].begin == begin1) && 
				(itre2[j].end == end1)){
				I_Polyhedron itemp = ipoly_intersection(itre1[i],itre2[j]);
				if(!itemp.phed.is_empty())
					ires.push_back(itemp);
				j--;
			}
		}
	}
	return ires;
}

vector<I_Polyhedron> itre_concatenation(vector<I_Polyhedron> &itre1, vector<I_Polyhedron> &itre2){
	std::sort(itre2.begin(), itre2.end(), compareBegin);
	vector<I_Polyhedron> ires;
	for(int i=0;i<itre1.size();i++){
		int end1;
		end1 = itre1[i].end;
		int i2 = findBegin(itre2, 0, itre2.size()-1, end1);
		if(i2!=-1){
			int j = i2;
			while(j < itre2.size() && end1 == itre2[j].begin){
				I_Polyhedron itemp = ipoly_sequential_composition(itre1[i], itre2[j]);
				if(!itemp.phed.is_empty())
					ires.push_back(itemp);
				j++;
			}
			j = i2-1;
			while(j>=0 && end1 == itre2[j].begin){
				I_Polyhedron itemp = ipoly_sequential_composition(itre1[i], itre2[j]);
				if(!itemp.phed.is_empty())
					ires.push_back(itemp);
				j--;
			}
		}
	}
	return ires;
}

vector<I_Polyhedron>& itre_duration_restriction(vector<I_Polyhedron> &itre, Linear_Expression &l1, 
		Linear_Expression &l2){
	for(int i=0;i<itre.size();i++){
		ipoly_duration_restriction(itre[i], l1, l2);
	}
	return itre;
}


int itre_contained(vector<I_Polyhedron> &plyNew, vector<I_Polyhedron> &plyOld){
	int all_contained = 1;
	for(int i=0;i<plyNew.size();i++){
		int contained = 0;	
		for(int j=0;j<plyOld.size();j++){
			contained = (contained || plyOld[j].contains(plyNew[i]));
		}
		all_contained = (all_contained && contained);
	}
	return all_contained;
}


vector<I_Polyhedron> itre_KleenePlusSquaring(vector<I_Polyhedron> &plyR){
	vector<I_Polyhedron> plyQ,plyP;
	plyQ = plyR;
	plyP = itre_concatenation(plyR, plyR);
	int count =0;
	while(!itre_contained(plyP, plyQ)){
		plyQ = itre_union(itre_union(plyP,plyQ), itre_concatenation(plyP,plyQ));
		plyP = itre_concatenation(plyP,plyP);
		cout<<"1,2:"<<plyP.size()<<","<<plyQ.size()<<endl;
		count++;
	}
	cout<<"Kleene count:"<<count<<endl;
	return plyQ;
}

void itre_print(vector<I_Polyhedron> &itre){
	cout<<"itre size:"<<itre.size()<<endl;
	for(int i=0;i<itre.size();i++){
		cout<<"t="<<itre[i].begin<<","<<"t'="<<itre[i].end<<endl;
		cout<<itre[i].phed<<endl;	
	}
}