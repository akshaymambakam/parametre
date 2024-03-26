#include "helperFunctions.hpp"

vector<C_Polyhedron> labelZoneList(vector<C_Polyhedron> vLabel, vector<C_Polyhedron> vZone, int numParams){
    vector<C_Polyhedron> vRes;
    for(int i=0;i<vLabel.size();i++){
        for(int j=0;j<vZone.size();j++){
            vRes.push_back(labelZoneCombination(vLabel[i],vZone[j],numParams));
        }
    }
    return vRes;
}


C_Polyhedron labelZoneCombination(C_Polyhedron plabel, C_Polyhedron matchZone, int numParams){
    Constraint_System csLabel = plabel.constraints();
    Constraint_System csZone  = matchZone.constraints();
    Variable t3(3+numParams),t4(3+numParams+1);

    Constraint_System csnew;

    Variable t0(0),t1(1);
    // Swap t0,t1 with t3 and t4 respectively
    for(Constraint_System_const_iterator csit=csLabel.begin();csit!=csLabel.end();csit++){
        Constraint ctemp=*csit;
        ctemp.swap_space_dimensions(t0,t3);
        ctemp.swap_space_dimensions(t1,t4);
        csnew.insert(ctemp);
    }
    // Add all constraints from matchZone
    for(Constraint_System_const_iterator csit=csZone.begin();csit!=csZone.end();csit++){
        Constraint ctemp=*csit;
        csnew.insert(ctemp);
    }
    // Add constraints (t0 <= t3 <= t1) and (t0 <= t4 <= t1)
    csnew.insert(t0<=t3);
    csnew.insert(t3<=t1);
    csnew.insert(t0<=t4);
    csnew.insert(t4<=t1);
    C_Polyhedron plyRes(csnew);
    plyRes.unconstrain(t3);
    plyRes.unconstrain(t4);
    return plyRes;
}

vector<C_Polyhedron> projectParams(vector<C_Polyhedron> ptre){
    for(int i=0;i<ptre.size();i++){
        Variable t0(0),t1(1);
        ptre[i].unconstrain(t0);
        ptre[i].unconstrain(t1);    
    }
    return ptre;
}

void ptreUnconstrain(vector<C_Polyhedron> &ptre, int numParams){
	for(int i=0;i<ptre.size();i++){
		for(int j=0;j<numParams;j++){
			Variable v(j+3);
			ptre[i].unconstrain(v);
		}
	}
}

int placesAfterDot(string sf){
	int count=0;
	int i=0;
	int dotFlag=0;
	while(sf[i]!=0){
		if(sf[i]=='.'){
			count = 0;
			dotFlag = 1;
		}
		count++;
		i++;
	}
	if(dotFlag){
		return count-1;
	}else{
		return 0;
	}
}

int exponentOf10(int expo){
	int exp10=1;
	for(int i=0;i<expo;i++){
		exp10*=10;
	}
	return exp10;
}

int removeDot(string fnum){
	char dup[100];
	int i=0,j=0;
	while(fnum[i]!=0){
		if(fnum[i]!='.'){
			dup[j]=fnum[i];
			j++;
		}
		i++;
	}
	dup[i]=0;
	return atoi(dup);
}

void readLabelFile(char *fName, vector<pair<int,int>> &labellingList){
    ifstream infile(fName);
    string line;
    int lineNum=0;
    while(std::getline(infile,line)){
        string s=line;
        int pos = s.find(",");
        int listart = atoi(s.substr(0,pos).c_str());
        int liend   = atoi(s.substr(s.find(",")+1).c_str());
        labellingList.push_back(pair<int,int>(listart,liend));
    }
    infile.close();
}

void readEventFile(char *fName, vector<int> &timeVec, vector<int> &eventNumVec){
    ifstream infile(fName);
    string line;
    int lineNum=0;
    while(std::getline(infile,line)){
        string s=line;
        int pos = s.find(",");
        int timeStamp = atoi(s.substr(0,pos).c_str());
        int eventNum = atoi(s.substr(s.find("y")+1).c_str());
        timeVec.push_back(timeStamp);
        eventNumVec.push_back(eventNum);
    }
    infile.close();
}

void readDataFile(char *fName, vector<int> &timeVec, vector<vector<string>> &valueVec){
	ifstream infile(fName);
    string delimiter = ",";
    string line;
    int lineNum=0;
    while (std::getline(infile, line))
    {
        size_t pos = 0;
        std::string token;
        int num=0;
        string s=line;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            if(num==0){
                timeVec.push_back(atoi(s.c_str()));
            }else{
                if(lineNum==0){
                    vector<string> tcVec;
                    tcVec.push_back(token);
                    valueVec.push_back(tcVec);
                }else{
                    valueVec[num-1].push_back(token);
                }
                
            }
            s.erase(0, pos + delimiter.length());
            num++;
        }
        if(lineNum==0){
            vector<string> tcVec;
            tcVec.push_back(s);
            valueVec.push_back(tcVec);
        }else{
            valueVec[num-1].push_back(s);
        }
        lineNum++;
    }
    infile.close();
}

vector<C_Polyhedron> labellingToZones(vector<Constraint> paramS, vector<pair<int,int>> &labellingList){
    Constraint_System zoneTemp;
    vector<C_Polyhedron> labZones;
    Variable t0(0);
    Variable t1(1);
    
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    for(int i=0;i<labellingList.size();i++){
        Constraint_System zonePush = zoneTemp;
        int ztstart = labellingList[i].first;
        int ztend   = labellingList[i].second;
        zonePush.insert(ztstart<=t0);
        zonePush.insert(t0<=ztend);
        zonePush.insert(ztstart<=t1);
        zonePush.insert(t1<=ztend);
        zonePush.insert(0<=(t1-t0));
        zonePush.insert((t1-t0)<=(ztend-ztstart));
        C_Polyhedron cpTemp(zonePush);
        labZones.push_back(cpTemp);
    }
    return labZones;
}

void readBoolZones(char *fName, vector<vector<C_Polyhedron>> &bZonelist, 
        vector<vector<C_Polyhedron>> &bNegZonelist){
    ifstream infile(fName);
    string line;
    int lineNum=0;
    while (std::getline(infile, line)){
        vector<int> timeVec;
        vector<vector<string>> valueVec;
        vector<Constraint> paramS;
        readDataFile(strdup(line.c_str()), timeVec, valueVec);
        vector<C_Polyhedron> ptreRet = boolToZone(paramS, timeVec, valueVec[0], 1);
        vector<C_Polyhedron> ptreRetNeg = boolToZone(paramS, timeVec, valueVec[0], 0);
        bZonelist.push_back(ptreRet);
        bNegZonelist.push_back(ptreRetNeg);
    }
    infile.close();
}

vector<C_Polyhedron> eventToZone(vector<Constraint> paramS, vector<int> timeVec,
    vector<int> eventVec, int eventNum){
    Variable t0(0);
    Variable t1(1);

    int prevts = 0;
    
    vector<C_Polyhedron> ptreRet;
    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    for(int i=0;i<timeVec.size();i++){
        if(eventVec[i]==eventNum){
            Constraint_System zonePush = zoneTemp;
            zonePush.insert(t0==prevts);
            zonePush.insert(t1==timeVec[i]);
            
            C_Polyhedron cpTemp(zonePush);
            ptreRet.push_back(cpTemp);    
        }
        prevts = timeVec[i];
    }

    return ptreRet;
}

vector<I_Polyhedron> eventToIntervalPoly(vector<Constraint> paramS, vector<int> timeVec,
    vector<int> eventVec, int eventNum){
    Variable t0(0);
    Variable t1(1);

    int prevts = 0;
    
    vector<I_Polyhedron> itreRet;
    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    for(int i=0;i<timeVec.size();i++){
        if(eventVec[i]==eventNum){
            I_Polyhedron cpTemp(C_Polyhedron(zoneTemp), prevts, timeVec[i]);
            itreRet.push_back(cpTemp);    
        }
        prevts = timeVec[i];
    }

    return itreRet;
}

vector<I_Polyhedron> epsToIntervalPoly(vector<Constraint> paramS, vector<int> timeVec,
    vector<int> eventVec){
    Variable t0(0);
    Variable t1(1);

    int prevts = 0;
    
    vector<I_Polyhedron> itreRet;
    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    for(int i=0;i<timeVec.size();i++){
        if(prevts == timeVec[i]){
            I_Polyhedron cpTemp(C_Polyhedron(zoneTemp), prevts, timeVec[i]);
            itreRet.push_back(cpTemp);
        }
        prevts = timeVec[i];
    }

    return itreRet;
}


vector<C_Polyhedron> boolToEdgeZone(vector<Constraint> paramS, vector<int> timeVec,
    vector<string> fvec, int prevOneIn){
    Variable t0(0);
    Variable t1(1);

    int prevOne = 0;
    
    vector<C_Polyhedron> ptreRet;
    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    for(int i=0;i<timeVec.size()-1;i++){
        int exp10 = exponentOf10(placesAfterDot(fvec[i]));
        int fvectemp = removeDot(fvec[i]);
        if(i==0){
            prevOne = fvectemp;
        }else{
            if(fvectemp != prevOne){
                if( (prevOneIn == 1 && fvectemp > prevOne) || 
                    (prevOneIn == 0 && fvectemp < prevOne) ){
                    Constraint_System zonePush = zoneTemp;
                    zonePush.insert(t0==timeVec[i]);
                    zonePush.insert(t1==timeVec[i]);
                    C_Polyhedron cpTemp(zonePush);
                    ptreRet.push_back(cpTemp);
                }
                prevOne = fvectemp;
            }
        }
    }
    return ptreRet;
}

vector<C_Polyhedron> boolToZone(vector<Constraint> paramS, vector<int> timeVec, 
    vector<string> fvec, int prevOneIn){
	
    Variable t0(0);
    Variable t1(1);
    
    vector<C_Polyhedron> ptreRet;
    int ztstart = 0, ztend = 0;
    int prevOne = 1-prevOneIn;
    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }
    
    for(int i=0;i<timeVec.size()-1;i++){
        int exp10 = exponentOf10(placesAfterDot(fvec[i]));
        int fvectemp = removeDot(fvec[i]);
        if(fvectemp==1-prevOneIn){
            if(prevOne == prevOneIn){
                Constraint_System zonePush = zoneTemp;
                zonePush.insert(ztstart<=t0);
                zonePush.insert(t0<=ztend);
                zonePush.insert(ztstart<=t1);
                zonePush.insert(t1<=ztend);
                zonePush.insert(0<=(t1-t0));
                zonePush.insert((t1-t0)<=(ztend-ztstart));
                C_Polyhedron cpTemp(zonePush);
                ptreRet.push_back(cpTemp);
            }
            prevOne = 1-prevOneIn;
        }else if(fvectemp==prevOneIn){
            if(prevOne == prevOneIn){
                ztend = timeVec[i+1];
            }else{
                prevOne = prevOneIn;
                ztstart = timeVec[i];
                ztend   = timeVec[i+1];
            }
        }else{
            cerr<<"Error:Non-Boolean predicate"<<endl;
            exit(0);
        }
    }
    if(prevOne == prevOneIn){
        Constraint_System zonePush = zoneTemp;
        zonePush.insert(ztstart<=t0);
        zonePush.insert(t0<=ztend);
        zonePush.insert(ztstart<=t1);
        zonePush.insert(t1<=ztend);
        zonePush.insert(0<=(t1-t0));
        zonePush.insert((t1-t0)<=(ztend-ztstart));
        C_Polyhedron cpTemp(zonePush);
        ptreRet.push_back(cpTemp);
    }
    return ptreRet;
}

vector<C_Polyhedron> aggrToZone(vector<Constraint> paramS, vector<int> timeVec, 
    vector<string> fvec){
    int prevOneIn = 1;
    
    Variable t0(0);
    Variable t1(1);

    vector<int> startVec,endVec;
    
    vector<C_Polyhedron> ptreRet;
    int ztstart = 0, ztend = 0;
    int prevOne = 1-prevOneIn;
    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }
    
    for(int i=0;i<timeVec.size()-1;i++){
        int exp10 = exponentOf10(placesAfterDot(fvec[i]));
        int fvectemp = removeDot(fvec[i]);
        if(fvectemp==1-prevOneIn){
            if(prevOne == prevOneIn){
                startVec.push_back(ztstart);
                endVec.push_back(ztend);
            }
            prevOne = 1-prevOneIn;
        }else if(fvectemp==prevOneIn){
            if(prevOne == prevOneIn){
                ztend = timeVec[i+1];
            }else{
                prevOne = prevOneIn;
                ztstart = timeVec[i];
                ztend   = timeVec[i+1];
            }
        }else{
            cerr<<"Error:Non-Boolean predicate"<<endl;
            exit(0);
        }
    }
    if(prevOne == prevOneIn){
        startVec.push_back(ztstart);
        endVec.push_back(ztend);
    }

    int oldEnd=0;
    for(int i=0;i<startVec.size();i++){
        if(i==0 && startVec[i]==0){
            oldEnd = endVec[i];
        }else{
            Constraint_System zonePush = zoneTemp;
            zonePush.insert(0<=t0);
            zonePush.insert(t0<=startVec[i]);
            zonePush.insert(oldEnd<=t1);
            zonePush.insert(t1<=endVec[i]);
            zonePush.insert(0<=(t1-t0));
            C_Polyhedron cpTemp(zonePush);
            ptreRet.push_back(cpTemp);

            oldEnd = endVec[i];
        }
    }

    if(endVec.size()==0 || endVec[endVec.size()-1] < timeVec[timeVec.size()-1]){
        Constraint_System zonePush = zoneTemp;
        zonePush.insert(0<=t0);
        zonePush.insert(oldEnd<=t1);
        zonePush.insert(t1<=timeVec[timeVec.size()-1]);
        zonePush.insert(0<=(t1-t0));
        C_Polyhedron cpTemp(zonePush);
        ptreRet.push_back(cpTemp);
    }

    return ptreRet;
}

/* x <> p */
vector<C_Polyhedron> porvZones(vector<int> timeVec, vector<string> fvec, vector<Constraint> paramS,
        int coeff, Linear_Expression le, int begin, int end, int mxFlag){
    vector<C_Polyhedron> res;
    Variable t0(0),t1(1);
    int mxIndex;
    float mxValue;
    for(int i=begin;i<=end;i++){
        if(i==begin){
            mxIndex = begin;
            mxValue = atof(fvec[i].c_str());
        }else{
            float currValue = atof(fvec[i].c_str());
            if(mxFlag == 3){
                if(currValue < mxValue){
                    mxValue = currValue;
                    mxIndex = i;
                }
            }else if(mxFlag == 1){
                if(currValue > mxValue){
                    mxValue = currValue;
                    mxIndex = i;
                }
            }
            
        }
    }

    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    int exp10 = exponentOf10(placesAfterDot(fvec[mxIndex]));
    int fvectemp = removeDot(fvec[mxIndex]);

    Constraint_System zTemp = zoneTemp;
    zTemp.insert(timeVec[begin]<=t0);
    zTemp.insert(t0<=timeVec[mxIndex+1]);
    zTemp.insert(timeVec[mxIndex]<=t1);
    zTemp.insert(t1<=timeVec[end+1]);
    zTemp.insert(t0<=t1);

    if(mxFlag == 3)
        zTemp.insert(le*exp10 <= fvectemp*coeff);
    else
        zTemp.insert(le*exp10 >= fvectemp*coeff);

    vector<C_Polyhedron> subRes1,subRes2;

    if(mxIndex > begin)
        subRes1 = porvZones(timeVec, fvec, paramS, coeff, le, begin, mxIndex-1, mxFlag);
    
    res.insert(res.end(), subRes1.begin(), subRes1.end());
    res.push_back(C_Polyhedron(zTemp));
    
    if(mxIndex < end)
        subRes2 = porvZones(timeVec, fvec, paramS, coeff, le, mxIndex+1, end, mxFlag);
    
    res.insert(res.end(), subRes2.begin(), subRes2.end());

    return res;
}


/* (Max x >= c) or (Min x <= c) (Max = 1) (Min = 2) */
vector<C_Polyhedron> aggrZones(vector<int> timeVec, vector<string> fvec, vector<Constraint> paramS,
        int coeff, Linear_Expression le, int begin, int end, int mxFlag){
    vector<C_Polyhedron> res;
    Variable t0(0),t1(1);
    int mxIndex;
    float mxValue;
    for(int i=begin;i<=end;i++){
        if(i==begin){
            mxIndex = begin;
            mxValue = atof(fvec[i].c_str());
        }else{
            float currValue = atof(fvec[i].c_str());
            if(mxFlag == 2){
                if(currValue < mxValue){
                    mxValue = currValue;
                    mxIndex = i;
                }
            }else if(mxFlag == 1){
                if(currValue > mxValue){
                    mxValue = currValue;
                    mxIndex = i;
                }
            }
            
        }
    }

    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    int exp10 = exponentOf10(placesAfterDot(fvec[mxIndex]));
    int fvectemp = removeDot(fvec[mxIndex]);

    Constraint_System zTemp = zoneTemp;
    zTemp.insert(timeVec[begin]<=t0);
    zTemp.insert(t0<=timeVec[mxIndex+1]);
    zTemp.insert(timeVec[mxIndex]<=t1);
    zTemp.insert(t1<=timeVec[end+1]);
    zTemp.insert(t0<=t1);
    if(mxFlag == 1)
        zTemp.insert(le*exp10 <= fvectemp*coeff);
    else
        zTemp.insert(le*exp10 >= fvectemp*coeff);
    
    C_Polyhedron rhed(zTemp); 
    if(!rhed.is_empty())
        res.push_back(rhed);

    vector<C_Polyhedron> subRes1,subRes2;
    if(mxIndex > begin)
        subRes1 = aggrZones(timeVec, fvec, paramS, coeff, le, begin, mxIndex-1, mxFlag);
    if(mxIndex < end)
        subRes2 = aggrZones(timeVec, fvec, paramS, coeff, le, mxIndex+1, end, mxFlag);
    
    res.insert(res.end(), subRes1.begin(), subRes1.end());
    res.insert(res.end(), subRes2.begin(), subRes2.end());

    return res;
}

vector<C_Polyhedron> diffLeqZones(vector<int> timeVec, vector<string> fvec, vector<Constraint> paramS,
        int coeff, Linear_Expression le){
    vector<C_Polyhedron> res;
    Variable t0(0),t1(1);

    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    int minN, minD, maxN, maxD;
    for(int i=0;i<timeVec.size()-1;i++){
        int minN = removeDot(fvec[i]);
        int minD = exponentOf10(placesAfterDot(fvec[i]));

        int maxN = removeDot(fvec[i]);
        int maxD = exponentOf10(placesAfterDot(fvec[i]));

        int j;
        for(j=i;j<timeVec.size()-1;j++){
            int currN = removeDot(fvec[j]);
            int currD = exponentOf10(placesAfterDot(fvec[j]));

            if(currN*maxD > maxN*currD){ // Update max
                maxN = currN;
                maxD = currD;
            }
            if(currN*minD < minN*currD){ // Update min
                minN = currN;
                minD = currD;
            }

            Constraint_System zt;
            zt.insert(maxN*minD*coeff-minN*maxD*coeff <= le*maxD*minD);
            C_Polyhedron cpt(zt);
            if(cpt.is_empty()){
                break;
            }
        }
        Constraint_System zoneF = zoneTemp;
        zoneF.insert(timeVec[i]<=t0);
        zoneF.insert(t0<=timeVec[i+1]);
        zoneF.insert(t0<=t1);
        zoneF.insert(t1<=timeVec[j]);
        C_Polyhedron cpTemp(zoneF);
        res.push_back(cpTemp);
    }

    return res;
}


vector<C_Polyhedron> diffGeqZones(vector<int> timeVec, vector<string> fvec, vector<Constraint> paramS,
        int coeff, Linear_Expression le){
    vector<C_Polyhedron> res;
    Variable t0(0),t1(1);

    Constraint_System zoneTemp;
    for(int i=0;i<paramS.size();i++){
        zoneTemp.insert(paramS[i]);
    }

    int minN, minD, maxN, maxD;
    for(int i=0;i<timeVec.size()-1;i++){
        int minN = removeDot(fvec[i]);
        int minD = exponentOf10(placesAfterDot(fvec[i]));

        int maxN = removeDot(fvec[i]);
        int maxD = exponentOf10(placesAfterDot(fvec[i]));

        int j;
        for(j=i;j<timeVec.size()-1;j++){
            int currN = removeDot(fvec[j]);
            int currD = exponentOf10(placesAfterDot(fvec[j]));

            if(currN*maxD > maxN*currD){ // Update max
                maxN = currN;
                maxD = currD;
            }
            if(currN*minD < minN*currD){ // Update min
                minN = currN;
                minD = currD;
            }

            Constraint_System zt;
            zt.insert(maxN*minD*coeff-minN*maxD*coeff < le*maxD*minD);
            C_Polyhedron cpt(zt);
            if(cpt.is_empty()){
                break;
            }
        }

        if(j!=timeVec.size()-1){
            Constraint_System zoneF = zoneTemp;
            zoneF.insert(timeVec[i]<=t0);
            zoneF.insert(t0<=timeVec[i+1]);
            zoneF.insert(t0<=t1);
            zoneF.insert(timeVec[j]<=t1);
            zoneF.insert(t1<=timeVec[timeVec.size()-1]);
            C_Polyhedron cpTemp(zoneF);
            res.push_back(cpTemp);    
        }
        
    }

    return res;
}


vector<C_Polyhedron> poly_list_filter(vector<C_Polyhedron> plyIn){
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


vector<C_Polyhedron> poly_list_intersection(vector<C_Polyhedron> &ptre1, vector<C_Polyhedron> &ptre2){
    vector<C_Polyhedron> ptre_inter;
    for(int i=0;i<ptre1.size();i++){
        for(int j=0;j<ptre2.size();j++){
            C_Polyhedron ply1 = ptre1[i];
            C_Polyhedron ply2 = ptre2[j];

            ply1.intersection_assign(ply2);
            
            bool inter_empty = ply1.is_empty();
            if(!inter_empty){
                ptre_inter.push_back(ply1);
            }
        }
    }
    return poly_list_filter(ptre_inter);
}

vector<C_Polyhedron> ptre_label_project(vector<C_Polyhedron> ptre, vector<Constraint> paramS,
    vector<pair<int,int>> labelIntervals){
    vector<C_Polyhedron> plres;

    for(int i=0;i<labelIntervals.size();i++){
        Variable t0(0),t1(1);
        Constraint_System zoneTemp;
        for(int i=0;i<paramS.size();i++){
            zoneTemp.insert(paramS[i]);
        }
        zoneTemp.insert(t0==labelIntervals[i].first);
        zoneTemp.insert(t1==labelIntervals[i].second);
        C_Polyhedron ci(zoneTemp);
        vector<C_Polyhedron> tpl;
        tpl.push_back(ci);
        tpl = poly_list_intersection(tpl, ptre);

        vector<C_Polyhedron> tprojl = projectParams(tpl);
        if(i==0){
            plres = tprojl;
        }else{
            plres = poly_list_intersection(plres, tprojl);
        }
    }
    
    return plres;
}

void compute_subsignal(pair<int,int> interval, vector<int> orig_tvec, vector<string> orig_fvec, vector<int> &new_tvec,
    vector<string> &new_fvec){
    for(int i=0;i<orig_tvec.size();i++){
        int l = interval.first;
        int u = interval.second;
        if((i-2>0) && (i+1)<orig_tvec.size() && (orig_tvec[i+1]>=l) && (orig_tvec[i-2]<=u)){
            new_tvec.push_back(orig_tvec[i]);
            new_fvec.push_back(orig_fvec[i]);
        }
    }
}