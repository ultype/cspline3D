#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include"ppform.h"
#include"RK4.h"
#define max_size 100
int main(int argc,char* argv[]){
	std::cout<<"main in\n";
	std::ifstream ifs(argv[1]);
	std::cout<<"read file ok\n";
	std::string trs;
	std::vector< std::vector<double> > points;
	int num;
	double value;
	while(ifs>>num){
		std::vector<double> temp;
		std::cout<<num<<" ";
		ifs>>value;
		temp.push_back(value);
		std::cout<<std::fixed<<value<<" ";
		ifs>>value;
		temp.push_back(value);
		std::cout<<std::fixed<<value<<" ";
		ifs>>value;
		temp.push_back(value);
		std::cout<<std::fixed<<value<<"\n";
		points.push_back(temp);
	}
	ifs.close();
		
	PPform pp(points);
	
	RK4 rk;
	std::cout<<"initial ok\n";
	rk.arrange(0.01,10,&pp);
	std::cout<<"arrange ok \n";
	rk.RKrunner();
	std::cout<<"RKrunner ok\n";
	rk.saveData();
	std::cout<<"saveData ok\n";
	//std::ofstream ofs("map.in", std::ios_base::app);
	//ofs<<"lastline\n";
	//ofs.close();
	
	int k1,k2,k3;
	for(int i=0;i<pp.size-1;i++){
		k1=3*i;
		k2=k1+1;
		k3=k2+1;
		std::cout<<i<<": "<<std::fixed<<pp.coefA[k1]<<" "<<pp.coefB[k1]<<" "<<pp.coefC[k1]<<" "<<pp.coefD[k1]<<"\n";
		std::cout<<i<<": "<<std::fixed<<pp.coefA[k2]<<" "<<pp.coefB[k2]<<" "<<pp.coefC[k2]<<" "<<pp.coefD[k2]<<"\n";
		std::cout<<i<<": "<<std::fixed<<pp.coefA[k3]<<" "<<pp.coefB[k3]<<" "<<pp.coefC[k3]<<" "<<pp.coefD[k3]<<"\n";
	}

	for(int i=0;i<rk.totalpoints;i++){
		std::cout<<i<<" :";
		for(int j=0;j<6;j++){
			std::cout<<rk.RKy[i][j]<<" ";
		}
		std::cout<<"\n";
	}
	return 0;
}
