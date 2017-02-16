#include<iostream>
#include"cspline.c"

int main(int argc,char *argv[]){
	struct PPform pp0;
	int n=9;
	double** points=(double**)malloc(n*sizeof(double*));
	for(int i=0;i<n;i++){
		points[i]=(double*)malloc(3*sizeof(double));
	}
	
	points[0][0]=0; points[0][1]=0; points[0][2]=1;
	points[1][0]=1; points[1][1]=0; points[1][2]=0;
	points[2][0]=1; points[2][1]=1; points[2][2]=1;
	points[3][0]=0; points[3][1]=2; points[3][2]=0;
	points[4][0]=-1; points[4][1]=1; points[4][2]=1;
	points[5][0]=-1; points[5][1]=0; points[5][2]=0;
	points[6][0]=0; points[6][1]=-1; points[6][2]=1;
	points[7][0]=0; points[7][1]=-2; points[7][2]=0;
	points[8][0]=0; points[8][1]=0; points[8][2]=1;
	pp0.points=points;
	pp0.n=n;
	cspline3D_periodic(&pp0);
	
	for(int i=0;i<n-1;i++){
		std::cout<<i<<": "<<std::fixed<<pp0.coefA[i][0]<<" "<<pp0.coefB[i][0]<<" "<<pp0.coefC[i][0]<<" "<<pp0.coefD[i][0]<<"\n";
		std::cout<<i<<": "<<std::fixed<<pp0.coefA[i][1]<<" "<<pp0.coefB[i][1]<<" "<<pp0.coefC[i][1]<<" "<<pp0.coefD[i][1]<<"\n";
		std::cout<<i<<": "<<std::fixed<<pp0.coefA[i][2]<<" "<<pp0.coefB[i][2]<<" "<<pp0.coefC[i][2]<<" "<<pp0.coefD[i][2]<<"\n";
	}
	std::cout<<"==========================================\n";
	struct PPform pp1;
	n=8;
	pp1.n=8;
	pp1.points=points;

	pp1.cond_begin=(double*)malloc(3*sizeof(double));
	pp1.cond_end=(double*)malloc(3*sizeof(double));
	
	pp1.cond_begin[0]=0;
	pp1.cond_begin[1]=0;
	pp1.cond_begin[2]=0;

	pp1.cond_end[0]=0;
	pp1.cond_end[1]=0;
	pp1.cond_end[2]=0;
	

	cspline3D_var(&pp1);
	for(int i=0;i<n-1;i++){
		std::cout<<i<<": "<<std::fixed<<pp1.coefA[i][0]<<" "<<pp1.coefB[i][0]<<" "<<pp1.coefC[i][0]<<" "<<pp1.coefD[i][0]<<"\n";
		std::cout<<i<<": "<<std::fixed<<pp1.coefA[i][1]<<" "<<pp1.coefB[i][1]<<" "<<pp1.coefC[i][1]<<" "<<pp1.coefD[i][1]<<"\n";
		std::cout<<i<<": "<<std::fixed<<pp1.coefA[i][2]<<" "<<pp1.coefB[i][2]<<" "<<pp1.coefC[i][2]<<" "<<pp1.coefD[i][2]<<"\n";
	}
}
