#include"ppform.h"
PPform::PPform(std::vector<std::vector<double> > &vec)
{
	loadPoints(vec);
	if(periodic){
		cspline3D_periodic();	
	}else{
		cspline3D_var();
	}
}
PPform::PPform(){
	periodic=false;
	points=NULL;
	size=0;
	coefA=NULL;
	coefB=NULL;
	coefC=NULL;
	coefD=NULL;
	breaks=NULL;//ct
	dt=NULL;
}
void PPform::loadPoints(std::vector< std::vector<double> > &vec){
	if(vec.empty()) return;
	this->size=vec.size();
	this->periodic=false;
	if(vec.front()[0]==vec.back()[0] &&
	   vec.front()[1]==vec.back()[1] &&
	   vec.front()[2]==vec.back()[2]) periodic=true;
	this->points=(double*)malloc(3*this->size*sizeof(double*));
	int k;
	for(int i=0;i<this->size;i++){
		k=3*i;
		points[k]=vec[i][0];
		points[k+1]=vec[i][1];
		points[k+2]=vec[i][2];
	}
}

void PPform::loadConds(double *begin,double *end,bool periodic){
	if(periodic) this->periodic=true;
	if(begin==NULL || end==NULL) return;
	this->cond_begin[0]=begin[0];	
	this->cond_begin[1]=begin[1];
	this->cond_begin[2]=begin[2];

	this->cond_end[0]=end[0];	
	this->cond_end[1]=end[1];	
	this->cond_end[2]=end[2];	
}

void PPform::cspline3D_periodic(){
	int row=size;
	if(points==NULL) return;
	int i,j,k1,k2,k3;
	int row_1=row-1;
	int row_2=row-2;
	int row_3=row-3;

	double divdif[row_1][3];	
	double matrix[row][8];//for periodic
	double b[row][3];

	double sol[row][3];
	

	//initial matrix
	for(i=0;i<row;i++){
		for(j=0;j<8;j++){
			matrix[i][j]=0.0;
		}
	}
		//printf("matrix solve ok\n");
	//initial new coefs
	coefA=(double*)malloc(3*row_1*sizeof(double));
	coefB=(double*)malloc(3*row_1*sizeof(double));
	coefC=(double*)malloc(3*row_1*sizeof(double));
	coefD=(double*)malloc(3*row_1*sizeof(double));
	breaks=(double*)malloc(row*sizeof(double));
	dt=(double*)malloc(row_1*sizeof(double));
	//divdif and dt
	breaks[0]=0;
	double diff[3];
	for(i=0;i<row_1;i++){
		j=i+1;
		k1=3*i;
		k2=3*j;
		diff[0]=points[k2]-points[k1];
		diff[1]=points[k2+1]-points[k1+1];
		diff[2]=points[k2+2]-points[k1+2];

		dt[i]=pow(pow(diff[0],2)+pow(diff[1],2)+pow(diff[2],2),0.25);
		breaks[j]=breaks[i]+dt[i];


		divdif[i][0]=diff[0]/dt[i];
		divdif[i][1]=diff[1]/dt[i];
		divdif[i][2]=diff[2]/dt[i];
	}

	//lower tri seq0
	for(i=1;i<row_1;i++){
		matrix[i][0]=dt[i];
	}

	//upper tri seq2
	for(i=0;i<row_2;i++){
		matrix[i+1][2]=dt[i];
	}

	//mid tri seq1
	for(i=1;i<row_1;i++){
		matrix[i][1]=2*(dt[i]+dt[i-1]);
	}
	
	//for periodic condition matrix
	matrix[0][1]=1;
	matrix[0][3]=-1;
	matrix[0][4]=2*dt[row_2];
	
	matrix[1][4]=dt[row_2];

	matrix[row_1][0]=dt[0];
	matrix[row_1][1]=2*dt[0];
	
	//compute b	
	for(i=1;i<row_1;i++){
		j=i-1;
		matrix[i][5]=3*(dt[j]*divdif[i][0]+dt[i]*divdif[j][0]);
		matrix[i][6]=3*(dt[j]*divdif[i][1]+dt[i]*divdif[j][1]);
		matrix[i][7]=3*(dt[j]*divdif[i][2]+dt[i]*divdif[j][2]);
	}

	//periodic condiction for b
	matrix[0][5]=0;
	matrix[0][6]=0;
	matrix[0][7]=0;

	matrix[row_1][5]=3*(dt[row_2]*divdif[0][0]+dt[0]*divdif[row_2][0]);
	matrix[row_1][6]=3*(dt[row_2]*divdif[0][1]+dt[0]*divdif[row_2][1]);
	matrix[row_1][7]=3*(dt[row_2]*divdif[0][2]+dt[0]*divdif[row_2][2]);

	//solve matrix 
	double m;
	for(i=0;i<row_2;i++){//perodic row-2 variationsl row-1
		j=i+1;
		m= - matrix[j][0]/matrix[i][1];
		matrix[j][0]=0;
		matrix[j][1]=matrix[j][1]+m*matrix[i][2];
		//opr on sequece 3
		matrix[j][3]=m*matrix[i][3];
		//opr on b
		matrix[j][5]=matrix[j][5]+m*matrix[i][5];
		matrix[j][6]=matrix[j][6]+m*matrix[i][6];
		matrix[j][7]=matrix[j][7]+m*matrix[i][7];
	}

	//sequence 3 collide sequence 2
	matrix[row_2][2]=matrix[row_2][2]+matrix[row_2][3];

	//sequence 4
	for(i=0;i<row_2;i++){
		j=i+1;
		m= - matrix[i][4]/matrix[i][1];
		//matrix[i][7]=0;
		matrix[row_1][1]=matrix[row_1][1]+m*matrix[i][3];
		matrix[j][4]=matrix[j][4]+m*matrix[i][2];
		//opr on b
		matrix[row_1][5]=matrix[row_1][5]+m*matrix[i][5];
		matrix[row_1][6]=matrix[row_1][6]+m*matrix[i][6];
		matrix[row_1][7]=matrix[row_1][7]+m*matrix[i][7];
	}
	
	//sequence 4 collide sequence 0
	matrix[row_1][0]=matrix[row_1][0]+matrix[row_2][4];
	
	//last opr
	m= - matrix[row_1][0]/matrix[row_2][1];
	matrix[row_1][0]=0;
	matrix[row_1][1]=matrix[row_1][1]+m*matrix[row_2][2];
	
	//last opr on b
	matrix[row_1][5]=matrix[row_1][5]+m*matrix[row_2][5];
	matrix[row_1][6]=matrix[row_1][6]+m*matrix[row_2][6];
	matrix[row_1][7]=matrix[row_1][7]+m*matrix[row_2][7];
	
	//solve C*sol=b
	sol[row_1][0]=matrix[row_1][5]/matrix[row_1][1];
	sol[row_1][1]=matrix[row_1][6]/matrix[row_1][1];
	sol[row_1][2]=matrix[row_1][7]/matrix[row_1][1];
	
	sol[row_2][0]=(matrix[row_2][5]-matrix[row_2][2]*sol[row_1][0])/matrix[row_2][1];
	sol[row_2][1]=(matrix[row_2][6]-matrix[row_2][2]*sol[row_1][1])/matrix[row_2][1];
	sol[row_2][2]=(matrix[row_2][7]-matrix[row_2][2]*sol[row_1][2])/matrix[row_2][1];
	

	for(i=row_3;i>=0;i--){
		j=i+1;
		sol[i][0]=(matrix[i][5]-matrix[i][2]*sol[j][0]-matrix[i][3]*sol[row_1][0])/matrix[i][1];
		sol[i][1]=(matrix[i][6]-matrix[i][2]*sol[j][1]-matrix[i][3]*sol[row_1][1])/matrix[i][1];
		sol[i][2]=(matrix[i][7]-matrix[i][2]*sol[j][2]-matrix[i][3]*sol[row_1][2])/matrix[i][1];
	}
	
	//coef
	double c4[3];
	for(int i=0;i<row_1;i++){
		c4[0]=(sol[i][0]+sol[i+1][0]-2*divdif[i][0])/dt[i];
		c4[1]=(sol[i][1]+sol[i+1][1]-2*divdif[i][1])/dt[i];
		c4[2]=(sol[i][2]+sol[i+1][2]-2*divdif[i][2])/dt[i];
		
		k1=3*i;
		k2=k1+1;
		k3=k2+1;
		coefA[k1]=c4[0]/dt[i];
		coefA[k2]=c4[1]/dt[i];
		coefA[k3]=c4[2]/dt[i];


		coefB[k1]=(divdif[i][0]-sol[i][0])/dt[i]-c4[0];	
		coefB[k2]=(divdif[i][1]-sol[i][1])/dt[i]-c4[1];
		coefB[k3]=(divdif[i][2]-sol[i][2])/dt[i]-c4[2];
		
		
		coefC[k1]=sol[i][0];
		coefC[k2]=sol[i][1];
		coefC[k3]=sol[i][2];

		coefD[k1]=points[k1];
		coefD[k2]=points[k2];
		coefD[k3]=points[k3];
	}	
}

void PPform::cspline3D_var(){
	int row=size;
	if(points==NULL) return;
	int i,j,k1,k2,k3;
	int row_1=row-1;
	int row_2=row-2;
	int row_3=row-3;

	double divdif[row_1][3];	
	double matrix[row][8];//for periodic
	double b[row][3];

	double sol[row][3];
	

	//initial matrix
	for(i=0;i<row;i++){
		for(j=0;j<8;j++){
			matrix[i][j]=0.0;
		}
	}
		//printf("matrix solve ok\n");
	//initial new coefs
	coefA=(double*)malloc(3*row_1*sizeof(double));
	coefB=(double*)malloc(3*row_1*sizeof(double));
	coefC=(double*)malloc(3*row_1*sizeof(double));
	coefD=(double*)malloc(3*row_1*sizeof(double));
	breaks=(double*)malloc(row*sizeof(double));
	dt=(double*)malloc(row_1*sizeof(double));
	//divdif and dt
	breaks[0]=0;
	double diff[3];
	for(i=0;i<row_1;i++){
		j=i+1;
		k1=3*i;
		k2=3*j;
		diff[0]=points[k2]-points[k1];
		diff[1]=points[k2+1]-points[k1+1];
		diff[2]=points[k2+2]-points[k1+2];

		dt[i]=pow(pow(diff[0],2)+pow(diff[1],2)+pow(diff[2],2),0.25);
		breaks[j]=breaks[i]+dt[i];


		divdif[i][0]=diff[0]/dt[i];
		divdif[i][1]=diff[1]/dt[i];
		divdif[i][2]=diff[2]/dt[i];
	}

	//lower tri seq0
	for(i=1;i<row_1;i++){
		matrix[i][0]=dt[i];
	}

	//upper tri seq2
	for(i=0;i<row_2;i++){
		matrix[i+1][2]=dt[i];
	}

	//mid tri seq1
	for(i=1;i<row_1;i++){
		matrix[i][1]=2*(dt[i]+dt[i-1]);
	}
	
	//for var condition matrix
	matrix[0][1]=2;
	matrix[0][2]=1;

	matrix[row_1][0]=1;
	matrix[row_1][1]=2;
	
	//compute b	
	for(i=1;i<row_1;i++){
		j=i-1;
		matrix[i][5]=3*(dt[j]*divdif[i][0]+dt[i]*divdif[j][0]);
		matrix[i][6]=3*(dt[j]*divdif[i][1]+dt[i]*divdif[j][1]);
		matrix[i][7]=3*(dt[j]*divdif[i][2]+dt[i]*divdif[j][2]);
	}
	
	//var condiction for b
	matrix[0][5]=3*divdif[0][0]-(dt[0]*cond_begin[0]/2);
	matrix[0][6]=3*divdif[0][1]-(dt[0]*cond_begin[1]/2);
	matrix[0][7]=3*divdif[0][2]-(dt[0]*cond_begin[2]/2);

	matrix[row_1][5]=3*divdif[row_2][0]+(dt[row_1]*cond_end[0]/2);
	matrix[row_1][6]=3*divdif[row_2][1]+(dt[row_1]*cond_end[1]/2);
	matrix[row_1][7]=3*divdif[row_2][2]+(dt[row_1]*cond_end[2]/2);

	//solve matrix 
	double m;
	for(i=0;i<row_1;i++){//perodic row-2 variationsl row-1
		j=i+1;
		m= - matrix[j][0]/matrix[i][1];
		matrix[j][0]=0;
		matrix[j][1]=matrix[j][1]+m*matrix[i][2];
		//opr on b
		matrix[j][5]=matrix[j][5]+m*matrix[i][5];
		matrix[j][6]=matrix[j][6]+m*matrix[i][6];
		matrix[j][7]=matrix[j][7]+m*matrix[i][7];
	}
//	printf("matrix operation ok\n");
	
	//solve C*sol=b
	sol[row_1][0]=matrix[row_1][5]/matrix[row_1][1];
	sol[row_1][1]=matrix[row_1][6]/matrix[row_1][1];
	sol[row_1][2]=matrix[row_1][7]/matrix[row_1][1];
	
	for(int i=row_2;i>=0;i--){
		j=i+1;
		sol[i][0]=(matrix[i][5]-matrix[i][2]*sol[j][0])/matrix[i][1];
		sol[i][1]=(matrix[i][6]-matrix[i][2]*sol[j][1])/matrix[i][1];
		sol[i][2]=(matrix[i][7]-matrix[i][2]*sol[j][2])/matrix[i][1];
	}
	

	
	//coef
	double c4[3];
	for(int i=0;i<row_1;i++){
		c4[0]=(sol[i][0]+sol[i+1][0]-2*divdif[i][0])/dt[i];
		c4[1]=(sol[i][1]+sol[i+1][1]-2*divdif[i][1])/dt[i];
		c4[2]=(sol[i][2]+sol[i+1][2]-2*divdif[i][2])/dt[i];
		
		k1=3*i;
		k2=k1+1;
		k3=k2+1;
		coefA[k1]=c4[0]/dt[i];
		coefA[k2]=c4[1]/dt[i];
		coefA[k3]=c4[2]/dt[i];


		coefB[k1]=(divdif[i][0]-sol[i][0])/dt[i]-c4[0];	
		coefB[k2]=(divdif[i][1]-sol[i][1])/dt[i]-c4[1];
		coefB[k3]=(divdif[i][2]-sol[i][2])/dt[i]-c4[2];
		
		
		coefC[k1]=sol[i][0];
		coefC[k2]=sol[i][1];
		coefC[k3]=sol[i][2];

		coefD[k1]=points[k1];
		coefD[k2]=points[k2];
		coefD[k3]=points[k3];
	}	
}

void PPform::showCoefs(){
	int k1,k2,k3;
	for(int i=0;i<size-1;i++){
		k1=3*i;
		k2=k1+1;
		k3=k2+1;
		std::cout<<std::fixed<<coefA[k1]<<" "<<coefB[k1]<<" "<<coefC[k1]<<" "<<coefD[k1]<<"\n";
		std::cout<<std::fixed<<coefA[k2]<<" "<<coefB[k2]<<" "<<coefC[k2]<<" "<<coefD[k2]<<"\n";
		std::cout<<std::fixed<<coefA[k3]<<" "<<coefB[k3]<<" "<<coefC[k3]<<" "<<coefD[k3]<<"\n";
	}
}
