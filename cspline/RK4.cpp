#include"RK4.h"


void RK4::arrange(double h,double t,PPform *p){
	this->pp=p;
	dt=pp->dt;
	ct=pp->breaks;
	A=pp->coefA;
	B=pp->coefB;
	C=pp->coefC;
	D=pp->coefD;

	int size=pp->size;
	int size_1=size-1;
	int size_2=size_1-1;

	this->h=h;
	totaltime=t;
	rat=ct[size-1]/t;
	interval=(int*)malloc(size_1*sizeof(int));
	
	this->totalpoints=size;
	for(int i=0;i<size_1;i++){
		interval[i]=floor(dt[i]/h);
		totalpoints+=interval[i];
	}
	//build space
	RKy=(double**)malloc(totalpoints*sizeof(double*));
	RKe=(double*)malloc(totalpoints*sizeof(double));
	pt=(double*)malloc(totalpoints*sizeof(double));//extend points
	
	for(int i=0;i<totalpoints;i++){
		RKy[i]=(double*)malloc(6*sizeof(double));
	}

	for(int i=0;i<4;i++){
		k[i]=(double*)malloc(6*sizeof(double));
	}
	//inpolate t points for p
	int q=0;
	double step=0;
	for(int i=0;i<size_1;i++){
		pt[q]=step;
		q++;
		for(int j=0;j<interval[i];j++){
			step+=h;
			pt[q]=step;
			q++;
		}	
		step=ct[i+1];
	}
	pt[totalpoints-1]=ct[size_1];
}

double* RK4::kFunction(double* y,double* x,double step,double dh,int index){
	int k1=index*3;
	int k2=k1+1;
	int k3=k2+1;
	double t=step+dh;
	//accelrate
	y[3]=2*B[k1]+6*A[k1]*t;
	y[4]=2*B[k2]+6*A[k2]*t;
	y[5]=2*B[k3]+6*A[k3]*t;
	
	//velocity
	y[0]=x[3]+dh*y[3];
	y[1]=x[4]+dh*y[4];
	y[2]=x[5]+dh*y[5];

	return y;
}
//RK4 table
/*
 *  0  |  0	0     0    0
 * 1/3 | 1/3	0     0	   0
 * 2/3 |-1/3	1     0    0
 *  1  |  1    -1     1    0
 *  --------------------------
 *     | 1/8   3/8   3/8  1/8
 */

void RK4::RKrunner(){

RKy[0][0]=D[0];
RKy[0][1]=D[1];
RKy[0][2]=D[2];
RKy[0][3]=C[0];
RKy[0][4]=C[1];
RKy[0][5]=C[2];

int index=0;
double step;
double dh;
int counter=interval[0];
int i_1;
int totalpoints_1=totalpoints-1;
for(int i=0;i<totalpoints_1;i++){
	i_1=i+1;
	dh=pt[i_1]-pt[i];
	step=pt[i]-ct[index];
	double temp[6];
	//RK
	kFunction(k[0],RKy[i],step,0,index);
	//================================================================
	temp[0]=RKy[i][0]+dh*k[0][0]/3;	
	temp[1]=RKy[i][1]+dh*k[0][1]/3;
	temp[2]=RKy[i][2]+dh*k[0][2]/3;
	temp[3]=RKy[i][3]+dh*k[0][3]/3;
	temp[4]=RKy[i][4]+dh*k[0][4]/3;
	temp[5]=RKy[i][5]+dh*k[0][5]/3;
	kFunction(k[1],temp,step,dh/3,index);
	//=================================================================
	temp[0]=RKy[i][0]+dh*(-k[0][0]/3+k[1][0]);	
	temp[1]=RKy[i][1]+dh*(-k[0][1]/3+k[1][1]);
	temp[2]=RKy[i][2]+dh*(-k[0][2]/3+k[1][2]);
	temp[3]=RKy[i][3]+dh*(-k[0][3]/3+k[1][3]);
	temp[4]=RKy[i][4]+dh*(-k[0][4]/3+k[1][4]);
	temp[5]=RKy[i][5]+dh*(-k[0][5]/3+k[1][5]);
	kFunction(k[2],temp,step,2*dh/3,index);
	//====================================================================
	temp[0]=RKy[i][0]+dh*(k[0][0]-k[1][0]+k[2][0]);	
	temp[1]=RKy[i][1]+dh*(k[0][1]-k[1][1]+k[2][1]);
	temp[2]=RKy[i][2]+dh*(k[0][2]-k[1][2]+k[2][2]);
	temp[3]=RKy[i][3]+dh*(k[0][3]-k[1][3]+k[2][3]);
	temp[4]=RKy[i][4]+dh*(k[0][4]-k[1][4]+k[2][4]);
	temp[5]=RKy[i][5]+dh*(k[0][5]-k[1][5]+k[2][5]);
	kFunction(k[3],temp,step,dh,index);

	//==========================================================================
	RKy[i_1][0]=RKy[i][0]+dh*(k[0][0]+3*k[1][0]+3*k[2][0]+k[3][0])/8; 
	RKy[i_1][1]=RKy[i][1]+dh*(k[0][1]+3*k[1][1]+3*k[2][1]+k[3][1])/8;  
	RKy[i_1][2]=RKy[i][2]+dh*(k[0][2]+3*k[1][2]+3*k[2][2]+k[3][2])/8; 
	RKy[i_1][3]=RKy[i][3]+dh*(k[0][3]+3*k[1][3]+3*k[2][3]+k[3][3])/8; 
	RKy[i_1][4]=RKy[i][4]+dh*(k[0][4]+3*k[1][4]+3*k[2][4]+k[3][4])/8;
	RKy[i_1][5]=RKy[i][5]+dh*(k[0][5]+3*k[1][5]+3*k[2][5]+k[3][5])/8;  
	//===================================================================================	
	if(counter==-1){
		index++; 
		counter=interval[index];
	}else{
		counter--;
	}

}
std::cout<<"RKrunner 196 ok\n";
}


void RK4::saveData(){
	std::ofstream ofs("outcome0.in",std::ofstream::out);
	if(!ofs.is_open()){ 
		std::cout<<"RK45 open file failed\n"; 
		return;
	}
	for(int i=0;i<totalpoints;i++){
	ofs<<std::setprecision(11)<<RKy[i][0]<<" "
				  <<RKy[i][1]<<" " 
				  <<RKy[i][2]<<" "
				  <<RKy[i][3]<<" "
				  <<RKy[i][4]<<" "
				  <<RKy[i][5]<<"\n";
	}
	ofs.close();
}
