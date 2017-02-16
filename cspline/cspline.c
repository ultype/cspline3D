#include<stdio.h>
#include<stdlib.h>
#include<math.h>
struct PPform{
	double **points;
	int n;
	double **coefA;
	double **coefB;
	double **coefC;
	double **coefD;
	double *breaks;//ct
	double *dt;
	double *cond_begin;
	double *cond_end;
};

void cspline3D_periodic(struct PPform* pp){
	int row=pp->n;
	double **points=pp->points;
	if(points==NULL) return;
	int i,j;
	int row_1=row-1;
	int row_2=row-2;
	int row_3=row-3;
	

	double divdif[row_1][3];
	
	double dt[row_1];
	double ct[row];
	double matrix[row][8];//for periodic
	double b[row][3];

	double sol[row][3];
	
	
	//initial matrix
	for(i=0;i<row;i++){
		for(j=0;j<8;j++){
			matrix[i][j]=0.0;
		}
	}
	//divdif and dt
	ct[0]=0;
	double diff[3];
	for(i=0;i<row_1;i++){
		j=i+1;
		diff[0]=points[j][0]-points[i][0];
		diff[1]=points[j][1]-points[i][1];
		diff[2]=points[j][2]-points[i][2];

		dt[i]=pow(pow(diff[0],2)+pow(diff[1],2)+pow(diff[2],2),0.25);
		ct[i+1]=ct[i]+dt[i];


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
	
	//initial new coefs
	pp->coefA=(double**)malloc(row_1*sizeof(double*));
	pp->coefB=(double**)malloc(row_1*sizeof(double*));
	pp->coefC=(double**)malloc(row_1*sizeof(double*));
	pp->coefD=(double**)malloc(row_1*sizeof(double*));
	pp->breaks=(double*)malloc(row*sizeof(double));
	pp->dt=(double*)malloc(row_1*sizeof(double));

	for(i=0;i<row_1;i++){
		pp->coefA[i]=(double*)malloc(3*sizeof(double));
		pp->coefB[i]=(double*)malloc(3*sizeof(double));
		pp->coefC[i]=(double*)malloc(3*sizeof(double));
		pp->coefD[i]=(double*)malloc(3*sizeof(double));
		pp->dt[i]=dt[i];
		pp->breaks[i]=ct[i];
	}
	ct[row_1]=ct[row_2]+dt[row_2];

	//coef
	double c4[3];
	for(int i=0;i<row_1;i++){
		c4[0]=(sol[i][0]+sol[i+1][0]-2*divdif[i][0])/dt[i];
		c4[1]=(sol[i][1]+sol[i+1][1]-2*divdif[i][1])/dt[i];
		c4[2]=(sol[i][2]+sol[i+1][2]-2*divdif[i][2])/dt[i];

		pp->coefA[i][0]=c4[0]/dt[i];
		pp->coefA[i][1]=c4[1]/dt[i];
		pp->coefA[i][2]=c4[2]/dt[i];


		pp->coefB[i][0]=(divdif[i][0]-sol[i][0])/dt[i]-c4[0];	
		pp->coefB[i][1]=(divdif[i][1]-sol[i][1])/dt[i]-c4[1];
		pp->coefB[i][2]=(divdif[i][2]-sol[i][2])/dt[i]-c4[2];
		
		
		pp->coefC[i][0]=sol[i][0];
		pp->coefC[i][1]=sol[i][1];
		pp->coefC[i][2]=sol[i][2];

		pp->coefD[i][0]=points[i][0];
		pp->coefD[i][1]=points[i][1];
		pp->coefD[i][2]=points[i][2];
	}
}

void cspline3D_var(struct PPform* pp){
	int row=pp->n;
	double **points=pp->points;
	if(points==NULL) return;
	int i,j;
	int row_1=row-1;
	int row_2=row-2;
	int row_3=row-3;
	

	double divdif[row_1][3];
	
	double dt[row_1];
	double ct[row];
	double matrix[row][8];//for periodic
	double b[row][3];

	double sol[row][3];
	
	
	//initial matrix
	for(i=0;i<row;i++){
		for(j=0;j<8;j++){
			matrix[i][j]=0.0;
		}
	}
	//divdif and dt
	ct[0]=0;
	double diff[3];
	for(i=0;i<row_1;i++){
		j=i+1;
		diff[0]=points[j][0]-points[i][0];
		diff[1]=points[j][1]-points[i][1];
		diff[2]=points[j][2]-points[i][2];

		dt[i]=pow(pow(diff[0],2)+pow(diff[1],2)+pow(diff[2],2),0.25);
		ct[i+1]=ct[i]+dt[i];


		divdif[i][0]=diff[0]/dt[i];
		divdif[i][1]=diff[1]/dt[i];
		divdif[i][2]=diff[2]/dt[i];
	}
//	printf("divdif ok\n");
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
//	printf("matrix initial ok\n");
	//compute b	
	for(i=1;i<row_1;i++){
		j=i-1;
		matrix[i][5]=3*(dt[j]*divdif[i][0]+dt[i]*divdif[j][0]);
		matrix[i][6]=3*(dt[j]*divdif[i][1]+dt[i]*divdif[j][1]);
		matrix[i][7]=3*(dt[j]*divdif[i][2]+dt[i]*divdif[j][2]);
	}
//	printf("line 281 ok\n");
	//var condiction for b
	matrix[0][5]=3*divdif[0][0]-(dt[0]*pp->cond_begin[0]/2);
	matrix[0][6]=3*divdif[0][1]-(dt[0]*pp->cond_begin[1]/2);
	matrix[0][7]=3*divdif[0][2]-(dt[0]*pp->cond_begin[2]/2);

	matrix[row_1][5]=3*divdif[row_2][0]+(dt[row_1]*pp->cond_end[0]/2);
	matrix[row_1][6]=3*divdif[row_2][1]+(dt[row_1]*pp->cond_end[1]/2);
	matrix[row_1][7]=3*divdif[row_2][2]+(dt[row_1]*pp->cond_end[2]/2);
	
//	printf("vector b ok\n");
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
	//printf("matrix solve ok\n");
	//initial new coefs
	pp->coefA=(double**)malloc(row_1*sizeof(double*));
	pp->coefB=(double**)malloc(row_1*sizeof(double*));
	pp->coefC=(double**)malloc(row_1*sizeof(double*));
	pp->coefD=(double**)malloc(row_1*sizeof(double*));
	pp->breaks=(double*)malloc(row*sizeof(double));
	pp->dt=(double*)malloc(row_1*sizeof(double));

	for(i=0;i<row_1;i++){
		pp->coefA[i]=(double*)malloc(3*sizeof(double));
		pp->coefB[i]=(double*)malloc(3*sizeof(double));
		pp->coefC[i]=(double*)malloc(3*sizeof(double));
		pp->coefD[i]=(double*)malloc(3*sizeof(double));
		pp->dt[i]=dt[i];
		pp->breaks[i]=ct[i];
	}
	ct[row_1]=ct[row_2]+dt[row_2];

	//coef
	double c4[3];
	for(int i=0;i<row_1;i++){
		c4[0]=(sol[i][0]+sol[i+1][0]-2*divdif[i][0])/dt[i];
		c4[1]=(sol[i][1]+sol[i+1][1]-2*divdif[i][1])/dt[i];
		c4[2]=(sol[i][2]+sol[i+1][2]-2*divdif[i][2])/dt[i];

		pp->coefA[i][0]=c4[0]/dt[i];
		pp->coefA[i][1]=c4[1]/dt[i];
		pp->coefA[i][2]=c4[2]/dt[i];


		pp->coefB[i][0]=(divdif[i][0]-sol[i][0])/dt[i]-c4[0];	
		pp->coefB[i][1]=(divdif[i][1]-sol[i][1])/dt[i]-c4[1];
		pp->coefB[i][2]=(divdif[i][2]-sol[i][2])/dt[i]-c4[2];
		
		
		pp->coefC[i][0]=sol[i][0];
		pp->coefC[i][1]=sol[i][1];
		pp->coefC[i][2]=sol[i][2];

		pp->coefD[i][0]=points[i][0];
		pp->coefD[i][1]=points[i][1];
		pp->coefD[i][2]=points[i][2];
	}
	//printf("coef ok\n");
}

