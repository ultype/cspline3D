#include<iostream>
#include<stdlib.h>
#include<math.h>

int main(int argc,char *argv[]){
	int row=9;
	int col=8;

	double divdif[row-1][3];
	double points[9][3];
	double matrix[row][8];
	double dt[row-1];
	double ct[row];
	double b[row][3];
	


	
	double sol[row][3];
	//double v_cond0[3]={0,0,0};
	//double v_cond1[3]={0,0,0};
	double c4[row-1][3];
	double c3[row-1][3];
	double coefA[row-1][3];
	double coefB[row-1][3];
	double coefC[row-1][3];
	double coefD[row-1][3];

	points[0][0]=0; points[0][1]=0; points[0][2]=1;
	points[1][0]=1; points[1][1]=0; points[1][2]=0;
	points[2][0]=1; points[2][1]=1; points[2][2]=1;
	points[3][0]=0; points[3][1]=2; points[3][2]=0;
	points[4][0]=-1; points[4][1]=1; points[4][2]=1;
	points[5][0]=-1; points[5][1]=0; points[5][2]=0;
	points[6][0]=0; points[6][1]=-1; points[6][2]=1;
	points[7][0]=0; points[7][1]=-2; points[7][2]=0;
	points[8][0]=0; points[8][1]=0; points[8][2]=1;
	
	//diff
	double diff[8][3];
	for(int i=0;i<8;i++){
		diff[i][0]=points[i+1][0]-points[i][0];
		diff[i][1]=points[i+1][1]-points[i][1];
		diff[i][2]=points[i+1][2]-points[i][2];
	}
	
	//init matrix and dt 
	for(int i=0;i<row;i++){
		for(int j=0;j<col;j++){
			matrix[i][j]=0.0;
		}
		dt[i]=0.0;
	}
	
	//divdif difference
	for(int i=0;i<row-1;i++){
		dt[i]=pow(pow(diff[i][0],2)+pow(diff[i][1],2)+pow(diff[i][2],2),0.25);	
		divdif[i][0]=diff[i][0]/dt[i];
		divdif[i][1]=diff[i][1]/dt[i];
		divdif[i][2]=diff[i][2]/dt[i];
	}
	std::cout<<"divdif==================\n";
	for(int i=0;i<row-1;i++){
		std::cout<<divdif[i][0]<<" "<<divdif[i][1]<<" "<<divdif[i][2]<<"\n";
	}

	//lower tri
	for(int i=1;i<row-1;i++){
		matrix[i][0]=dt[i];
	}
	
	//upper tri
	for(int i=0;i<row-2;i++){
		matrix[i+1][2]=dt[i];
	}

	//mid tri
	for(int i=1;i<row-1;i++){
		matrix[i][1]=2*(dt[i]+dt[i-1]);
	}

	/*for variative
	matrix[0][1]=2;
	matrix[0][2]=1;
	

	matrix[row-1][0]=1;
	matrix[row-1][1]=2;
	*/
		
	//for pariodic
	matrix[0][1]=1;
	matrix[0][6]=-1;
	matrix[0][7]=2*dt[row-2];

	matrix[1][7]=dt[row-2];

	matrix[row-1][0]=dt[0];
	matrix[row-1][1]=2*dt[0];
	//==========================
	std::cout<<"dt=====================\n";
	for(int i=0;i<row-1;i++){
		std::cout<<i<<": "<<dt[i]<<" ";
	}
	
	std::cout<<"\nmatrix=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<std::fixed<<matrix[i][0]<<" "<<matrix[i][1]<<" "<<matrix[i][2]<<" "<<matrix[i][6]<<" "<<matrix[i][7]<<"\n";
	}
	std::cout<<"=======================\n";
	for(int i=1;i<row-1;i++){
		matrix[i][3]=3*(dt[i-1]*divdif[i][0]+dt[i]*divdif[i-1][0]);
		matrix[i][4]=3*(dt[i-1]*divdif[i][1]+dt[i]*divdif[i-1][1]);
		matrix[i][5]=3*(dt[i-1]*divdif[i][2]+dt[i]*divdif[i-1][2]);
	}
	/* for var both ends condition 
	matrix[0][3]=3*divdif[0][0]-(dt[0]*v_cond0[0]/2);
	matrix[0][4]=3*divdif[0][1]-(dt[0]*v_cond0[1]/2);
	matrix[0][5]=3*divdif[0][2]-(dt[0]*v_cond0[2]/2);

	matrix[row-1][3]=3*divdif[row-2][0]+(dt[row-1]*v_cond1[0]/2);
	matrix[row-1][4]=3*divdif[row-2][1]+(dt[row-1]*v_cond1[1]/2);
	matrix[row-1][5]=3*divdif[row-2][2]+(dt[row-1]*v_cond1[2]/2);
	*/
	
	//for periodic conditions
	matrix[0][3]=0;
	matrix[0][4]=0;
	matrix[0][5]=0;

	matrix[row-1][3]=3*(dt[row-2]*divdif[0][0]+dt[0]*divdif[row-2][0]);
	matrix[row-1][4]=3*(dt[row-2]*divdif[0][1]+dt[0]*divdif[row-2][1]);
	matrix[row-1][5]=3*(dt[row-2]*divdif[0][2]+dt[0]*divdif[row-2][2]);


	std::cout<<"\nmatrix rhs=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<std::fixed<<matrix[i][3]<<" "<<matrix[i][4]<<" "<<matrix[i][5]<<"\n";
	}
	
	//[0] sequence0
	//[1] sequence1
	//[2] sequence2
	//[6] sequence3
	//[7] sequence4
	double m;	
	for(int i=0;i<row-2;i++){//perodic row-2 variationsl row-1
		m= - matrix[i+1][0]/matrix[i][1];
		matrix[i+1][0]=0;
		matrix[i+1][1]=matrix[i+1][1]+m*matrix[i][2];
		//opr on sequece 3
		matrix[i+1][6]=m*matrix[i][6];
		//opr on b
		matrix[i+1][3]=matrix[i+1][3]+m*matrix[i][3];
		matrix[i+1][4]=matrix[i+1][4]+m*matrix[i][4];
		matrix[i+1][5]=matrix[i+1][5]+m*matrix[i][5];
	}
	std::cout<<"\nmatrix step1=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<std::fixed<<matrix[i][0]<<" "<<matrix[i][1]<<" "<<matrix[i][2]<<" "<<matrix[i][6]<<" "<<matrix[i][7]<<"\n";
	}	

	//sequence 3 collide sequence 2
	matrix[row-2][2]=matrix[row-2][2]+matrix[row-2][6];
	
	std::cout<<"\nmatrix step2=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<std::fixed<<matrix[i][0]<<" "<<matrix[i][1]<<" "<<matrix[i][2]<<" "<<matrix[i][6]<<" "<<matrix[i][7]<<"\n";
	}	


	//sequence 4
	for(int i=0;i<row-2;i++){
		m= - matrix[i][7]/matrix[i][1];
		//matrix[i][7]=0;
		matrix[row-1][1]=matrix[row-1][1]+m*matrix[i][6];
		matrix[i+1][7]=matrix[i+1][7]+m*matrix[i][2];
		//opr on b
		matrix[row-1][3]=matrix[row-1][3]+m*matrix[i][3];
		matrix[row-1][4]=matrix[row-1][4]+m*matrix[i][4];
		matrix[row-1][5]=matrix[row-1][5]+m*matrix[i][5];
	}

	std::cout<<"\nmatrix step3=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<std::fixed<<matrix[i][0]<<" "<<matrix[i][1]<<" "<<matrix[i][2]<<" "<<matrix[i][6]<<" "<<matrix[i][7]<<"\n";
	}	
	//sequence 4 collide sequence 0
	matrix[row-1][0]=matrix[row-1][0]+matrix[row-2][7];
	
	std::cout<<"\nmatrix step4=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<std::fixed<<matrix[i][0]<<" "<<matrix[i][1]<<" "<<matrix[i][2]<<" "<<matrix[i][6]<<" "<<matrix[i][7]<<"\n";
	}	
	//last opr
	m= - matrix[row-1][0]/matrix[row-2][1];
	matrix[row-1][0]=0;
	matrix[row-1][1]=matrix[row-1][1]+m*matrix[row-2][2];
	//opr on b
	matrix[row-1][3]=matrix[row-1][3]+m*matrix[row-2][3];
	matrix[row-1][4]=matrix[row-1][4]+m*matrix[row-2][4];
	matrix[row-1][5]=matrix[row-1][5]+m*matrix[row-2][5];



	std::cout<<"\nmatrix step5=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<std::fixed<<matrix[i][0]<<" "<<matrix[i][1]<<" "<<matrix[i][2]<<" "<<matrix[i][6]<<" "<<matrix[i][7]<<"\n";
	}
	std::cout<<"=======================\n";
	// solve c sol = b
	sol[row-1][0]=matrix[row-1][3]/matrix[row-1][1];
	sol[row-1][1]=matrix[row-1][4]/matrix[row-1][1];
	sol[row-1][2]=matrix[row-1][5]/matrix[row-1][1];
	
	sol[row-2][0]=(matrix[row-2][3]-matrix[row-2][2]*sol[row-1][0])/matrix[row-2][1];
	sol[row-2][1]=(matrix[row-2][4]-matrix[row-2][2]*sol[row-1][1])/matrix[row-2][1];
	sol[row-2][2]=(matrix[row-2][5]-matrix[row-2][2]*sol[row-1][2])/matrix[row-2][1];
	for(int i=row-3;i>=0;i--){
		sol[i][0]=(matrix[i][3]-matrix[i][2]*sol[i+1][0]-matrix[i][6]*sol[row-1][0])/matrix[i][1];
		sol[i][1]=(matrix[i][4]-matrix[i][2]*sol[i+1][1]-matrix[i][6]*sol[row-1][1])/matrix[i][1];
		sol[i][2]=(matrix[i][5]-matrix[i][2]*sol[i+1][2]-matrix[i][6]*sol[row-1][2])/matrix[i][1];
	}
	std::cout<<"sol:========================\n";
	for(int i=0;i<row;i++){
		std::cout<<std::fixed<<sol[i][0]<<" "<<sol[i][1]<<" "<<sol[i][2]<<"\n";
	}
	
	for(int i=0;i<row-1;i++){
		c4[i][0]=(sol[i][0]+sol[i+1][0]-2*divdif[i][0])/dt[i];
		c4[i][1]=(sol[i][1]+sol[i+1][1]-2*divdif[i][1])/dt[i];
		c4[i][2]=(sol[i][2]+sol[i+1][2]-2*divdif[i][2])/dt[i];

		c3[i][0]=(divdif[i][0]-sol[i][0])/dt[i]-c4[i][0];	
		c3[i][1]=(divdif[i][1]-sol[i][1])/dt[i]-c4[i][1];
		c3[i][2]=(divdif[i][2]-sol[i][2])/dt[i]-c4[i][2];
	}

	for(int i=0;i<row-1;i++){
		coefA[i][0]=c4[i][0]/dt[i];
		coefA[i][1]=c4[i][1]/dt[i];
		coefA[i][2]=c4[i][2]/dt[i];

		coefB[i][0]=c3[i][0];
		coefB[i][1]=c3[i][1];
		coefB[i][2]=c3[i][2];

		coefC[i][0]=sol[i][0];
		coefC[i][1]=sol[i][1];
		coefC[i][2]=sol[i][2];

		coefD[i][0]=points[i][0];
		coefD[i][1]=points[i][1];
		coefD[i][2]=points[i][2];
	}
	
	for(int i=0;i<row-1;i++){
		std::cout<<i<<": "<<coefA[i][0]<<" "<<coefB[i][0]<<" "<<coefC[i][0]<<" "<<coefD[i][0]<<"\n";
		std::cout<<i<<": "<<coefA[i][1]<<" "<<coefB[i][1]<<" "<<coefC[i][1]<<" "<<coefD[i][1]<<"\n";
		std::cout<<i<<": "<<coefA[i][2]<<" "<<coefB[i][2]<<" "<<coefC[i][2]<<" "<<coefD[i][2]<<"\n";
	}

	/*for periodic=================================================================
	// 1-st row
	matrix[0][0]=0.0;
	matrix[0][1]=1.0;
	matrix[0][2]=0.0;
	matrix[0][3]=-1.0;
	matrix[0][4]=2*dt[row-2];//n-1 th
	//matrix[0][5]=
	
	//2nd row
	matrix[1][4]=dt[row-1];
	

	//2-st~n-1 row
	matrix[row-2][0]=dt[0];
	matrix[row-1][1]=2*dt[0];

	matrix[row-2][3]=dt[row-2];
	matrix[row-1][3]=2*(dt[row-1]+dt[row-2]);

	matrix[row-2][4]=matrix[row-2][4]=dt[0];
	matrix[row-1][4]=matrix[row-1][4]=2*dt[0];


	int m;
	for(int i=0;i<row-2;row++){
		m= - matrix[i+1][0]/matrix[i][1];
		
		matrix[i+1][0]=0;
		matrix[i+1][1]=matrix[i+1][1]+m*matrix[i][2];
		matrix[i+1][3]=matrix[i+1][3]+m*matrix[i][3];
		//matrix[i+1][5]=matrix[i+1][5]+m*matrix[i][5];
	}
	matrix[row-2][2]=matrix[row-2][3];//last row shift
	
	//last row opr
	m= - matrix[row-1][0]/matrix[row-2][1];
	matrix[row-1][0]=0;
	matrix[row-1][1]=matrix[row-1][1]+m*matrix[row-2][2];
	//matrix[row-1][5]=matrix[row-1][5]+m*matrix[row-2][5];

	
	for(int i=0;i<row-1;i++){
		m= - matrix[i][4]/matrix[i][1];
		matrix[i+1][4]=matrix[i+1][4]+m*matrix[i][2];
		matrix[row-1][1]=matrix[row-1][1]+m*matrix[i][3];
		//matrix[row-1][5]=matrix[row-1][5]+m*matrix[i][5];
	}
	*/

}
