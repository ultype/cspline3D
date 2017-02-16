#include<iostream>
#include<stdlib.h>
#include<math.h>

int main(int argc,char *argv[]){
	double divided[8][3];
	double points[9][3];
	double matrix[8][8];
	double dt[7];
	double ct[8];
	double b[8][3];
	

	int row=8;
	int col=6;
	
	double sol[row][3];
	double v_cond0[3]={0,0,0};
	double v_cond1[3]={0,0,0};
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
	double pos[8][3];
	for(int i=0;i<8;i++){
		pos[i][0]=points[i+1][0]-points[i][0];
		pos[i][1]=points[i+1][1]-points[i][1];
		pos[i][2]=points[i+1][2]-points[i][2];
	}
	

	for(int i=0;i<row;i++){
		for(int j=0;j<col;j++){
			matrix[i][j]=0.0;
		}
		dt[i]=0.0;
	}
	
	//divided difference
	for(int i=0;i<row-1;i++){
		dt[i]=pow(pow(pos[i][0],2)+pow(pos[i][1],2)+pow(pos[i][2],2),0.25);	
		divided[i][0]=pos[i][0]/dt[i];
		divided[i][1]=pos[i][1]/dt[i];
		divided[i][2]=pos[i][2]/dt[i];
	}
	std::cout<<"divided==================\n";
	for(int i=0;i<row-1;i++){
		std::cout<<divided[i][0]<<" "<<divided[i][1]<<" "<<divided[i][2]<<"\n";
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
	//for var
	matrix[0][1]=2;
	matrix[0][2]=1;

	matrix[row-1][0]=1;
	matrix[row-1][1]=2;
	std::cout<<"dt=====================\n";
	for(int i=0;i<row-1;i++){
		std::cout<<i<<": "<<dt[i]<<" ";
	}
	
	std::cout<<"\nmatrix=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<matrix[i][0]<<" "<<matrix[i][1]<<" "<<matrix[i][2]<<"\n";
	}
	std::cout<<"=======================\n";
	for(int i=1;i<row-1;i++){
		matrix[i][3]=3*(dt[i-1]*divided[i][0]+dt[i]*divided[i-1][0]);
		matrix[i][4]=3*(dt[i-1]*divided[i][1]+dt[i]*divided[i-1][1]);
		matrix[i][5]=3*(dt[i-1]*divided[i][2]+dt[i]*divided[i-1][2]);
	}

	matrix[0][3]=3*divided[0][0]-(dt[0]*v_cond0[0]/2);
	matrix[0][4]=3*divided[0][1]-(dt[0]*v_cond0[1]/2);
	matrix[0][5]=3*divided[0][2]-(dt[0]*v_cond0[2]/2);

	matrix[row-1][3]=3*divided[row-2][0]+(dt[row-1]*v_cond1[0]/2);
	matrix[row-1][4]=3*divided[row-2][1]+(dt[row-1]*v_cond1[1]/2);
	matrix[row-1][5]=3*divided[row-2][2]+(dt[row-1]*v_cond1[2]/2);

	std::cout<<"\nmatrix rhs=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<matrix[i][3]<<" "<<matrix[i][4]<<" "<<matrix[i][5]<<"\n";
	}
	

	double m;	
	for(int i=0;i<row-1;i++){
		m= - matrix[i+1][0]/matrix[i][1];
		
		matrix[i+1][0]=0;
		matrix[i+1][1]=matrix[i+1][1]+m*matrix[i][2];
		
		matrix[i+1][3]=matrix[i+1][3]+m*matrix[i][3];
		matrix[i+1][4]=matrix[i+1][4]+m*matrix[i][4];
		matrix[i+1][5]=matrix[i+1][5]+m*matrix[i][5];
	}
	std::cout<<"\nmatrix=================\n";
	for(int i=0;i<row;i++){
		std::cout<<i<<": "<<matrix[i][0]<<" "<<matrix[i][1]<<" "<<matrix[i][2]<<"\n";
	}
	std::cout<<"=======================\n";
	// solve c sol = b
	sol[row-1][0]=matrix[row-1][3]/matrix[row-1][1];
	sol[row-1][1]=matrix[row-1][4]/matrix[row-1][1];
	sol[row-1][2]=matrix[row-1][5]/matrix[row-1][1];
	for(int i=row-2;i>=0;i--){
		sol[i][0]=(matrix[i][3]-matrix[i][2]*sol[i+1][0])/matrix[i][1];
		sol[i][1]=(matrix[i][4]-matrix[i][2]*sol[i+1][1])/matrix[i][1];
		sol[i][2]=(matrix[i][5]-matrix[i][2]*sol[i+1][2])/matrix[i][1];
	}
	std::cout<<"sol:========================\n";
	for(int i=0;i<row;i++){
		std::cout<<sol[i][0]<<" "<<sol[i][1]<<" "<<sol[i][2]<<"\n";
	}
	
	for(int i=0;i<row-1;i++){
		c4[i][0]=(sol[i][0]+sol[i+1][0]-2*divided[i][0])/dt[i];
		c4[i][1]=(sol[i][1]+sol[i+1][1]-2*divided[i][1])/dt[i];
		c4[i][2]=(sol[i][2]+sol[i+1][2]-2*divided[i][2])/dt[i];

		c3[i][0]=(divided[i][0]-sol[i][0])/dt[i]-c4[i][0];	
		c3[i][1]=(divided[i][1]-sol[i][1])/dt[i]-c4[i][1];
		c3[i][2]=(divided[i][2]-sol[i][2])/dt[i]-c4[i][2];
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
