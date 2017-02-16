#ifndef RK4
#define RK4_3D
#include<math.h>
#include"ppform.h"
#include<fstream>
#include <iomanip>
class RK4
{
	public:
	//RK
	double* k[6];
	double **RKy;//4 order
	double *RKe;

	//polynomial
	double *pt;
	double **acc;
	double **vel;
	double **pos;
	
	//cotrol para
	double *dt;
	double *ct;
	int *interval;//num points in the order interval
	int totalpoints;
	double rat;// distance cm/time sec;
	double  totaltime;
	
	//coef
	double *A;
	double *B;
	double *C;
	double *D;

	class PPform *pp;
	double h;//step size
	double* kFunction(double* y,double* x,double step,double dh,int index);
	//double pfunction();
	//double initialTable();
	void RKrunner();
	void PLrunner();
	void saveData();
	void arrange(double h,double t,PPform *p);//t is total time
	RK4(){};
};

#endif

//RK4 table
//
