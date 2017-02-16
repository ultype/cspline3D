#ifndef PPFORM
#define PPFORM
#include<iostream>
#include<math.h>
#include<vector>
#include<stdlib.h>


class PPform
{
	public:
	bool periodic;
	double *points;
	int size;
	double *coefA;
	double *coefB;
	double *coefC;
	double *coefD;
	double *breaks;//ct cumulative time
	double *dt;
	double cond_begin[3];
	double cond_end[3];
	//===========================================
	PPform();
	PPform(std::vector< std::vector<double> > &vec);
	void loadPoints(std::vector< std::vector<double> > &vec);
	void loadConds(double *begin,double *end,bool periodic);
	void cspline3D_periodic();
	void cspline3D_var();
	void showCoefs();
};


#endif



