
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstring>
using namespace std;


int main(int argc,char* argv[]){

// Euler algorithm for the motion of a particle subject to an external
// force f(x) = -x. An equal time step DT is used. The position,
// velocity and energy of the particle are written out.

	if(argc<4){
		cout << "<dT (in unit of pi)> <Totstep> <print skip>\n";
		exit(1);
	}

	unsigned int N       = atoi(argv[2]);
	const int print_step = atoi(argv[3]);
	double dT            = atof(argv[1])*M_PI;
	vector<double> T,V,X,E;
	T.resize(N);
	V = X = E = T;

	//Assign constants, initial position, and initial velocity
	//dT = 0.0025*M_PI; 
	V[0] = 1.0;
	X[0] = 0.0;
	T[0] = 0.0;
	
	//layout :
	printf("# dT = %8.8lf\n",dT/M_PI);
	printf("# Vi = %8.8lf\n",V[0]);
	printf("# Xi = %8.8lf\n",X[0]);
	printf("# tot_Step = %7d\n",N);
	printf("# print_skip = %7d\n",print_step);

	//Recursion for position and velocity at later time
	for(unsigned int i=0;i<N-1;i++){
		T[i+1] = (i+1)*dT;
		X[i+1] = X[i] + V[i]*dT;
		V[i+1] = V[i] - X[i]*dT;
	}
	for(unsigned int i=0;i<N;i++){
		E[i] = 0.5*(pow(V[i],2) + pow(X[i],2));
	}
	
	//Write the position, velocity and energy every IN steps
	printf("#%20s %20s %20s %20s\n","t","vel","pos","E");
	for(unsigned int i=0;i<N;i+=print_step){
		printf("%12.8lf %12.8lf %12.8lf %12.8lf\n",i,T[i]/M_PI,X[i],V[i],E[i]);
	}


	return 0;
}

