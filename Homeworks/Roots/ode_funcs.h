#include <gsl/gsl_vector.h>


void rkstep12g (
        void g(int n, double x, double* y, double* dydx,double E),
        int n,
        double x,
        double* yx,
        double h,
        double* yh,
        double* dy,
	double E
);

void rkstep12Tc (
        void f(int n, double x, double* y, double* dydx, double Tc, double Tr, int N),
        int n,
        double x,
        double* yx,
        double h,
        double* yh,
        double* dy,
        double Tc,
        double Tr,
        int N
);



void f(int n, double x, double* y, double* dydx,double Tc, double Tr,int N);

void g(int n, double x, double* y, double* dydx);



void driver ( 
        void f(int n, double x, double* y, double* dydx, double Tc, double Tr, int N),
        int n,
        double a,
        double* ya,
        double b,
        double h,
        double acc,
        double eps,
	double Tc,
	double Tr,
	int N,
	char* name
);

void driverg ( 
        void g(int n, double x, double* y, double* dydx,double E),
        int n,
        double a,
        double* ya,
        double b,
        double h,
        double acc,
        double eps,
        char* name,
	gsl_vector* res,
	double E
);







