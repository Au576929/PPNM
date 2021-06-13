#include <unistd.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h>

#include "lin_funcs.h"
#include "rf_funcs.h"
#include "ode_funcs.h"

int main(){

unsigned int rand_seed=3;

void func1(gsl_vector* X, gsl_vector* fX){
	gsl_blas_dcopy(X,fX);
}


int dim=1;
printf("Starting with the %i dimensional identity function:\n",dim);

gsl_vector* X1=gsl_vector_alloc(dim);
gsl_vector* fX1=gsl_vector_alloc(dim);
for (int i=0;i<dim;i++){
	gsl_vector_set(X1,i,rand_r(&rand_seed)*1.0/RAND_MAX*15);
}




printf("Start point:\n");
vector_print(X1);
printf("gives func value:\n");
func1(X1,fX1);
vector_print(fX1);

double eps=0.00005;
printf("running algoritm \n");
Newton(func1,X1,eps);

printf("found root::\n");
vector_print(X1);
func1(X1,fX1);
printf("function value of root:\n");
vector_print(fX1);


printf("used tolerance: %10f\n",eps);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void func2(gsl_vector* X, gsl_vector* fX){
        double x0=gsl_vector_get(X,0);
        double x1=gsl_vector_get(X,1);
        gsl_vector_set(fX,0,50*sin(x0)*cos(x1));
        gsl_vector_set(fX,1,50*cos(x0)*sin(x0));
}

dim=2;
printf("\n\n now onto the %i dimensional function:\n f(x,y)=(f1,f2)    f1=50*sin(x)*cos(y)    f2=50*cos(y)*sin(x)\n",dim);

gsl_vector* X=gsl_vector_alloc(dim);
gsl_vector* fX=gsl_vector_alloc(dim);


for (int i=0;i<dim;i++){
        gsl_vector_set(X,i,rand_r(&rand_seed)*1.0/RAND_MAX*0.5);
}




printf("Start point:\n");
vector_print(X);
printf("gives func value:\n");
func2(X,fX);
vector_print(fX);

eps=0.00005;
printf("running algoritm \n");

Newton(func2,X,eps);
printf("found root:\n");
vector_print(X);
func2(X,fX);
printf("function value of root:\n");
vector_print(fX);


printf("used tolerance: %10f\n",eps);



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void func3(gsl_vector* X, gsl_vector* fX){
        double x0=gsl_vector_get(X,0);
        double x1=gsl_vector_get(X,1);
        gsl_vector_set(fX,0,sin(cos(x0)+sin(x1)));
        gsl_vector_set(fX,1,cos(sin(x0)+cos(x1)));
}

printf("\n\n now onto the %i dimensional function:\n f(x,y)=(f1,f2)    f1=sin(cos(x)+sin(y))    f2=cos(sin(x)+cos(y))\n",dim);

for (int i=0;i<dim;i++){
        gsl_vector_set(X,i,rand_r(&rand_seed)*1.0/RAND_MAX*10);
}



printf("Start point:\n");
vector_print(X);
printf("gives func value:\n");
func3(X,fX);
vector_print(fX);

eps=0.00005;
printf("running algoritm \n");
Newton(func3,X,eps);
printf("found root::\n");
vector_print(X);
func3(X,fX);
printf("function value of root:\n");
vector_print(fX);


printf("used tolerance: %10f\n",eps);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void func4(gsl_vector* X, gsl_vector* fX){
        double x0=gsl_vector_get(X,0);
        double x1=gsl_vector_get(X,1);
        gsl_vector_set(fX,0,2*(1-x0)*(-1)+200*(x1-x0*x0)*(-2*x0));
        gsl_vector_set(fX,1,200*(x1-x0*x0));
}

printf("\n\n now onto the gradient of the Rosenbrock's valley function\n");

for (int i=0;i<dim;i++){
        gsl_vector_set(X,i,1+rand_r(&rand_seed)*1.0/RAND_MAX*0.001);
}




printf("Start point:\n");
vector_print(X);
printf("gives func value:\n");
func4(X,fX);
vector_print(fX);

eps=0.000005;
printf("running algoritm \n");

Newton(func4,X,eps);
printf("found root::\n");
vector_print(X);
func4(X,fX);
printf("function value of root:\n");
vector_print(fX);


printf("used tolerance: %10f\n",eps);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void g(int n, double x, double* y, double* dydx, double E){

	dydx[1]=(E+1.0/x)*y[0]*(-2);
	dydx[0]=y[1];
}



void func8newton(gsl_vector* x,gsl_vector* fx){
	char name[80]="out.Psitest.txt";

	double E=gsl_vector_get(x,0);
	int n=2; 
	double a=1.0/1000;
	double yA[n];
	yA[0]=a-a*a;
	yA[1]=1-2*a;
	double* ya=yA;
	double b=8;
	double h=0.01;
	double acc=0.05;
	double eps=0.05;
	gsl_vector* res=gsl_vector_alloc(n);
	driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
	gsl_vector_set(fx,0,gsl_vector_get(res,0));
	gsl_vector_free(res);
}
void func8print(gsl_vector* x,gsl_vector* fx){

        char name[80]="out.Psi8.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=8;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        gsl_vector_set(fx,0,gsl_vector_get(res,0));
        gsl_vector_free(res);
}



printf("\n\n now onto the groundstate of hydrogen:\n");
dim=1;

gsl_vector* E=gsl_vector_alloc(dim);
gsl_vector* Psib=gsl_vector_alloc(dim);

for (int i=0;i<dim;i++){
        gsl_vector_set(E,i,-0.6);
}




printf("Start point:\n");
vector_print(E);
printf("gives func value:\n");
func8newton(E,Psib);
vector_print(Psib);

eps=0.01;
printf("running algoritm \n");

Newton(func8newton,E,eps);
printf("found root::\n");
vector_print(E);
printf("The analytical result is -0.5\n");
func8print(E,Psib);
printf("function value of root:\n");
vector_print(Psib);
printf("used tolerance: %10f\n",eps);

printf("Testing for different endvalues for the =0 boundary condition, above the function was set to zero at x=8\n");
printf("same initial coniditions, tolerance etc. is used:\n");
printf("f(x)=0 for x=5\n");
void func5print(gsl_vector* x,gsl_vector* fx){

        char name[80]="out.Psi5.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=5;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        gsl_vector_set(fx,0,gsl_vector_get(res,0));
        gsl_vector_free(res);
}

void func5newton(gsl_vector* x,gsl_vector* fx){
        char name[80]="out.Psitest.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=5;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        gsl_vector_set(fx,0,gsl_vector_get(res,0));
        gsl_vector_free(res);
}

for (int i=0;i<dim;i++){
        gsl_vector_set(E,i,-0.6);
}



printf("running algoritm \n");

Newton(func5newton,E,eps);
printf("found root::\n");
vector_print(E);
printf("The analytical result is -0.5\n");
func5print(E,Psib);
printf("function value of root:\n");
vector_print(Psib);


void func65print(gsl_vector* x,gsl_vector* fx){

        char name[80]="out.Psi65.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=6.5;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        gsl_vector_set(fx,0,gsl_vector_get(res,0));
        gsl_vector_free(res);
}

void func65newton(gsl_vector* x,gsl_vector* fx){
        char name[80]="out.Psitest.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=6.5;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        gsl_vector_set(fx,0,gsl_vector_get(res,0));
        gsl_vector_free(res);
}

for (int i=0;i<dim;i++){
        gsl_vector_set(E,i,-0.6);
}




printf("running algoritm f(x)=0 for x=6.5\n");

Newton(func65newton,E,eps);
printf("found root::\n");
vector_print(E);
printf("The analytical result is -0.5\n");
func65print(E,Psib);
printf("function value of root:\n");
vector_print(Psib);

void func65print2(gsl_vector* x,gsl_vector* fx){

        char name[80]="out.Psi65_2.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=6.5;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
	double bound_value=b*exp(-b);
        gsl_vector_set(fx,0,gsl_vector_get(res,0)-bound_value);
        gsl_vector_free(res);
}

void func65newton2(gsl_vector* x,gsl_vector* fx){
        char name[80]="out.Psitest.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=6.5;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
	double bound_value=b*exp(-b);
        gsl_vector_set(fx,0,gsl_vector_get(res,0)-bound_value);
        gsl_vector_free(res);
}

void func5print2(gsl_vector* x,gsl_vector* fx){

        char name[80]="out.Psi5_2.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=5;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        double bound_value=b*exp(-b);
        gsl_vector_set(fx,0,gsl_vector_get(res,0)-bound_value);
        gsl_vector_free(res);
}

void func5newton2(gsl_vector* x,gsl_vector* fx){
        char name[80]="out.Psitest.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=5;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        double bound_value=b*exp(-b);
        gsl_vector_set(fx,0,gsl_vector_get(res,0)-bound_value);
        gsl_vector_free(res);
}

void func8print2(gsl_vector* x,gsl_vector* fx){

        char name[80]="out.Psi8_2.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=8;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        double bound_value=b*exp(-b);
        gsl_vector_set(fx,0,gsl_vector_get(res,0)-bound_value);
        gsl_vector_free(res);
}

void func8newton2(gsl_vector* x,gsl_vector* fx){
        char name[80]="out.Psitest.txt";

        double E=gsl_vector_get(x,0);
        int n=2;
        double a=1.0/1000;
        double yA[n];
        yA[0]=a-a*a;
        yA[1]=1-2*a;
        double* ya=yA;
        double b=8;
        double h=0.01;
        double acc=0.05;
        double eps=0.05;
        gsl_vector* res=gsl_vector_alloc(n);
        driverg(&g,n,a,ya,b,h,acc,eps,name,res,E);
        double bound_value=b*exp(-b);
        gsl_vector_set(fx,0,gsl_vector_get(res,0)-bound_value);
        gsl_vector_free(res);
}

printf("Now with the better boundary condition:\n");

for (int i=0;i<dim;i++){
        gsl_vector_set(E,i,-0.6);
}




printf("running algoritm f(x)=x*exp(-x) for x=5\n");

Newton(func5newton2,E,eps);
printf("found root::\n");
vector_print(E);
printf("The analytical result is -0.5\n");
func5print2(E,Psib);
printf("diference between function value and x*exp(-x) at integration end:\n");
vector_print(Psib);

for (int i=0;i<dim;i++){
        gsl_vector_set(E,i,-0.6);
}




printf("running algoritm f(x)=x*exp(-x) for x=6.5\n");

Newton(func65newton2,E,eps);
printf("found root::\n");
vector_print(E);
printf("The analytical result is -0.5\n");
func65print2(E,Psib);
printf("diference between function value and x*exp(-x) at integration end:\n");
vector_print(Psib);

printf("running algoritm f(x)=x*exp(-x) for x=8\n");

Newton(func8newton2,E,eps);
printf("found root::\n");
vector_print(E);
printf("The analytical result is -0.5\n");
func8print2(E,Psib);
printf("diference between function value and x*exp(-x) at integration end:\n");
vector_print(Psib);






gsl_vector_free(E);
gsl_vector_free(Psib);
gsl_vector_free(fX);
gsl_vector_free(X);
return 0;
}

