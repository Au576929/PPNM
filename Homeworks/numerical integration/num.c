#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>

#include "num_funcs.h"

int main () {

double xmin=0;
double xmax=1;

double acc=0.0001;
double eps=0.0001;

int nrec=0;
int calls=0;

double Err=0;



double func(double x){
calls++;
return sqrt(x);
}

double I=Integ(func,xmin,xmax,acc,eps,nrec,&Err);
printf("integral of sqrt(x) from 0 to 1\n");
printf("Nummerical result=%f \n Analytical result 2/3 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);
printf("The error estimate is somewhat pessimistic by the nature of the calculation method, this goes for the entire work.\n");

Err=0;
calls=0;
nrec=0;
acc=1e-9;
eps=1e-9;

double func2(double x){
calls++;
return 4*sqrt(1-x*x);
}

I=Integ(func2,xmin,xmax,acc,eps,nrec,&Err);
printf("\n\n\nintegral of 4*sqrt(1-x²) from 0 to 1\n");
printf("Nummerical result=%f \n Analytical result pi \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);



double func3(double x){
calls++;
return 1.0/sqrt(x);
}


xmin=0;
xmax=1;
nrec=0;
Err=0;
calls=0;
acc=1e-4;
eps=1e-4;

I=Integ(func3,xmin,xmax,acc,eps,nrec,&Err);
printf("\n\n\nintegral of 1/sqrt(x) from 0 to 1:\n Without Clenshaw-Curtis trnaformation:\n");
printf("Nummerical result=%f \n Analytical result 2 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);


nrec=0;
Err=0;
calls=0;

double func4(double x){
calls++;
return 1.0/sqrt(x/2.0+1.0/2)*0.5;
}




I=IntegCC(func4,acc,eps,nrec,&Err);

printf("\n\n With the Clenshaw-Curtis transformation:\n");
printf("Nummerical result=%f \n Analytical result 2 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);


double func5(double x){
calls++;
return log(x)/sqrt(x);
}


xmin=0;
xmax=1;
nrec=0;
Err=0;
calls=0;




I=Integ(func5,xmin,xmax,acc,eps,nrec,&Err);
printf("\n\n\nintegral of log(x)/sqrt(x) from 0 to 1:\n Without Clenshaw-Curtis trnaformation:\n");
printf("Nummerical result=%f \n Analytical result -4 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);


nrec=0;
Err=0;
calls=0;


double func6(double x){
calls++;
return log(x/2+1.0/2)/sqrt(x/2+1.0/2)*0.5;
}




I=IntegCC(func6,acc,eps,nrec,&Err);

printf("\n\n With the Clenshaw-Curtis transformation:\n");
printf("Nummerical result=%f \n Analytical result -4 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);



nrec=0;
Err=0;
calls=0;
acc=1e-9;
eps=acc;



double func7(double x){
calls++;
return 4.0*sqrt(1.0-(x/2+1.0/2)*(x/2+1.0/2))*0.5;
}

printf("The final integral of part B (without transform) is the same as te final one in part A...");


I=IntegCC(func7,acc,eps,nrec,&Err);

printf("\n\n integral of 4*sqrt(1-x*x) from 0 to 1 with the Clenshaw-Curtis transformation:\n");
printf("Nummerical result=%g \n Analytical result 3.14159265358979324 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);


double call_lim1=1e8;
size_t limit=1e7;
gsl_integration_workspace* spacefunc7 = gsl_integration_workspace_alloc(call_lim1);
calls=0;
double gslFunc7(double x,void* p){
	calls++;
	return 4.0*sqrt(1.0-x*x);
}


gsl_function gslfunc7;
gslfunc7.function=gslFunc7;
double res=0;
double abserr;
printf("with gsl's gsl_integration_qags method:\n");
gsl_integration_qags(&gslfunc7,0,1,acc,eps,limit, spacefunc7,&res,&abserr);
printf("with gsl method res= %g\n with error %g\n function called %i times\n",res,abserr,calls);


acc=0.0001;
eps=acc;
xmin=0;
xmax=INFINITY;
nrec=0;
Err=0;
calls=0;


double func8(double x){
calls++;
return exp(-1.0*x);
}

I=Integ(func8,xmin,xmax,acc,eps,nrec,&Err);

printf("\n\n integral of exp(-x) from 0 to infty:\n");
printf("Nummerical result=%f \n Analytical result 1 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);

double func9(double x){
calls++;
return exp(-x*x);
}


acc=0.00001;
eps=0.00001;

xmin=-INFINITY;
xmax=INFINITY;
nrec=0;
Err=0;
calls=0;





I=Integ(func9,xmin,xmax,acc,eps,nrec,&Err);

printf("\n\n integral of exp(-x²) from -infty to infty:\n");
printf("Nummerical result=%f \n Analytical result sqrt(pi)=1.7724538509 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);



xmin=-INFINITY;
xmax=0;
nrec=0;
Err=0;
calls=0;

I=Integ(func9,xmin,xmax,acc,eps,nrec,&Err);

printf("\n\n integral of exp(-x²) from -infty to 0:\n");
printf("Nummerical result=%f \n Analytical result sqrt(pi)/2=0.886226925 \nFunction called %i times\n",I,calls);
printf("Integration error %f\n",Err);


printf("\n\nFrom here on gsl routines are used!\n\n");

acc=0.001;
eps=0.001;

xmin=-5;
xmax=5;
nrec=0;
Err=0;
calls=0;
size_t call_lim=70000;
//size_t callgsl;

double Res=0;

gsl_integration_workspace* test = gsl_integration_workspace_alloc(call_lim);

double test2 (double x,void* p){
	calls++;
	return exp(-x*x);
}


gsl_function func10;
func10.function=test2;


gsl_integration_qags(&func10,xmin,xmax,acc,eps,call_lim,test,&Res,&Err);

printf("\n\n integral of exp(-x²) from -5 to 5:\n");
printf("Nummerical result=%f \n Wolfram alpha result 1.77245\nFunction called %i times\n",Res,calls);
printf("Integration error %f\n",Err);

calls=0;

gsl_integration_qagil(&func10,xmax,acc,eps,call_lim,test,&Res,&Err);


printf("\n\n integral of exp(-x²) from -Infinity to 5:\n");
printf("Nummerical result=%f \n Wolfram alpha result 1.77245\nFunction called %i times\n",Res,calls);
printf("Integration error %f\n",Err);

calls=0;
gsl_integration_qagiu(&func10,xmin,acc,eps,call_lim,test,&Res,&Err);


printf("\n\n integral of exp(-x²) from -5 to infinity:\n");
printf("Nummerical result=%f \n Wolfram alpha result 1.77245\nFunction called %i times\n",Res,calls);
printf("Integration error %f\n",Err);
calls=0;
gsl_integration_qagi(&func10,acc,eps,call_lim,test,&Res,&Err);


printf("\n\n integral of exp(-x²) from -infinity to infinity:\n");
printf("Nummerical result=%f \n Wolfram alpha result 1.77245\nFunction called %i times\n",Res,calls);
printf("Integration error %f\n",Err);




gsl_integration_workspace_free(test);
return 0;
}
