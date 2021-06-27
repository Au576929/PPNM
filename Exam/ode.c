#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "ode_funcs.h"


int main(){

int n=2;
double complex x0[n];
x0[0]=0+0*I;
x0[1]=1+0*I;


void harm_osc (int dim, double complex t, double complex* x, double complex* dx){

	dx[1]=-x[0];
	dx[0]=x[1];

}

double complex a=0;
double complex b=20;

double acc=0.00005;
double eps=0.00005;
double dt=0.1;

char name_osc1[80]="out.osc1.txt";

driver(&harm_osc,n,a,x0,b,dt,acc,eps,name_osc1);

a=0;
b=5*I;

x0[0]=1;
x0[1]=0;

char name_osc2[80]="out.osc2.txt";

driver(&harm_osc,n,a,x0,b,dt,acc,eps,name_osc2);

a=0;
b=5+5*I;

x0[0]=0;
x0[1]=1;

char name_osc3[80]="out.osc3.txt";

driver(&harm_osc,n,a,x0,b,dt,acc,eps,name_osc3);

return 0;
}
