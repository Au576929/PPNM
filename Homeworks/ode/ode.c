#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ode_funcs.h"


int main(){


int n=3;
double a;
double yA[n];
int N=5.5e6;

a=0.0;
yA[1]=661; //Infectious 
yA[2]=5e5; //Removed
yA[0]=N-yA[2]-yA[1]; //susceptible

double* ya=yA;
double b=100; //days simulated

double h=0.333;
double acc=0.0005;
double eps=0.00005;
double Tr=10;
char pre_name[20]="out.Tc";
char post_name[20]=".txt";
char name[80];
double Tc=5;
for (int i=1;i<=5;i++){
sprintf(name,"%s%i%s",pre_name,i,post_name);
driver(&f,n,a,ya,b,h,acc,eps,Tc,Tr,N,name);
Tc/=2;
}


return  0;
}
