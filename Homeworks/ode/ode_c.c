#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ode_funcs.h"


int main(){


int n=12;
double a;
double yA[n];

a=0.0;

yA[0]=0.97000436;	 //x-coordinate of particle 1
yA[1]=-0.24308753;	 //y-coordinate of particle 1
yA[4]=-yA[0];		 //x-coordinate of particle 2 
yA[5]=-yA[1];            //y-coordinate of particle 2 
yA[8]=0;                 //x-coordinate of particle 3 
yA[9]=0;                 //y-coordinate of particle 3 

yA[10]=-0.93240737;        //x-velocity of particle 3
yA[11]=-0.86473146;         //y-velocity of particle 3
yA[6]=-yA[10]/2;            //x-velocity of particle 2 
yA[7]=-yA[11]/2;            //y-velocity of particle 2 
yA[2]=-yA[10]/2;            //x-velocity of particle 1 
yA[3]=-yA[11]/2;             //y-velocity of particle 1 

double* ya=yA;
double b=10; 

double h=0.001;
double acc=0.01;
double eps=0.005;

char name[80]="out.gravity.txt";
driverg(&g,n,a,ya,b,h,acc,eps,name);


return  0;
}
