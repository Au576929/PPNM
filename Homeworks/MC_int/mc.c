#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mc_funcs.h"

int main(){

//(int dim, double f(int dim, double* x), double* a, double* b, int N, double res, double eps)


double func1(int dim, double* p){
	double x=p[0];
	return exp(-x);

}
double p[]={-1};
int dim=sizeof(p)/sizeof(p[0]);
printf("Main: checking func1: %f\n",func1(dim,p));

printf("integral of exp(-x) from 0 to 5, with 10000 points:");
double a[]={0.0};
double b[]={5.0};
int N=10000;
double res=0;
double err=0;
printf("Main: plainmc not called. res=%f\n",res);

plainmc(dim,func1,a,b,N,&res,&err);
printf("Main: plainmc called. \n res=%f\nerr=%f\nactual result: 0.99326\n",res,err);

printf("integral of cos(x)*cos(y) from 0 to 3.14159265 in both x and y, with 25000 points:");


double func2(int dim, double* p){
        double x=p[0];
	double y=p[1];
        return cos(x)*cos(y);

}


double a2[]={0.0,0.0};
double b2[]={3.14159265,3.14159265};
dim=sizeof(a2)/sizeof(a2[0]);

N=25000;
plainmc(dim,func2,a2,b2,N,&res,&err);
printf("Main: plainmc called. \n res=%f\nerr=%f\nactual result: 0\n",res,err);


printf("integral of 1/(1-cos(x)*cos(y)*cos(z)) from 0 to 3.14159265 in both x, y and z, with 100000 points:");


double func3(int dim, double* p){
        double x=p[0];
        double y=p[1];
	double z=p[2];
	double pi=3.14159265;
        return 1.0/(1-cos(x)*cos(y)*cos(z))*1.0/(pi*pi*pi);

}


double a3[]={0.0,0.0,0.0};
double b3[]={3.14159265,3.14159265, 3.14159265};
dim=sizeof(a3)/sizeof(a3[0]);

N=100000;
plainmc(dim,func3,a3,b3,N,&res,&err);
printf("Main: plainmc called. \n res=%f\nerr=%f\nactual result: 1.39320392968\n",res,err);

printf("Corput algorithm implemented.\n");
int base=10;
printf("test of corput algorithm:\n base is: %i\n",base);
for (int i=1;i<11;i++){
	printf("n=%i  corput returns: %f\n",i,corput(i,base));
}


dim=sizeof(a2)/sizeof(a2[0]);

printf("Testing quasi integration on integral of cos(x)*cos(y) from 0 to 3.14159265 in both x and y, with 25000 points:");
N=25000;
quasimc(dim,func2,a2,b2,N,&res,&err);
printf("Quasimc called. \n res=%f\nerr=%f\nactual result: 0\n",res,err);

FILE* errStream=fopen("out.err.txt","w");
for (N=10;N<1e5;N+=1000){
	quasimc(dim,func2,a2,b2,N,&res,&err);
	fprintf(errStream,"%10i  %10f",N,err);
	plainmc(dim,func2,a2,b2,N,&res,&err);
	fprintf(errStream,"%10f\n",err);
}

printf("Testing stratified mc algorithm:\n");

double acc=0.005;
double eps=0;
int n_reuse=0;
double mean_reuse=0;
dim=sizeof(a2)/sizeof(a2[0]);

res=stratamc(dim,func2,a2,b2,acc,eps,n_reuse,mean_reuse);

printf("Quasimc called. \n res=%f\nactual result: 0\naccuracy %f\n eps %f\n",res,acc,eps);





printf("\n");
return 0;
}
