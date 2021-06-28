#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>

void rkstep12 (
	void f(int n, double complex x, double complex* y, double complex* dydx),
	int n,
	double complex x,
	double complex* yx,
	double complex h,
	double complex* yh,
	double complex* dy
){
double complex k0[n];
double complex k12[n];

double complex yt[n];

f(n,x,yx,k0);
for (int i=0; i<n; i++){
	yt[i]=yx[i]+k0[i]*h/2;
}
f(n,x+h/2,yt,k12);
for (int i=0; i<n; i++){
	yh[i]=yx[i]+k12[i]*h;
}
for (int i=0;i<n;i++){
	dy[i]=(k0[i]-k12[i])*h/2;
}

}

void driver ( 
        void f(int n, double complex x, double complex* y, double complex* dydx),
        int n,
        double complex a,
        double complex* ya,
        double complex b,
        double dt,
        double acc,
        double eps,
        char* name
){
// The idea is to describe the function input via a real variable t going from 0 to 1


double t=0;
double complex x=a;
assert(x==a);

double complex* y=ya;

double complex dhdt=b-a; //this is then what dt (the step) needs to be multiplied with to produce h from the real setup

double s, normy, tol, err;
double complex yh[n], dy[n];
FILE* list_stream=fopen(name,"w");

while (t<1){ 
        if(t+dt>1){
                dt=1-t;
        }

        rkstep12(f,n,x,y,dt*dhdt,yh,dy);
        s=0;
        for (int i=0;i<n;i++){
                s+=cabs(dy[i])*cabs(dy[i]);
        }
        err=sqrt(s);
	
        s=0;
        for (int i=0;i<n;i++){
                s+=cabs(yh[i])*cabs(yh[i]);
        }
        normy=sqrt(s);
        tol=(normy*eps+acc)*sqrt(dt);
//        assert(err<100*tol);
	if (err<tol){
                x+=dt*dhdt;
		t+=dt;
		y=yh;
                fprintf(list_stream,"%10g %10g ",creal(x),cimag(x));
                for (int i=0;i<n;i++){
                        fprintf(list_stream,"%10g %10g ",creal(yh[i]),cimag(yh[i]));
                }
                fprintf(list_stream,"\n");
        } 
        dt*=pow(tol/err,0.25)*0.95; 
	assert(dt>1e-13);
}
}

