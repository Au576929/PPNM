#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_vector.h>

void rkstep12g ( //This stepper will function in with pretty much any function in place of g, used for part c
	void g(int n, double x, double* y, double* dydx,double E),
	int n,
	double x,
	double* yx,
	double h,
	double* yh,
	double* dy,
	double E
){
double k0[n];
double k12[n];

double yt[n];

g(n,x,yx,k0,E);
for (int i=0; i<n; i++){
	yt[i]=yx[i]+k0[i]*h/2;
}
g(n,x+h/2,yt,k12,E);
for (int i=0; i<n; i++){
	yh[i]=yx[i]+k12[i]*h;
}
for (int i=0;i<n;i++){
	dy[i]=(k0[i]-k12[i])*h/2;
}

}


void rkstep12Tc ( //this stepper was created in order to pass the variables Tc,Tr and N to the function f
        void f(int n, double x, double* y, double* dydx,double Tc, double Tr,int N),
        int n,
        double x,
        double* yx,
        double h,
        double* yh,
        double* dy,
	double Tc,
	double Tr,
	int N
){
double k0[n];
double k12[n];

double yt[n];

f(n,x,yx,k0,Tc,Tr,N);
for (int i=0; i<n; i++){
        yt[i]=yx[i]+k0[i]*h/2;
}
f(n,x+h/2,yt,k12,Tc,Tr,N);
for (int i=0; i<n; i++){
        yh[i]=yx[i]+k12[i]*h;
}
for (int i=0;i<n;i++){
        dy[i]=(k0[i]-k12[i])*h/2;
}
}


void f(int n, double x, double* y, double* dydx,double Tc,double Tr,int N){
//double N=5.5e6;
//double Tc=4.5;
//double Tr=9;
dydx[0]=-1.0*(y[1]*y[0])/(N*Tc);
assert(dydx[0]<0);
dydx[1]=(y[1]*y[0])/(N*Tc)-1.0*y[1]/Tr;
dydx[2]=(y[1])/Tr;
}


void g(int n, double x, double* y, double* dydx){

double r12=(y[0]-y[4])*(y[0]-y[4])+(y[1]-y[5])*(y[1]-y[5]); //distance from particle 1 to particle 2 SQUARED!
double r13=(y[0]-y[8])*(y[0]-y[8])+(y[1]-y[9])*(y[1]-y[9]); //---------//----------- 1 -----//---- 3 ---//---!
double r23=(y[4]-y[8])*(y[4]-y[8])+(y[5]-y[9])*(y[5]-y[9]); // you know

//acceleration on particle 1 along x
dydx[2]=-(y[0]-y[4])/pow(r12,3.0/2)-(y[0]-y[8])/pow(r13,3.0/2);
//acceleration on particle 1 along y
dydx[3]=-(y[1]-y[5])/pow(r12,3.0/2)-(y[1]-y[9])/pow(r13,3.0/2);

//acceleration on particle 2 along x
dydx[6]=-(y[4]-y[0])/pow(r12,3.0/2)-(y[4]-y[8])/pow(r23,3.0/2);
//acceleration on particle 2 along y
dydx[7]=-(y[5]-y[1])/pow(r12,3.0/2)-(y[5]-y[9])/pow(r23,3.0/2);

//acceleration on particle 3 along x
dydx[10]=-(y[8]-y[0])/pow(r13,3.0/2)-(y[8]-y[4])/pow(r23,3.0/2);
//acceleration on particle 3 along y
dydx[11]=-(y[9]-y[1])/pow(r13,3.0/2)-(y[9]-y[5])/pow(r23,3.0/2);


//velocity of particle 1 along x
dydx[0]=y[2];
//velocity of particle 1 along y
dydx[1]=y[3];


//velocity of particle 2 along x
dydx[4]=y[6];
//velocity of particle 2 along y
dydx[5]=y[7];

//velocity of particle 3 along x
dydx[8]=y[10];
//velocity of particle 3 along y
dydx[9]=y[11];


}




/*void f(int n, double x, double* y, double* dydx){
dydx[0]=y[1];
dydx[1]=-y[0];
}*/



void driver ( //this driver passes Tc, Tr, and N so therefore wont be used for part c
	void f(int n, double x, double* y, double* dydx,double Tc,double Tr, int N),
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
){
double x=a;
double* y=ya;

double s, normy, yh[n], dy[n], tol, err;
FILE* list_stream=fopen(name,"w");

while (x<b){ /*if b<a in your is of interest, it is simple to change x->-x in diff lign and then use that instead  */

	if(x+h>b){
		h=b-x;
	}

	rkstep12Tc(f,n,x,y,h,yh,dy,Tc,Tr,N);
	s=0;
	for (int i=0;i<n;i++){
		s+=dy[i]*dy[i];
	}
	err=sqrt(s);
	s=0;
	for (int i=0;i<n;i++){
        	s+=yh[i]*yh[i];
	}
	normy=sqrt(s);

	tol=(normy*eps+acc)*sqrt(h/(b-a));
	if (err<tol){
		x+=h;
		y=yh;
		fprintf(list_stream,"%10g ",x);
		for (int i=0;i<n;i++){
			fprintf(list_stream,"%10g ",yh[i]);
		}
		fprintf(list_stream,"\n");
	} 
        h*=pow(tol/err,0.25)*0.95;
 
}
}




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
){
double x=a;
double* y=ya;

double s, normy, yh[n], dy[n], tol, err;
FILE* list_stream=fopen(name,"w");

while (x<b){ /*if b<a in your is of interest, it is simple to change x->-x in diff lign and then use that instead  */

        if(x+h>b){
                h=b-x;
        }

        rkstep12g(g,n,x,y,h,yh,dy,E);
        s=0;
        for (int i=0;i<n;i++){
                s+=dy[i]*dy[i];
        }
        err=sqrt(s);
        s=0;
        for (int i=0;i<n;i++){
                s+=yh[i]*yh[i];
        }
        normy=sqrt(s);

        tol=(normy*eps+acc)*sqrt(h/(b-a));
        if (err<tol){
                x+=h;
                y=yh;
		fprintf(list_stream,"%10f ",x);
		for (int i=0; i<n;i++){
			fprintf(list_stream,"%10f ",y[i]);
		}
		fprintf(list_stream,"\n");
	} 
        h*=pow(tol/err,0.25)*0.95;
 
}

for (int i=0;i<n;i++){
	gsl_vector_set(res,i,yh[i]);
}


fprintf(list_stream,"\n \n \n \n \n"); //creating new index in out file
}

