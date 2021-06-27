#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void plainmc(int dim, double f(int dim, double* x), double* a, double* b, int N, double* res, double* err){


double V=1;
for (int i=0;i<dim;i++){
	V*=b[i]-a[i];
}

double sum=0;
double sum2=0;
double x[dim];
unsigned int rand_seed=3;
double fx;

for (int i=0;i<N;i++){
	for (int j=0;j<dim;j++){
		x[j]=a[j]+rand_r(&rand_seed)*1.0/RAND_MAX*(b[j]-a[j]);
	}
	
	fx=f(dim,x);
	sum+=fx;
	sum2+=fx*fx;
}
double avg=sum/N;
*res=avg*V;
*err=sqrt(sum2/N-avg*avg)*V/sqrt(N);

}

double corput (int n, int base){
double q=0;
double bk=1.0/base;
while(n>0){
	q+=(n % base)*bk;
	n/=base;
	bk/=base;
}
return q;
}



void quasimc(int dim, double f(int dim, double* x), double* a, double* b, int N, double* res, double* err){


double V=1;
for (int i=0;i<dim;i++){
        V*=b[i]-a[i];
}

double Nhalf=round(N/2);

//double up on sum variables used for error estimates
double sum1=0;
double sum2=0;
double x[dim];

double fx;

//list of primes, sets a maximum dimension on the integral of 29 dimensions. more primes should be added if more dimensions needed. 
int bases[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271};

for (int i=0;i<Nhalf;i++){
        for (int j=0;j<dim;j++){
                x[j]=a[j]+corput(i,bases[j])*(b[j]-a[j]);
        }

        fx=f(dim,x);
        sum1+=fx;

	for (int j=0;j<dim;j++){
                x[j]=a[j]+corput(i,bases[j+dim])*(b[j]-a[j]);
        }

	fx=f(dim,x);
	sum2+=fx;
}
double avg=(sum2+sum1)/(2*Nhalf); //the rounding may cause 2*Nhalf=/=N,however for reasonably large N wont matter much.
*res=avg*V;
*err=fabs(sum1-sum2)/(2*Nhalf)*V;

}






double stratamc(int dim, double f(int dim, double* x), double* a, double* b,double acc, double eps, int n_reuse, double mean_reuse){

int N=16*dim;
double V=1;
for (int i=0;i<dim;i++){
	V*=b[i]-a[i];
}

int n_l[dim];
int n_r[dim];

double x[dim];
double mean_l[dim];
double mean_r[dim];
double mean=0;

unsigned int rand_seed=3;


for (int i=0;i<N;i++){
	for(int j=0;j<dim;j++){
		x[j]=a[j]+rand_r(&rand_seed)*1.0/RAND_MAX*(b[j]-a[j]);
	}

	double fx=f(dim,x);
	for (int j=0;j<dim;j++){
		if(x[j]>(a[j]+b[j])/2){
			n_r[j]++;
			mean_r[j]+=fx;
		} else {
			n_l[j]++;
			mean_l[j]+=fx;
		}
	}

	mean+=fx;
}
mean/=N;

for(int i=0;i<dim;i++){
	mean_l[i]/=n_l[i];
	mean_r[i]/=n_r[i];
}

int idiv=0;
double maxvar=0;

for(int i =0;i<dim;i++){
	double var=fabs(mean_r[i]-mean_l[i]);
	if(var>maxvar){
		maxvar=var;
		idiv=i;
	}
}


double integ=(mean*N+mean_reuse*n_reuse)/(N+n_reuse)*V;
double err=fabs(mean_reuse-mean)*V;
double tol=acc+fabs(integ)*eps;
if (err<tol){
	return integ;
}

double a2[dim];
double b2[dim];

for (int i=0;i<dim;i++){
	a2[i]=a[i];
	b2[i]=b[i];
}

a2[idiv]=(a[idiv]+b[idiv])/2;
b2[idiv]=(a[idiv]+b[idiv])/2;
double integ_l=stratamc(dim,f,a,b2,acc/sqrt(2.0),eps,n_l[idiv],mean_l[idiv]);
double integ_r=stratamc(dim,f,a2,b,acc/sqrt(2.0),eps,n_r[idiv],mean_r[idiv]);
return integ_l+integ_r;

}



