#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "mini_funcs.h"


int main(){

void task_doer(double func(gsl_vector* x),gsl_vector* x,double eps,char* func_description){
	int dim=x->size;
	char string1[80];
	sprintf(string1, "%s%s\n","Function to be minimized: f(x)=",func_description);
	puts(string1);
	printf("starting x=%10f\n",gsl_vector_get(x,0));
        for (int i =1;i<dim;i++){
                printf("           %10f\n",gsl_vector_get(x,i));
        }
	printf("Gives: f(x)=%10f\n",func(x));

	
	int steps=minimize(func,x,eps);
	printf("minimization done:\n x=%10f\n",gsl_vector_get(x,0));
	for (int i =1;i<dim;i++){
		printf("   %10f\n",gsl_vector_get(x,i));
	}
	printf("Accomplished in %i steps\n",steps);
	printf("gives: f(x)=%10f",func(x));
	//rintf("done in %i steps\n",*steps);
	printf("\n\n\n");
}


int dim1=1;
gsl_vector* x1=gsl_vector_alloc(dim1);
gsl_vector_set(x1,0,3);
double eps1=0.5;

double func1(gsl_vector* x){
	double X=gsl_vector_get(x,0);
	return X*X;
}


char func_description1[50]="x²";

task_doer(func1,x1,eps1,func_description1);




int dim2=2;
gsl_vector* x2=gsl_vector_alloc(dim2);
gsl_vector_set(x2,0,5);
gsl_vector_set(x2,1,5);
double eps2=0.005;

double func2(gsl_vector* x){
	double X=gsl_vector_get(x,0);
	double Y=gsl_vector_get(x,1);
	return (1-X)*(1-X)+100*(Y-X*X)*(Y-X*X);
}

char func_description2[50]="(1-x)²+100(y-x²)²";

task_doer(func2,x2,eps2,func_description2);


double func3(gsl_vector* x){
	double X=gsl_vector_get(x,0);
	double Y=gsl_vector_get(x,1);
	return (X*X+Y-11.0)*(X*X+Y-11.0)+(X+Y*Y-7.0)*(X+Y*Y-7.0);
}

int dim3=2;
gsl_vector* x3=gsl_vector_alloc(dim3);
gsl_vector_set(x3,0,3.5);
gsl_vector_set(x3,1,2.9);
double eps3=0.01;

char func_description3[50]="(x²+y-11)²+(x+y²-7)²";

task_doer(func3,x3,eps3,func_description3);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double breit_wigner(double E,gsl_vector* X){
                double A=gsl_vector_get(X,0);
                double m=fabs(gsl_vector_get(X,1));
                double Gamma=fabs(gsl_vector_get(X,2));

                return A/((E-m)*(E-m)+Gamma*Gamma/4.0);
        }


double func4(gsl_vector* X){
	int ndata=20;
	double Es[]={101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 151, 153, 155, 157, 159};
	double c_sec[]={-0.25, -0.3, -0.15, -1.71, 0.81, 0.65, -0.91, 0.91, 0.96, -2.52, -1.01, 2.01, 4.83, 4.58, 1.26, 0.45, 0.15, -0.9,1 -0.81, -1.41, 1.36, 0.5, -0.45, 1.61, -2.21, -1.86, 1.76, -0.5};
	double errs[]={2, 2, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 0.9, 0.9, 0.9};

	double s=0;
	double bw;
	
	for (int i=0;i<ndata;i++){
	bw=breit_wigner(Es[i],X);
	double c=c_sec[i];
	double err=errs[i];
	s+=(bw-c)*(bw-c)/(err*err);
	}
	return s;
}
	
	int dim4=3;
	gsl_vector* x4=gsl_vector_alloc(dim4);
	gsl_vector_set(x4,0,1);
	gsl_vector_set(x4,1,130);
	gsl_vector_set(x4,2,0.003);
	double eps4=1e-10;

	char func_description4[50]="Deviation function for ex b";
	task_doer(func4,x4,eps4,func_description4);
	printf("resultant mass: %f GeV\n According to wikipedia: m=125.1 GeV\n resultant width: %f\n",fabs(gsl_vector_get(x4,1)),fabs(gsl_vector_get(x4,2)));
	FILE* stream=fopen("out.bwfit.txt","w");
	double bw;
	for(int i=0;i<60*10;i++){
	bw=breit_wigner(100+i*1.0/10,x4);
	fprintf(stream,"%10f  %10f\n",100+i*1.0/10,bw);
	}


printf("\n\n\nSIMPLEX METHOD FROM HERE ON OUT:\n");
printf("Function used: (1-x)²+100(y-x²)²\n");


int dim=2;
gsl_matrix* simplex=gsl_matrix_calloc(dim,dim+1);
gsl_vector* f_values=gsl_vector_calloc(dim+1);
gsl_vector* holder=gsl_vector_calloc(dim);
unsigned int rand_seed=3;

for (int i=0;i<dim+1;i++){
	for (int j=0;j<dim;j++){
		gsl_matrix_set(simplex,j,i,rand_r(&rand_seed)*5.0/RAND_MAX);
	}
	vec_extract(simplex,holder,i);
	gsl_vector_set(f_values,i,func2(holder));
}

printf("Simplex before minimization:\n");
for (int i=0;i<dim;i++){
	for (int j=0;j<dim+1;j++){
		printf("%10f  ",gsl_matrix_get(simplex,i,j));
	}
	printf("\n");
}
printf("\n\n");

printf("Function values before minimization:\n");

for (int i=0;i<dim+1;i++){
	printf("%10f  ",gsl_vector_get(f_values,i));
}
printf("\n\n");

double size_goal=0.01;
int ring_count=simplex_ringdown(func2,simplex,size_goal);

printf("number of iterations used for minimization=%i\n",ring_count);

printf("Simplex after minimization:\n");

for (int i=0;i<dim;i++){
        for (int j=0;j<dim+1;j++){
                printf("%10f  ",gsl_matrix_get(simplex,i,j));
        }
        printf("\n");
}
printf("\n\n");

printf("Function values after minimization:\n");

for (int i=0;i<dim+1;i++){
	vec_extract(simplex,holder,i);
        printf("%10f  ",func2(holder));
}
printf("\n\n");






gsl_matrix_free(simplex);
gsl_vector_free(f_values);
gsl_vector_free(holder);
gsl_vector_free(x4);
gsl_vector_free(x3);
gsl_vector_free(x2);
gsl_vector_free(x1);

return 0;
}

