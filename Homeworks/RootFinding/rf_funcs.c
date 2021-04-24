#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

#include "lin_funcs.h"


void Newton(void f(gsl_vector* X, gsl_vector* fX), gsl_vector* X, double eps){
double dx=sqrt(3e-16); //close enough to machine epsilon 
//double dx=0.000001;

int dim = X->size;
double lambda=1;
double lambda_min=dx;
int counter=0;

gsl_vector* fdX=gsl_vector_alloc(dim);
gsl_vector* fX=gsl_vector_alloc(dim);
gsl_vector* dX=gsl_vector_alloc(dim);
gsl_matrix* J=gsl_matrix_alloc(dim,dim);
gsl_matrix* R=gsl_matrix_alloc(dim,dim);

for (int i=0;i<dim;i++){
	gsl_vector_set(dX,i,0);
}


do{
	counter++;
	assert(counter<10000);
	f(X,fX);
	gsl_blas_dscal(0,dX);
	gsl_vector_set(dX,0,dx);
	gsl_blas_daxpy(1.0,dX,X);
	f(X,fdX);
	gsl_blas_daxpy(-1.0,dX,X);
	for (int i=0; i<dim;i++){
		gsl_matrix_set(J,i,0,(gsl_vector_get(fdX,i)-gsl_vector_get(fX,i))/dx);
	}

	for (int i=1;i<dim;i++){
		gsl_vector_set(dX,i-1,0);
		gsl_vector_set(dX,i,dx);
		
		gsl_blas_daxpy(1.0,dX,X);
		f(X,fdX);
		gsl_blas_daxpy(-1.0,dX,X);
		for (int j=0; j<dim;j++){
        	        gsl_matrix_set(J,j,i,(gsl_vector_get(fdX,j)-gsl_vector_get(fX,j))*1.0/dx);
	        }
	}

	GS_decomp(J,R);
	gsl_blas_daxpy(-2.0,fX,fX);
	GS_solve(J,R,fX,dX);
	gsl_blas_daxpy(-2.0,fX,fX);
	gsl_blas_daxpy(1.0,dX,X);
	f(X,fdX);
	gsl_blas_daxpy(-1.0,dX,X);
	while(gsl_blas_dnrm2(fdX)>(1-lambda/2)*gsl_blas_dnrm2(fX) && lambda>lambda_min){
		lambda/=2;
	}

	gsl_blas_daxpy(lambda,dX,X);
	f(X,fX);

} while (gsl_blas_dnrm2(fX)>eps);

gsl_vector_free(dX);
gsl_vector_free(fX);
gsl_vector_free(fdX);
gsl_matrix_free(J);
gsl_matrix_free(R);
}

void vector_print(gsl_vector* x){
for (int i=0;i<x->size;i++){
	printf("%10g\n",gsl_vector_get(x,i));
}
printf("\n");
}


void maxtrix_print(gsl_matrix* X){
for (int i=0;i<X->size1;i++){
	for (int j=0;j<X->size2;j++){
		printf("%10g",gsl_matrix_get(X,j,i));
	}
	printf("\n");
}
}
