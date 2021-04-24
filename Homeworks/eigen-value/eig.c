#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <time.h>

#include "eig_funcs.h"


int main(){
int n=5;

gsl_matrix* A=gsl_matrix_alloc(n,n);
gsl_matrix* V=gsl_matrix_alloc(n,n);
gsl_matrix* A_copy=gsl_matrix_alloc(n,n);
gsl_matrix* Prodsi=gsl_matrix_alloc(n,n);
gsl_matrix* Prodsii=gsl_matrix_alloc(n,n);
double r=0;
unsigned int rand_seed=1;
for (int i=0; i<n;i++){
	for (int j=i;j<n;j++){
		r=1.0*rand_r(&rand_seed)/RAND_MAX;
		gsl_matrix_set(A,i,j,r);
		gsl_matrix_set(A,j,i,r);
		gsl_matrix_set(V,i,j,0);
	}
	gsl_matrix_set(V,i,i,1);
}

//A is copied, for use after diagonalisation
gsl_matrix_memcpy(A_copy,A);

//abvoe the A matrix is made symmetrically and V is set to identity matrix.
//there is a bit of double setting og entrencies, e.g. the diagonal of V, but it should not take much time.
printf("Before diag is done:\n A=\n");
matrixprint(A);

printf("V=\n");
matrixprint(V); 

jacobi_diag(A,V);
printf("diagonlisation is done:\n");






printf("D=\n");
matrixprint(A);

printf("V=\n");
matrixprint(V);

printf("V^T A V=\n");
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A_copy,V,0.0,Prodsi);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,Prodsi,0.0,Prodsii);
matrixprint(Prodsii);



printf("V D V^T=\n");
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,A,V,0.0,Prodsi);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,Prodsi,0.0,Prodsii);
matrixprint(Prodsii);

printf("V^T V=\n");
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,V,0.0,Prodsi);
matrixprint(Prodsi);

printf("BEGIN INFINITE WELL PART:\n");

int N=20;
//double s=1/N;

gsl_matrix* H=gsl_matrix_alloc(N,N);
gsl_matrix* Psi=gsl_matrix_alloc(N,N);


for (int i=0;i<N;i++){
	gsl_matrix_set(Psi,i,i,1);
}

for (int i=0;i<N-1;i++){
		gsl_matrix_set(H,i,i,2*N*N);
		gsl_matrix_set(H,i,i+1,-N*N);
		gsl_matrix_set(H,i+1,i,-N*N);
}
gsl_matrix_set(H,N-1,N-1,2*N*N);

printf("Before diag:\n H=\n");
matrixprint(H);
printf("States=\n");
matrixprint(Psi);

jacobi_diag(H,Psi);

printf("After diag:\n H=\n");
matrixprint(H);
printf("States=\n");
matrixprint(Psi);
//printf("Only the first 10 energies are printed, as the matricies are chossen to be %ix%i, and such too large to print.\n",N,N);



printf("Energies:\n");
printf("Nummerical   Exact\n");
double exact_E;
for (int i=0; i<N;i++){
exact_E=M_PI*M_PI*(i+1)*(i+1);
printf("%10g %10g\n",gsl_matrix_get(H,i,i),exact_E);
}


FILE* States_Stream=fopen("out.states.txt","w");
for (int i=0;i<3;i++){
	fprintf(States_Stream,"%10g %10g\n",0.0,0.0);
	for (int j=0;j<N;j++){
		fprintf(States_Stream,"%10g %10g\n",1.0*(j+1)/(N+1),gsl_matrix_get(Psi,j,i));
	}
	fprintf(States_Stream,"%10g %10g",1.0,0.0);
	fprintf(States_Stream,"\n\n\n\n\n");
}

printf("BEGIN PART C ON TIMING\n");

clock_t t=clock();
double time_taken;

FILE* my_times_stream=fopen("out.my.times1.txt","w");


for (int n=20;n<200;n+=20){
gsl_matrix* M=gsl_matrix_alloc(n,n);
gsl_matrix* Mv=gsl_matrix_alloc(n,n);

for (int i=0; i<n;i++){
        for (int j=i;j<n;j++){
                r=1.0*rand_r(&rand_seed)/RAND_MAX;
                gsl_matrix_set(M,i,j,r);
                gsl_matrix_set(M,j,i,r);
                gsl_matrix_set(Mv,i,j,0);
        }
        gsl_matrix_set(Mv,i,i,1);
}

t=clock();

jacobi_diag(M,Mv);

t=clock()-t;

time_taken=((double) t)/CLOCKS_PER_SEC;
fprintf(my_times_stream,"%10i  %10g\n",n,time_taken);
gsl_matrix_free(M);
gsl_matrix_free(Mv);
}

FILE* my_times2_stream=fopen("out.my.times2.txt","w");


for (int n=20;n<200;n+=20){
gsl_matrix* M=gsl_matrix_alloc(n,n);
gsl_matrix* Mv=gsl_matrix_alloc(n,n);

for (int i=0; i<n;i++){
        for (int j=i;j<n;j++){
                r=1.0*rand_r(&rand_seed)/RAND_MAX;
                gsl_matrix_set(M,i,j,r);
                gsl_matrix_set(M,j,i,r);
                gsl_matrix_set(Mv,i,j,0);
        }
        gsl_matrix_set(Mv,i,i,1);
}

t=clock();

jacobi_diag2(M,Mv);

t=clock()-t;

time_taken=((double) t)/CLOCKS_PER_SEC;
fprintf(my_times2_stream,"%10i  %10g\n",n,time_taken);
gsl_matrix_free(M);
gsl_matrix_free(Mv);
}


FILE* gsl_times_stream=fopen("out.gsl.times.txt","w");


for (int n=20;n<200;n+=20){
gsl_matrix* M=gsl_matrix_alloc(n,n);
gsl_matrix* Mv=gsl_matrix_alloc(n,n);
gsl_vector* Vm=gsl_vector_alloc(n);


for (int i=0; i<n;i++){
        for (int j=i;j<n;j++){
                r=1.0*rand_r(&rand_seed)/RAND_MAX;
                gsl_matrix_set(M,i,j,r);
                gsl_matrix_set(M,j,i,r);
                gsl_matrix_set(Mv,i,j,0);
        }
        gsl_matrix_set(Mv,i,i,1);
}

t=clock();

gsl_linalg_SV_decomp_jacobi(M,Mv,Vm);

t=clock()-t;

time_taken=((double) t)/CLOCKS_PER_SEC;
fprintf(gsl_times_stream,"%10i  %10g\n",n,time_taken);
gsl_matrix_free(M);
gsl_matrix_free(Mv);
gsl_vector_free(Vm);
}




gsl_matrix_free(Psi);
gsl_matrix_free(H);

gsl_matrix_free(Prodsii);
gsl_matrix_free(A);
gsl_matrix_free(V);
gsl_matrix_free(A_copy);
gsl_matrix_free(Prodsi);
return 0; 
}
