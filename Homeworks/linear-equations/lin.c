#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_linalg.h>


#include "lin_funcs.h"

int main(){

int size[2];

size[0]=7;
size[1]=4;
assert(size[0]>=size[1]); //checking that variable number
//			    does not exceed the number of equations

gsl_matrix* A=gsl_matrix_alloc(size[0],size[1]);
gsl_matrix* R=gsl_matrix_alloc(size[1],size[1]);
gsl_vector* b=gsl_vector_alloc(size[1]);
//gsl_vector* x=gsl_vector_alloc(size[1]);

unsigned int rand_seed=2;

for (int i = 0; i<size[1]; i++){
	for (int j=0; j<size[0]; j++){
		gsl_matrix_set(A,j,i,1.0*rand_r(&rand_seed)/RAND_MAX);
		
	}
	gsl_vector_set(b,i,1.0*rand_r(&rand_seed)/RAND_MAX);
}

printf("QR-decomposition.\n b=\n");
for (int i =0; i<size[1]; i++){ 
printf("%10g \n",gsl_vector_get(b,i));
}
printf("A=\n");
for (int i =0; i<size[0];i++){
	for(int j=0; j<size[1];j++){
		printf("%10g   ",gsl_matrix_get(A,i,j));
	}
	printf("\n");
}

printf("QR decomposition is done\n");
GS_decomp(A,R);

printf("After decomp:\n Q=A=\n");

for (int i =0; i<size[0];i++){
	for(int j=0; j<size[1];j++){
		printf("%10g   ",gsl_matrix_get(A,i,j));
	}
	printf("\n");
}

printf("R=\n");

for (int i =0; i<size[1];i++){
	for(int j=0; j<size[1];j++){
		printf("%10g   ",gsl_matrix_get(R,i,j));
	}
	printf("\n");
}

printf("testing that Q'*Q=I :\n Q'*Q=\n ");

for (int i =0; i<size[1];i++){
	for(int j=0; j<size[1];j++){
		gsl_vector_view ai=gsl_matrix_column(A,i);
		gsl_vector_view aj=gsl_matrix_column(A,j);
		double dot;
		gsl_blas_ddot(&ai.vector,&aj.vector,&dot);
		printf("%10g   ",dot);
	}
	printf("\n");
}

printf("testing that Q*R=A=\n");

gsl_matrix* a_new=gsl_matrix_alloc(size[0],size[1]);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,R,0.0,a_new);
for (int i =0; i<size[0];i++){
	for(int j=0; j<size[1];j++){
		printf("%10g   ",gsl_matrix_get(a_new,i,j));
	}
	printf("\n");
}

gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(a_new);

printf("in order to test the solution algorithm, square matricies are created:");

gsl_matrix* As=gsl_matrix_alloc(size[1],size[1]);
gsl_matrix* Ascopy=gsl_matrix_alloc(size[1],size[1]);
gsl_matrix* Rs=gsl_matrix_alloc(size[1],size[1]);
gsl_vector* bs=gsl_vector_alloc(size[1]);
gsl_vector* bstest=gsl_vector_alloc(size[1]);
gsl_vector* xs=gsl_vector_alloc(size[1]);

int n=size[1];

for (int i = 0; i<n; i++){
	for (int j=0; j<n; j++){
		gsl_matrix_set(As,j,i,1.0*rand_r(&rand_seed)/RAND_MAX);
		
	}
	gsl_vector_set(bs,i,1.0*rand_r(&rand_seed)/RAND_MAX);
}
gsl_matrix_memcpy(Ascopy,As);
printf("new matricies:\n A=\n");

for (int i =0; i<n;i++){
	for(int j=0; j<n;j++){
		printf("%10g   ",gsl_matrix_get(As,i,j));
	}
	printf("\n");
}

printf("b=\n");


for(int j=0; j<n;j++){
	printf("%10g \n  ",gsl_vector_get(bs,j));
}



printf("QR decomposition is done\n");
GS_decomp(As,Rs);
printf("Solve the equation QRx=b\n");
GS_solve(As,Rs,bs,xs);

printf("Checking that Ax=b\n Ax=\n");

gsl_blas_dgemv(CblasNoTrans,1.0,Ascopy,xs,0.0,bstest);


for(int j=0; j<n;j++){
	printf("%10g \n  ",gsl_vector_get(bstest,j));
}


printf("\n\n\n\n\n\n\n Inverse matrix production:\n");

gsl_matrix* Ai=gsl_matrix_alloc(n,n);
gsl_matrix* ident=gsl_matrix_alloc(n,n);

invert(As,Rs,Ai);

printf("checking that A^-1*A=I\n A^-1*A=\n ");

gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Ai,Ascopy,0.0,ident);

for (int i =0; i<n;i++){
	for(int j=0; j<n;j++){
		printf("%10g   ",gsl_matrix_get(ident,i,j));
	}
	printf("\n");
}


gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Ascopy,Ai,0.0,ident);
printf("A*A^-1=\n");
for (int i =0; i<n;i++){
	for(int j=0; j<n;j++){
		printf("%10g   ",gsl_matrix_get(ident,i,j));
	}
	printf("\n");
}


gsl_matrix_free(ident);
gsl_matrix_free(Ascopy);
gsl_matrix_free(Ai);
gsl_matrix_free(Rs);
gsl_vector_free(bs);
gsl_vector_free(bstest);
gsl_vector_free(xs);

clock_t t=clock();

FILE* stream_time_mine=fopen("out.my.times.txt","w");
FILE* stream_time_gsl=fopen("out.gsl.times.txt","w");
double time_taken;

int nmin=300;
int nmax=800;
int nstep=10;

for(int n=nmin;n<=nmax;n+=nstep){
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* R=gsl_matrix_alloc(n,n);
	for (int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			gsl_matrix_set(A,i,j,rand_r(&rand_seed)/RAND_MAX);
		 }
	}

	t=clock();
	GS_decomp(A,R);
	t=clock()-t;
	time_taken=((double) t)/CLOCKS_PER_SEC;
	fprintf(stream_time_mine,"%10i  %10g\n",n,time_taken);
	gsl_matrix_free(A);
	gsl_matrix_free(R);
}

nmin=2000;
nmax=4000;
nstep=50;

for(int n=nmin;n<=nmax;n+=nstep){
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_vector* R=gsl_vector_alloc(n);
	for (int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			gsl_matrix_set(A,i,j,rand_r(&rand_seed)/RAND_MAX);
		 }
	}

	t=clock();
	gsl_linalg_QR_decomp(A,R);
	t=clock()-t;
	time_taken=((double) t)/CLOCKS_PER_SEC;
	fprintf(stream_time_gsl,"%10i  %10g\n",n,time_taken);
	gsl_matrix_free(A);
	gsl_vector_free(R);
}




return 0;
}

