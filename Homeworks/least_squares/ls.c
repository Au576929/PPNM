#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "ls_funcs.h"


int main(){

int n=7;
int m=4;

gsl_matrix* A=gsl_matrix_alloc(n,m);
gsl_matrix* R=gsl_matrix_alloc(m,m);
gsl_matrix* prodsi=gsl_matrix_alloc(m,m);
unsigned int rand_seed=1;

for (int i=0; i<n; i++){
	for (int j=0; j<m; j++){
		gsl_matrix_set(A,i,j,1.0*rand_r(&rand_seed)/RAND_MAX);
	}
}

FILE* data_stream=fopen("data.txt","r");
int i=0;
int N=9; //number of data points

gsl_matrix* data=gsl_matrix_alloc(N,3);

double X, Y;
//Data data[N]; //implementing struct

while(2== fscanf(data_stream, "%lg %lg", &X,&Y)){

gsl_matrix_set(data,i,0,X);
gsl_matrix_set(data,i,1,Y);
i++;
}

printf("X-data: Y-data:\n");
for(i=0;i<N;i++){
printf("%10g %10g\n",gsl_matrix_get(data,i,0),gsl_matrix_get(data,i,1));
}






printf("A=\n");
matrixprint(A);

GS_decomp(A,R);
printf("Q=\n");

matrixprint(A);


printf("R=\n");
matrixprint(R);

printf("Q^T*Q=\n");
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,prodsi);
matrixprint(prodsi);

printf("Thus QR decomposition seems to work on tall matricies.\n");

int Nf=2; //number of functions fitted to

gsl_vector* c=gsl_vector_alloc(Nf);

//As data doesn't have  uncertanties (YET!), and will be further modified below for fitting purposes the matrix will be printed to a seperate fil here
FILE* data_stream2=fopen("out.data.txt","w");

for (int i=0; i<N; i++){
        for (int j=0; j<2; j++){
                fprintf(data_stream2,"%10g", gsl_matrix_get(data,i,j));
        }
	fprintf(data_stream2," %10g \n",gsl_matrix_get(data,i,1)/20);
}





       	for (int j=0; j<N; j++){
           	gsl_matrix_set(data,j,1,log(gsl_matrix_get(data,j,1)));
		gsl_matrix_set(data,j,2,1./20);
	}


matrixprint(data);

gsl_matrix* Covy=gsl_matrix_alloc(Nf,Nf);


fit(data,&func,Nf,c,Covy);
printf("Covariance matrix=\n");
matrixprint(Covy);

printf("C=\n");

       for (int j=0; j<Nf; j++){
                printf("%10g   ",gsl_vector_get(c,j));
        }

printf("\nDelta C=\n");

       for (int j=0; j<Nf; j++){
                printf("%10g   ",sqrt(gsl_matrix_get(Covy,j,j)));
        }

double mht= -1*log(1./2)/gsl_vector_get(c,1); 
double mhtM=-1*log(1./2)/(gsl_vector_get(c,1)-sqrt(gsl_matrix_get(Covy,1,1)));
double mhtm=-1*log(1./2)/(gsl_vector_get(c,1)+sqrt(gsl_matrix_get(Covy,1,1))); 
 


printf("\n My halftime:\n %10g Max: %10g min: %10g",mht,mhtM,mhtm);
printf("\n Online: 3.63\n Thus they dont agree. However it should benoted that the other isotopes of Radium has markedly higher half times (14 days and 1600 years respectively) and so contamination could lead to the higher half-time observed.");

FILE* fit_stream = fopen("out.fit.txt","w");

       for (int j=1; j<=15; j++){
                fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)-gsl_vector_get(c,1)*j));
        }
 	fprintf(fit_stream,"\n\n\n\n\n\n\n\n");
      for (int j=1; j<=15; j++){
          fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)+sqrt(gsl_matrix_get(Covy,0,0))-(gsl_vector_get(c,1)+sqrt(gsl_matrix_get(Covy,1,1)))*j));
      }
        fprintf(fit_stream,"\n\n\n\n\n\n\n\n");

	for (int j=1; j<=15; j++){
          fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)-sqrt(gsl_matrix_get(Covy,0,0))-(gsl_vector_get(c,1)-sqrt(gsl_matrix_get(Covy,1,1)))*j));
      }
        printf("\n\n\n\n\n\n\n\n");








/*
gsl_matrix* Test=gsl_matrix_alloc(10,10);
gsl_matrix* Test2=gsl_matrix_alloc(10,10);
for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
                gsl_matrix_set(Test,i,j,1.0*rand_r(&rand_seed)/RAND_MAX);
        }
}

invert(Test,Test2);
*/



gsl_matrix_free(Covy);
gsl_vector_free(c);
gsl_matrix_free(prodsi);
gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(data);



return 0;
}

