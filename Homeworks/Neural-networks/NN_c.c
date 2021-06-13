#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "NN_funcs.h"
#include "NN_funcs_partc.h"
#include "print_funcs.h"


int main(){

double Phi(double x,double y, double dy,double ddy){
return ddy+y;
}

int params_count=9;
gsl_vector* params=gsl_vector_calloc(params_count);

gsl_vector_set(params,0,1.63);
gsl_vector_set(params,1,2.47);
gsl_vector_set(params,2,3.27);

gsl_vector_set(params,3,7.71);
gsl_vector_set(params,4,2.38);
gsl_vector_set(params,5,3.02);

gsl_vector_set(params,6,3.54);
gsl_vector_set(params,7,-11.38);
gsl_vector_set(params,8,-2.33);

unsigned int rand_seed=3;
double factor;
for (int i=0;i<params_count;i++){
factor=1+rand_r(&rand_seed)*0.075/RAND_MAX;
gsl_vector_set(params,i,factor*gsl_vector_get(params,i));

}

int neuron_count=3;

double xmin=0;
double xmax=3*M_PI;
double x_boundary=5*M_PI/4;
double y_boundary=-1/sqrt(2);
double dy_boundary=-1/sqrt(2);
double acc=1e-5;


double res=NeuralNetwork_c(Phi,params,neuron_count,xmin,xmax,x_boundary,y_boundary,dy_boundary);
printf("\n\n\n\n PART C:\n\n");
printf("Using neural network to solve harmonic oscillator\n");
printf("Boundary conditions are placed at x=5*pi/4, y=-1/sqrt(2), y'=-1/sqrt(2), gives sinus solution.\n");
printf("solved from x=0, to x=3*pi, using 3 neurons, as in the previous parts.\n");

printf("initial parameters:\n");
vector_print(params);

printf("givees cost function: (using the square version from the homework description)\n");
printf("\n\n%10f\n",res);


printf("trainining the neural network.\n");
NeuralNetworkTrainer_c(Phi,params,neuron_count,xmin,xmax,x_boundary,y_boundary,dy_boundary,acc);



res=NeuralNetwork_c(Phi,params,neuron_count,xmin,xmax,x_boundary,y_boundary,dy_boundary);
printf("new parameters:\n");
vector_print(params);
printf("gives costfunction value:\n");
printf("\n\n%10f\n",res);


int xs_print=300;
gsl_vector* Xs=gsl_vector_calloc(xs_print);
double dx=(xmax-xmin)/(xs_print-1);
double x=xmin;
for (int i=0;i<xs_print;i++){
	gsl_vector_set(Xs,i,x);
	x+=dx;
}

char name[80]="out.harm.osc.txt";




interp_print(Xs,params,neuron_count,name);




gsl_vector_free(params);


return 0;
}
