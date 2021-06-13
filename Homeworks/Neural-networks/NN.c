#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "NN_funcs.h"
#include "print_funcs.h"



int main(){
//Generating a data set for interpolation
//function,. its derivative and integral functions are made here, so they are palced together
double func (double x){
return sin(x);
}

double dfunc(double x){
return cos(x);
}
double Func(double x){
return -cos(x);
}




double xmin=0;
double xmax=3*M_PI;


int x_count=50;
gsl_vector* Xs=gsl_vector_calloc(x_count);
gsl_vector* Ys=gsl_vector_calloc(x_count);
double dx=(xmax-xmin)/(x_count-1);
double x=xmin;
FILE* stream=fopen("out.tab.func.txt","w");

for (int i=0;i<x_count;i++){
	gsl_vector_set(Xs,i,x);
	gsl_vector_set(Ys,i,func(x));
	fprintf(stream,"%10f  %10f\n",x,gsl_vector_get(Ys,i));
	x+=dx;
}



int neuron_count=3;
int params_pr_neuron=3;


gsl_vector* params=gsl_vector_calloc(neuron_count*params_pr_neuron);

gsl_vector_set(params,0,1.5);
gsl_vector_set(params,1,1.5);
gsl_vector_set(params,2,1.5);

gsl_vector_set(params,3,4.5);
gsl_vector_set(params,4,1.5);
gsl_vector_set(params,5,-4.5);

gsl_vector_set(params,6,7.5);
gsl_vector_set(params,7,1.5);
gsl_vector_set(params,8,7.5);






/*
unsigned int rand_seed=3;

for (int i=0; i<neuron_count*params_pr_neuron; i++){
	gsl_vector_set(params,i,rand_r(&rand_seed)*5.0/RAND_MAX;
}*/
printf("Starting with data from a sinus wave, in the interval 0 to 3 pi.\n");
printf("The activation function used is a gaussian\n");
printf("Three neurons will be used to begin with\n");

printf("3 parameters are used per neuron: a,b and w in f((x-a)/b)*w\n");
printf("the parameters are ordered so that all parameters for the first neuron comes first, then the second and so on.\n");
printf("for each neuron the first parameter is a then b and finally w\n");


//double neural_network(gsl_vector* Xs, gsl_vector* Ys, gsl_vector* params, int neuron_count)
double cost=0;

printf("starting parameters:\n");
vector_print(params);

printf("Which gives a cost function of:\n");
cost=NeuralNetwork(Xs,Ys,params,neuron_count);
printf("%10f\n",cost);

double acc=0.05; //accuracy of minimization of costfunction in neural netowork
printf("neural network is trained (minimization over parameters is done):\n");
NeuralNetworkTrainer(Xs,Ys,params,neuron_count,acc);

printf("new parameters:\n");
vector_print(params);
printf("cost function with new params:\n");
cost=NeuralNetwork(Xs,Ys,params,neuron_count);
printf("%10f\n",cost);

int Xs_print_count=300;
gsl_vector* Xs_print=gsl_vector_calloc(Xs_print_count);

dx=(xmax-xmin)/(Xs_print_count-1);

x=xmin;

for (int i=0;i<Xs_print_count;i++){
	gsl_vector_set(Xs_print,i,x);
	x+=dx;
}

char name[80]="out.interp.txt";

interp_print(Xs_print,params,neuron_count,name);

// integration and differentation - part B

char name2[50]="out.interp.integ.txt";

interp_integ_print(Xs_print,params,neuron_count,name2);

char name3[50]="out.interp.deriv.txt";

interp_deriv_print(Xs_print,params,neuron_count,name3);

FILE* stream_deriv=fopen("out.tab.func.deriv.txt","w");
double res;
for (int i=0;i<Xs_print_count;i++){
	x=gsl_vector_get(Xs_print,i);
	res=dfunc(x);
	fprintf(stream_deriv,"%10f  %10f\n",x,res);
}



FILE* stream_integ=fopen("out.tab.func.integ.txt","w");

xmin=gsl_vector_get(Xs_print,0);
for (int i=0;i<Xs_print_count;i++){
        x=gsl_vector_get(Xs_print,i);
        res=Func(x)-Func(xmin);
        fprintf(stream_integ,"%10f  %10f\n",x,res);
}
fclose(stream);









gsl_vector_free(params);
gsl_vector_free(Xs);
gsl_vector_free(Ys);
gsl_vector_free(Xs_print);

return 0;
}
