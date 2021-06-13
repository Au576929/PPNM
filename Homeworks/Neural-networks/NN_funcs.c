#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "mini_funcs.h"

double act_func(double x){

	return exp(-x*x);
}



void hidden_neuron(gsl_vector* Xs, gsl_vector* params, gsl_vector* Res){
	
	double a=gsl_vector_get(params,0);
	double b=gsl_vector_get(params,1);
	double w=gsl_vector_get(params,2);

	int Xs_count=Xs->size;
	gsl_vector* act_func_Xs=gsl_vector_alloc(Xs_count);
	gsl_vector_set_all(act_func_Xs,-1.0*a);
	
	gsl_blas_daxpy(1.0,Xs,act_func_Xs);
	gsl_vector_scale(act_func_Xs,1.0/b);
	double res=0;
	for (int i=0;i<Xs_count;i++){
		res=act_func(gsl_vector_get(act_func_Xs,i));
		gsl_vector_set(Res,i,res);
	}
	gsl_vector_scale(Res,w);

	gsl_vector_free(act_func_Xs);
}


double cost_func(gsl_vector* Res, gsl_vector* Ys){
        int Ys_count=Ys->size;
        double res=0;
        double holder=0;
        for (int i=0; i<Ys_count; i++){
 		assert(i<100);
                holder=(gsl_vector_get(Res,i)-gsl_vector_get(Ys,i));
                res+=holder*holder;
	}


return res;
}





double NeuralNetwork(gsl_vector* Xs, gsl_vector* Ys, gsl_vector* params, int neuron_count){
	
	int Xs_count=Xs->size;
	int param_count=params->size;
	int hidden_neuron_param_count=param_count/neuron_count; //int devision is not a problem as the two numbers have to be divisible
	

	gsl_vector* hidden_neuron_params=gsl_vector_calloc(hidden_neuron_param_count);
	gsl_vector* hidden_neuron_Res=gsl_vector_calloc(Xs_count);
	gsl_vector* Res=gsl_vector_calloc(Xs_count);
	
	for (int i=0; i<neuron_count; i++){
		for (int j=0;j<hidden_neuron_param_count; j++){
			gsl_vector_set(hidden_neuron_params, j ,gsl_vector_get(params, j+i*hidden_neuron_param_count));
		}
		hidden_neuron(Xs,hidden_neuron_params,hidden_neuron_Res);
		gsl_blas_daxpy(1.0,hidden_neuron_Res,Res);

	}
	
	double res=cost_func(Res,Ys);

	gsl_vector_free(hidden_neuron_params);
	gsl_vector_free(hidden_neuron_Res);
	gsl_vector_free(Res);

return res;
}




void NeuralNetworkTrainer(gsl_vector* Xs, gsl_vector* Ys, gsl_vector* params, int neuron_count, double acc){
//As my minimization algorithms minimize a function on the form double (gsl_vector) a new function is created on that form:

	double f_mini(gsl_vector* params_mini){

	return NeuralNetwork(Xs,Ys,params_mini,neuron_count);
	}

	minimize(f_mini,params,acc);

}


void interp_print(gsl_vector* Xs, gsl_vector* params, int neuron_count, char* name){

FILE* stream=fopen(name,"w");

int params_count=params->size;
int Xs_count=Xs->size;

gsl_vector* hidden_neuron_params=gsl_vector_alloc(params_count/neuron_count);
gsl_vector* Res_total=gsl_vector_calloc(Xs_count);
gsl_vector* Res=gsl_vector_calloc(Xs_count);


for (int i=0;i<neuron_count;i++){
	for (int j=0;j<params_count/neuron_count;j++){
		gsl_vector_set(hidden_neuron_params,j,gsl_vector_get(params,j+i*neuron_count));
	}
	hidden_neuron(Xs,hidden_neuron_params,Res);
	gsl_blas_daxpy(1.0,Res,Res_total);
}

for (int i=0;i<Xs_count;i++){
	fprintf(stream,"%10f  %10f\n",gsl_vector_get(Xs,i),gsl_vector_get(Res_total,i));
}



gsl_vector_free(hidden_neuron_params);
gsl_vector_free(Res_total);
gsl_vector_free(Res);

}



void interp_integ_print(gsl_vector* Xs, gsl_vector* params, int neuron_count, char* name){

FILE* stream=fopen(name,"w");

int params_count=params->size;
int Xs_count=Xs->size;

double xmin=gsl_vector_get(Xs,0);

double a=0;
double b=0;
double w=0;
int index=0;
double res=0;
double x=0;
for (int i=0;i<Xs_count;i++){
	x=gsl_vector_get(Xs,i);
	for (int j=0;j<neuron_count;j++){
		index=params_count/neuron_count*j;
        	a=gsl_vector_get(params,index);
		b=fabs(gsl_vector_get(params,index+1)); //b onlyenters the activation function as b², so it can turn out negative, it is however the uncertanty on the gaussian and as a result has to be positive. This is the only place where b factors not squared and so tha abs value ius just taken here.
		w=gsl_vector_get(params,index+2);
		res+=w*b*(erf((x-a)/b)-erf((xmin-a)/b))*sqrt(M_PI)/2;
	}
	fprintf(stream,"%10f  %10f\n",gsl_vector_get(Xs,i),res);
	res=0;
}

}


void interp_deriv_print(gsl_vector* Xs, gsl_vector* params, int neuron_count, char* name){

FILE* stream=fopen(name,"w");

int params_count=params->size;
int Xs_count=Xs->size;

double a=0;
double b=0;
double w=0;
int index=0; 
double res=0;
double x=0;
for (int i=0;i<Xs_count;i++){
	x=gsl_vector_get(Xs,i);
        for (int j=0;j<neuron_count;j++){
                index=params_count/neuron_count*j;
                a=gsl_vector_get(params,index);
                b=gsl_vector_get(params,index+1); 
                w=gsl_vector_get(params,index+2);
                res+=w*exp(-(x-a)*(x-a)/(b*b))*(-2*(x-a)/(b*b));
        }
        fprintf(stream,"%10f  %10f\n",x,res);
        res=0;
}

}
