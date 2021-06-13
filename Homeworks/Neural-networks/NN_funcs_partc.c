#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "mini_funcs.h"
#include "NN_funcs.h"
#include "num_funcs.h"

double hidden_neuron_c(double x, gsl_vector* params){
	
	double a=gsl_vector_get(params,0);
	double b=gsl_vector_get(params,1);
	double w=gsl_vector_get(params,2);

	double act_func_x=(x-a)/b;

	return w*act_func(act_func_x);
}


double interp_c(double x, gsl_vector* params, int neuron_count){

int params_count=params->size;

gsl_vector* hidden_neuron_params=gsl_vector_alloc(params_count/neuron_count);
double res=0;
for (int i=0;i<neuron_count;i++){
        for (int j=0;j<params_count/neuron_count;j++){
                gsl_vector_set(hidden_neuron_params,j,gsl_vector_get(params,j+i*neuron_count));
        }
        res+=hidden_neuron_c(x,hidden_neuron_params);
}

gsl_vector_free(hidden_neuron_params);

return res;
}



double interp_deriv_c(double x, gsl_vector* params, int neuron_count){


        int params_count=params->size;

        double a=0;
        double b=0;
        double w=0;
        int index=0;
        double res=0; 
        for (int j=0;j<neuron_count;j++){
                index=params_count/neuron_count*j;
                a=gsl_vector_get(params,index);
                b=gsl_vector_get(params,index+1); 
                w=gsl_vector_get(params,index+2);
                res+=w*exp(-(x-a)*(x-a)/(b*b))*(-2*(x-a)/(b*b));
        }

return res;
}



double interp_dderiv_c(double x, gsl_vector* params, int neuron_count){


        int params_count=params->size;

        double a=0;
        double b=0;
        double w=0;
        int index=0; 
        double res=0; 
        for (int j=0;j<neuron_count;j++){
                index=params_count/neuron_count*j;
                a=gsl_vector_get(params,index);
                b=gsl_vector_get(params,index+1); 
                w=gsl_vector_get(params,index+2);
                res+=w*exp(-(x-a)*(x-a)/(b*b))*(4.0*(x-a)*(x-a)/(b*b*b*b)-2.0/(b*b));
        }

return res;
}


double cost_func_c(double embed_func(double x),double xmin, double xmax, double acc, double eps
		, double x_boundary, double y_boundary, double dy_boundary, gsl_vector* params, int neuron_count){
	
	int nrec=0;
	double err=0;
	double holder=0;
	double res=0;
	holder=fabs(Integ(embed_func,xmin,xmax,acc,eps,nrec,&err));
	res+=holder;
	holder=fabs(interp_c(x_boundary,params,neuron_count)-y_boundary)*(xmax-xmin);
	res+=holder*holder;
	holder=fabs(interp_deriv_c(x_boundary,params,neuron_count)-dy_boundary)*(xmax-xmin);
	res+=holder*holder;
return res;
}





double NeuralNetwork_c(double Phi(double x, double y, double dy, double ddy),gsl_vector* params, int neuron_count, double xmin, double xmax,
			double x_boundary,double y_boundary,double dy_boundary){

	double embed_func(double x){
		double y=interp_c(x,params,neuron_count);
		double dy=interp_deriv_c(x,params,neuron_count);
		double ddy=interp_dderiv_c(x,params,neuron_count);
		
		return fabs(Phi(x,y,dy,ddy))*fabs(Phi(x,y,dy,ddy));
	}
	
	double acc=1e-2;
	double eps=1e-2;
	double res=cost_func_c(embed_func,xmin,xmax,acc,eps,x_boundary,y_boundary,dy_boundary,params,neuron_count);
return res;
}




void NeuralNetworkTrainer_c(double Phi(double x, double y, double dy, double ddy),gsl_vector* params, int neuron_count, double xmin, double xmax,
                        double x_boundary,double y_boundary,double dy_boundary, double acc){
//As my minimization algorithms minimize a function on the form double (gsl_vector) a new function is created on that form:

	double f_mini(gsl_vector* params_mini){

	return NeuralNetwork_c(Phi,params_mini,neuron_count,xmin,xmax,x_boundary,y_boundary,dy_boundary);
	}

	minimize(f_mini,params,acc);

}

