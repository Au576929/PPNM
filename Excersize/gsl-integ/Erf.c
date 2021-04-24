
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>

double f (double x,void* params) {
//  double z = *(double*)params; 
  double f = 2/sqrt(M_PI)*exp(-pow(x,2));
  return f;
}


double mygslinteg(double x){
	gsl_function F;
	F.function=&f;
	int limit=999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0,b=x,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){		
	for (double x=-2;x<=2;x+=0.1){
		printf("%10g  %10g\n",x,mygslinteg(x));
	}	
return 0;
}
