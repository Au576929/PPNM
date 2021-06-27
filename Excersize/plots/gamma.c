#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_gamma.h>


double Gamma(double x){
///single precision gamma function (Gergo Nemes, from Wikipedia)
if(x<0){
	return M_PI/sin(M_PI*x)/Gamma(1-x);
}
if(x<9){
	return Gamma(x+1)/x;
}
double lnGamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
return exp(lnGamma);
}




int main(){
	double x=-9.9;
	
	while (x<=10){
		printf("%10g %10g %10g %10g\n",x,Gamma(x),tgamma(x),gsl_sf_gamma(x));
		x+=0.25;
	}

return 0;
}
