#include "stdio.h"
#include "math.h"
#include "complex.h"

int main(){
complex z;


double Gamma=tgamma(5);
double Bessel=j1(0.5);

z=csqrt(-2);

printf("Gamafunction of 5: %f \n Besselfunction of 0.5: %f \n ",Gamma,Bessel);

printf("The squareroot of -2 is: %g+%gI \n",creal(z),cimag(z));

z=cexp(M_PI*I);

printf("exp(i*pi)=%g+%gI\n",creal(z),cimag(z));

z=cexp(I);

printf("exp(I)=%g+I%g\n",creal(z),cimag(z));

z=cpow(I,M_E);

printf("I^e=%g+%gI\n",creal(z),cimag(z));

z=cpow(I,I);

printf("I^I=%g+%gI\n",creal(z),cimag(z));

float x_float=1.f/9;
double x_double=1./9;
long double x_long_double=1.L/9;

printf("%.25g\n %.25lg \n %.25Lg",x_float,x_double,x_long_double);




return 0;
}
