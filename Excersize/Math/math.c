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

printf("1/9 as float:%.25g\n Thus float are precise to 9 decimal places. \n\n 1/9 as double: %.25lg \n Thus double are precise to 16 decimals.\n\n 1/9 as long double: %.25Lg\n Thus long doubles are precise to 19 decimals,\n",x_float,x_double,x_long_double);




return 0;
}
