#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int main(){

double func (double x){
return sin(x);
}


double xmin=0;
double xmax=3*M_PI;


int x_count=100;

double dx=(xmax-xmin)/(x_count-1);
double x=xmin;





for (int i=0;i<x_count;i++){
printf("%10f  %10f\n",x,func(x));
x+=dx;



}


return 0;
}
