#include<stdio.h>
#include<math.h>
#include<stdlib.h>



int equal(double a,double b, double tau, double eps){
if(abs(a-b)<tau || abs(a-b)/(abs(a)+abs(b))<eps/2){
return 1;
}
else{
return 0;
}
}


int main(){

double a=5;
double b=6;
double tau=1;
double epsilon=2;
int t=equal(a,b,tau,epsilon);


printf("for:\n a=%g\n b=%g\n tau=%g\n epsilon=%g\n",a,b,tau,epsilon);
printf("are they equal: %i\n",t);

b=7;
t=equal(a,b,tau,epsilon);
printf("for:\n a=%g\n b=%g\n tau=%g\n epsilon=%g\n",a,b,tau,epsilon);
printf("are they equal: %i\n",t);

tau=0;
epsilon=0;
t=equal(a,b,tau,epsilon);

printf("for:\n a=%g\n b=%g\n tau=%g\n epsilon=%g\n",a,b,tau,epsilon);
printf("are they equal: %i\n",t);

return 0;
}
