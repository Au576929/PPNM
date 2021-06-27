#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(){

int i=400000;


while(i+1>i){
i++;
}

printf("Maximum integer with while loop = %i\n",i);
int out=1;
for (int i=1;i<i+1;i++){
out=i;
}

printf("Max integer with for loop: %i\n",out+1);

do{i++;} while(i<i+1);


printf("Max integer with do while loop: %i\n",i);

printf("Max integer from limits.h: %i\n\n",INT_MAX);

while(i-1<i){
i--;
}

printf("Min integer with while loop= %i\n",i);

out=1;
for (int i=1;i>i-1;i--){
out=i;
}

printf("Min integer with for loop: %i\n",out-1);
i=1;
do{i--;} while(i-1<i);


printf("Min integer with do while loop: %i\n",i);
printf("Min integer from limits.h: %i\n\n", INT_MIN);









double xd=1;
long double xl=1;

while(1+xd!=1){xd/=2;}
xd*=2;


printf("min abs double while loop: %g\n DBL_EPSILON= %g\n",xd,DBL_EPSILON);

xd=1;

for(int i=0;1+xd!=1;xd/=2){}
printf("min abs double for loop= %g\n",xd*2);

xd=1;

do{xd/=2;}while(1+xd!=1);
printf("min abs double do while loop= %g\n\n",xd*2);

float xf=1;



while(1+xf!=1){xf/=2;}
xf*=2;
printf("Min abs float while loop: %g\nFLT_EPSILON=%g\n",xf,FLT_EPSILON);

xf=1;

for(int i=0;1+xf!=1;xf/=2){}
printf("min abs double for loop= %g\n",xf*2);

xf=1;

do{xf/=2;}while(1+xf!=1);
printf("min abs double do while loop= %g\n\n",xf*2);



while(1+xl!=1){xl/=2;}
xl*=2;

printf("Min abs long double : %Lg\nLDBL_EPSILON=%Lg\n",xl,LDBL_EPSILON);


xl=1;

for(int i=0;1+xl!=1;xl/=2){}
printf("min abs double for loop= %Lg\n",xl*2);

xl=1;

do{xl/=2;}while(1+xl!=1);
printf("min abs double do while loop= %Lg\n\n",xl*2);









int max=INT_MAX/2;
int counter=1;
float g=0;

while (counter<max+1){
g+=1.0f/counter;
counter++;
}

printf("sum up: %f\n",g);
g=0;

while (counter>0){
g+=1.0f/counter;
counter--;
}



printf("sum down: %f\n",g);
printf("The difference comes from when about due to 1/INT_MAX~1e-9<FLT_EPSILON. so when sum down is done the small numbers are calculated first, to which bigger and bigger numbers are added, thus keeping precicion and here getting to the bigger number. for sum up, the loss of percision causes the smaller number.\n\n");


double gd=0;

while (counter<max+1){
gd+=1.0/counter;
counter++;
}

printf("sum up: %g\n",gd);
gd=0;

while (counter>0){
gd+=1.0/counter;
counter--;
}


printf("sum down: %g\n",gd);
printf("Here 1/INT_MAX>DBL_EPSILON, so no loss of precision either way.");

return 0;
}
