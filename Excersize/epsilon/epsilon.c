#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(){

int i=400000;


while(i+1>i){
i++;
}

printf("Maximum integer = %i\n",i);
int out=1;
for (int i=1;i<i+1;i++){
out=i;
}

printf("Max integer: %i\n",out+1);

do{i++;} while(i<i+1);


printf("Max integer: %i\n",i);

while(i-1<i){
i--;
}

printf("Min integer = %i\n",i);

out=1;
for (int i=1;i>i-1;i--){
out=i;
}

printf("Min integer: %i\n",out-1);

do{i--;} while(i<i-1);


printf("Min integer: %i\n",i-1);


double xd=1;
long double xl=1;

while(1+xd!=1){xd/=2;}
xd*=2;


printf("min abs double: %g\n",xd);


float outf=1;

for (float xf=1;1+xf!=1;xf/=2){
outf=xf;
}
outf*=2;


printf("Min abs float: %g\n",outf);



do{xl/=2;} while(1+xl!=1);
xl*=2;

printf("Min abs long double : %lg\n",xl);










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


return 0;
}
