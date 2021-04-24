#include<stdio.h>
#include<pthread.h>
#include<math.h>
#include<stdlib.h>




void* bar(void* arg){
	double* count = (double*) arg;
	unsigned int seedx=1;
	unsigned int seedy=2;
	double x, y;
	for (int i=0;i<1e8;i++){ 
		x=rand_r(&seedx);
		y=rand_r(&seedy);
		x/=RAND_MAX;
		y/=RAND_MAX;
		if((pow(x,2)+pow(y,2))<=1){
			*count=*count+1;
		}
	}
return NULL;
}



int main(){
double x=0,y=0,z=0;

pthread_t threadx,thready;

pthread_attr_t* attributes = NULL;

pthread_create(&threadx,attributes,bar,(void*)&x);
pthread_create(&thready,attributes,bar,(void*)&y);
bar((void*)&z);

void* returnvalue=NULL;
pthread_join(threadx,returnvalue);
pthread_join(thready,returnvalue);

double pi_approx=4*(x+y+z)/3e8;

printf("pi is approx=%g\n",pi_approx);

return 0;
}
