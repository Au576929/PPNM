#include<stdio.h>
#include<pthread.h>
#include<math.h>
#include<stdlib.h>

struct params{
	int npoints;
	int count;
	unsigned int rand_seed;
};


void* bar(void* arg){
	struct params* p = (struct params*) arg;
	double x, y;
	for (int i=0;i<p->npoints;i++){ 
		x=rand_r(&(p->rand_seed));
		y=rand_r(&(p->rand_seed));
		x/=RAND_MAX;
		y/=RAND_MAX;
		if((pow(x,2)+pow(y,2))<=1){
			p->count++;
		}
	}
return NULL;
}



int main(){

int npoints=3e8;
int nthreads=3;

struct params thread1={.npoints=npoints/nthreads, .count=0, .rand_seed=1};
struct params thread2={.npoints=npoints/nthreads, .count=0, .rand_seed=2};
struct params thread3={.npoints=npoints/nthreads, .count=0, .rand_seed=3};


pthread_t threadx,thready,threadz;
pthread_create(&threadx,NULL,bar,(void*) &thread1);
pthread_create(&thready,NULL,bar,(void*) &thread2);
pthread_create(&threadz,NULL,bar,(void*) &thread3);

pthread_join(threadx,NULL);
pthread_join(thready,NULL);
pthread_join(threadz,NULL);

double counttot=thread1.count+thread2.count+thread3.count;
double pi_approx=4*counttot/npoints;
printf("Excersize A:\nPi has been approximated using 3e8 points, split in 3 threads:\n");
printf("pi is approx=%g\n\n\n",pi_approx);

return 0;
}
