#include<stdio.h>
#include<omp.h>
#include<math.h>
#include<stdlib.h>




void* bar(void* arg){
	int* count = (int*) arg;
	double count_temp=0;
	unsigned int seedx=1;
	unsigned int seedy=2;
	double x, y;
	for (int i=0;i<*count;i++){ 
		x=rand_r(&seedx);
		y=rand_r(&seedy);
		x/=RAND_MAX;
		y/=RAND_MAX;
		if((pow(x,2)+pow(y,2))<=1){
			count_temp++;
		}
	}
	*count=count_temp;
return NULL;
}



int main(){
int N;


for (N=1e3;N<1e8;N=N*2){
int x=N/3,y=N/3,z=N/3;

#pragma omp parallel sections
	{
	#pragma omp section
		{
		bar((void*)&x);
		}
	#pragma omp section
		{
		bar((void*)&y);
		}
	#pragma omp section
		{
		bar((void*)&z);
		}
	}

double pi_approx=4.0*(x+y+z)/N;
printf("For  %i random points:\n",N);
printf("pi is approx with omp =%g\n",pi_approx);
printf("Error is: %g\n\n",M_PI-pi_approx);
}
return 0;
}
