#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
//#include <gsl/gsl_integ.h>

int binsearch(int n, double* x, double input){
assert(n>1 && input<=x[n-1] && input>=x[0]);
int i=0, j=n-1, m=0; // setting up for binary search for intervals
while(j-i>1){
	m=(i+j)/2;

	if(input<x[m]){
		j=m;
	}else{
		i=m;
	}
}
return i;
}


double lininterp(int n, double* x,double* y,double input){

int i = binsearch(n,x,input);

assert(x[i+1]>x[i]);
double output=y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(input-x[i]);

return output;
}


double lininterp_integ(int n, double* x,double* y,double input){

int i = binsearch(n,x,input);

double output=0;
int p=1;
while(p<=i){
	output+=y[p-1]*(x[p]-x[p-1])+0.5*(y[p]-y[p-1])*(x[p]-x[p-1]); //calculating area between the points x[p-1] and x[p], first the rectangle under tyh line and then the triangle on top
	p++;
}


double yinput=lininterp(n,x,y,input);

output+=y[i]*(input-x[i])+0.5*(yinput-y[i])*(input-x[i]);

return output;
}





int main(){

int n=20;

double x[n],y[n],Y[n];

FILE* my_out_stream=fopen("out.xydata.txt","w");

for(int i=0;i<n;i++){
x[i]=i/2.0;
y[i]=sin(x[i]);
Y[i]=-cos(x[i])+1;
fprintf(my_out_stream,"%10g %10g \n",x[i],y[i]);
}

FILE* my_out_stream_integ_points=fopen("out.integ.exact.txt","w");
for (int i=0; i<n;i++){
	fprintf(my_out_stream_integ_points,"%10g	%10g\n",x[i],Y[i]);
}


int z=0;
double fineness=10;


while(z<=fineness*x[n-1]){
printf("%10g	%10g\n",z/fineness,lininterp(n,x,y,z/fineness));
z++;
}


FILE* my_out_stream_integ=fopen("out.xyinteg.txt","w");

z=0;

while(z<=fineness*x[n-1]){
	fprintf(my_out_stream_integ,"%10g	%10g\n",z/fineness,lininterp_integ(n,x,y,z/fineness));
	z++;
}

double xt[n/2], yt[n/2];
for (int i=0;i<n/2;i++){
	yt[i]=y[2*i];
	xt[i]=x[2*i];
}


FILE* my_out_stream_integ_half=fopen("out.xyinteg.half.txt","w");

z=0;

while(z<=fineness*xt[n/2-1]){
	fprintf(my_out_stream_integ_half,"%10g	%10g\n",z/fineness,lininterp_integ(n/2,xt,yt,z/fineness));
	
//	fprintf(my_out_stream_integ_half,"xt[0]=%10g  xt[n/2-1]=%10g z=%10g\n",xt[0],xt[n/2-1],z/fineness);

	z++;
}


return 0;
}
