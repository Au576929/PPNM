#include <math.h>
#include <assert.h>

//integral with open set.
//the comparrison is done between the 2. and 4. order 
// this is done because reusing points recure an even nyumber of points to be used...
//As this function reuses points another function is needed to calculate those the first time hence the "2" in the name
double Integ2 (double func(double),double xmin, double xmax, double acc, double eps, double f2, double f3,  int nrec,double* Err){
assert(nrec<70000);
double f1=func(xmin+(xmax-xmin)/6); //the function evaluations that could not be reused
double f4=func(xmin+5*(xmax-xmin)/6);
double Q=(2*f1+f2+f3+2*f4)/6*(xmax-xmin); 	//result of 4. order
double q=(f1+f2+f3+f4)/4*(xmax-xmin);		//result of 2. order
double tol=acc+eps*fabs(Q);
double err=fabs(Q-q);

if (err<tol){
	*Err=sqrt((*Err)*(*Err)+err*err);
	return Q;
} else {
	nrec++;
	double Qi =Integ2(func,xmin,(xmax+xmin)/2,acc/sqrt(2.0),eps,f1,f2,nrec,Err);
	double Qii=Integ2(func,(xmax+xmin)/2,xmax,acc/sqrt(2.0),eps,f3,f4,nrec,Err);
	return Qi+Qii;
}
}


////////////////////////////////////////////////
/*
double IntegCC2 (double funcCC(double func(double),double),double xmin, double xmax, double acc, double eps, double f2, double f3, int nrec){
assert(nrec<1000000);
double f1=func(xmin+(xmax-xmin)/6); //the function evaluations that could not be reused

double Q=(2*f1+f2+f3+2*f4)/6*(xmax-xmin);       //result of 4. order
double q=(f1+f2+f3+f4)/4*(xmax-xmin);           //result of 2. order
double tol=acc+eps*fabs(Q);
double err=fabs(Q-q);

if (err<tol){
        return Q;
} else {
        nrec++;
        double Qi=Integ2(funcCC(func(x),x),xmin,(xmax+xmin)/2,acc/sqrt(2.0),eps,f1,f2,nrec);
        double Qii=Integ2(funcCC(func(x),x),(xmax+xmin)/2,xmax,acc/sqrt(2),eps,f3,f4,nrec);
        return Qi+Qii;
}
}

*/



//Clenshaw-Curtis transformation integrator
double IntegCC(double func(double),double acc, double eps, int nrec, double* Err){
double xmin=0;
double xmax=3.14159265;
//Err=0; //ensuring that the error double is set to zero in the outset
double Func(double x){

return sin(x)*func(cos(x));
}

double f2=Func(xmin+(xmax-xmin)/3);
double f3=Func(xmin+2*(xmax-xmin)/3);


return Integ2(Func,xmin,xmax,acc,eps,f2,f3,nrec,Err);
}


////////////////////////7
double Integ(double func(double),double xmin, double xmax, double acc, double eps, int nrec,double* Err){


double f2=0;
double f3=0;
//int i=0;
double xMax=xmax;
double xMin=xmin;

double Res=0;

if (isinf(xmax)==1 && isinf(xmin)==-1){
//	i=1;
	xMax=1;
	xMin=-1;
	
	double Func (double x){
		return func(x/(1.0-x*x))*(1.0+x*x)/((1.0-x*x)*(1.0-x*x));
	}

	f2=Func(xMin+(xMax-xMin)/3);
	f3=Func(xMin+2*(xMax-xMin)/3);
	
	Res=Integ2(Func,xMin,xMax,acc,eps,f2,f3,nrec,Err);
}

if (isinf(xmax)==1 && isinf(xmin)==0) {
//	i=2;
	xMin=0;
	xMax=1;
	double Func (double x){
                return func(xmin+x/(1.0-x))/((1.0-x)*(1.0-x));
        }

	f2=Func(xMin+(xMax-xMin)/3);
        f3=Func(xMin+2*(xMax-xMin)/3);


        Res=Integ2(Func,xMin,xMax,acc,eps,f2,f3,nrec,Err);
}


if (isinf(xmin)==-1 && isinf(xmax)==0){
//	i=3;
	xMin=0;
	xMax=1;
	
	double Func (double x){
                return func(xmax-(1-x)/x)/(x*x);
        }

	f2=Func(xMin+(xMax-xMin)/3);
        f3=Func(xMin+2*(xMax-xMin)/3);


        Res=Integ2(Func,xMin,xMax,acc,eps,f2,f3,nrec,Err);
}

if (isinf(xmin)==0 && isinf(xmax)==0){
//        i=3;
        xMin=xmin;
        xMax=xmax;

        double Func (double x){
                return func(x);
        }
	
	f2=Func(xMin+(xMax-xMin)/3);
        f3=Func(xMin+2*(xMax-xMin)/3);


        Res=Integ2(Func,xMin,xMax,acc,eps,f2,f3,nrec,Err);
}

/*
double Func(double x){
	double Result;

	if (i==1){
		Result=func(x/(1.0-x*x))*(1.0+x*x)/((1.0-x*x)*(1.0-x*x));
	} else if (i==2){
		Result=func(xmin+x/(1.0-x))/((1.0-x)*(1.0-x));
	} else if (i==3){
		Result=func(xmax-(1-x)/x)/(x*x);
	} else {
		Result=func(x);
	}

return Result;
}
	
return Integ2(Func,xMin,xMax,acc,eps,f2,f3,nrec,Err);
*/

return Res;
}
