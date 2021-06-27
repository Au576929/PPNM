#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 3e-16
#endif
static const double DELTA=sqrt(DBL_EPSILON);

void gradient(double f(gsl_vector* X),gsl_vector* X, gsl_vector* df){

int n=X->size;
double fx=f(X);
double fxpdx;
double dx;
double xi;
for (int i=0;i<n;i++){
	xi=gsl_vector_get(X,i);
	if (fabs(xi)<sqrt(DELTA)){
		dx=DELTA;
	} else {
		dx=fabs(xi)*DELTA;
	}
        gsl_vector_set(X,i,gsl_vector_get(X,i)+dx);
        fxpdx=f(X);
        gsl_vector_set(X,i,gsl_vector_get(X,i)-dx);
        gsl_vector_set(df,i,(fxpdx-fx)/dx);
}

}






void minimize(double f(gsl_vector* X), gsl_vector* X, double acc){

int dim = X->size;
int nsteps=0;
int nbad=0;
int ngood=0;

gsl_matrix* B=gsl_matrix_alloc(dim,dim);
gsl_vector* dF=gsl_vector_alloc(dim);
gsl_vector* dX=gsl_vector_alloc(dim);
gsl_vector* XpdX=gsl_vector_alloc(dim);
gsl_vector* dXpdX=gsl_vector_alloc(dim);
gsl_vector* Y=gsl_vector_alloc(dim);
gsl_vector* U=gsl_vector_alloc(dim);
gsl_vector* A=gsl_vector_alloc(dim);

gsl_matrix_set_identity(B);

gradient(f,X,dF);
double fX=f(X);
double fXpdX;

while (nsteps<10000){
	nsteps++;
	gsl_blas_dgemv(CblasNoTrans,-1.0,B,dF,0.0,dX);
	if (gsl_blas_dnrm2(dX)<DELTA*gsl_blas_dnrm2(X) || gsl_blas_dnrm2(dF)<acc){
		break;
	}
	
	double lambda=1;
	
	while (1){
		gsl_vector_memcpy(XpdX,X);
		gsl_vector_add(XpdX,dX);
		fXpdX=f(XpdX);
		double SdotdF;
		gsl_blas_ddot(dX,dF,&SdotdF);
		if (fXpdX<fX+0.01*SdotdF){
		ngood++;
		break;
		}
		if (lambda<DELTA){
			nbad++;
			gsl_matrix_set_identity(B);
			break;
		}
		lambda/=2;
		gsl_vector_scale(dX,0.5);
	}

	gradient(f,XpdX,dXpdX);
	gsl_vector_memcpy(Y,dXpdX);
	gsl_blas_daxpy(-1,dF,Y);
	gsl_vector_memcpy(U,dX);
	gsl_blas_dgemv(CblasNoTrans,1.0,B,Y,1.0,U);
	double SdotY;
	double UdotY;
	gsl_blas_ddot(dX,Y,&SdotY);
	if (fabs(SdotY)>1e-12){
		gsl_blas_ddot(U,Y,&UdotY);
		double gamma=UdotY/2.0/SdotY;
		gsl_blas_daxpy(-1.0*gamma,dX,U);
		gsl_blas_dger(1.0/SdotY,U,dX,B);
		gsl_blas_dger(1.0/SdotY,dX,U,B);
	}
	gsl_vector_memcpy(X,XpdX);
	gsl_vector_memcpy(dF,dXpdX);
	fX=fXpdX;
}








gsl_matrix_free(B);
gsl_vector_free(dF);
gsl_vector_free(dX);
gsl_vector_free(XpdX);
gsl_vector_free(dXpdX);
gsl_vector_free(Y);
gsl_vector_free(U);
gsl_vector_free(A);
}



void vec_extract(gsl_matrix* simplex, gsl_vector* vec, int j){
int d=simplex->size1;
for (int i=0;i<d;i++){
        gsl_vector_set(vec,i,gsl_matrix_get(simplex,i,j));
}
}







void reflection(gsl_vector* highest, gsl_vector* centroid,gsl_vector* reflected){
	int d=highest->size;
	for (int i=0;i<d;i++){
		double value=2*gsl_vector_get(centroid,i)-gsl_vector_get(highest,i);
		gsl_vector_set(reflected,i,value);
	}
}

void expansion(gsl_vector* highest, gsl_vector* centroid,gsl_vector* expanded){
        int d=highest->size;
        for (int i=0;i<d;i++){
                double value=3*gsl_vector_get(centroid,i)-2*gsl_vector_get(highest,i);
                gsl_vector_set(expanded,i,value);
        }
}


void contraction(gsl_vector* highest, gsl_vector* centroid,gsl_vector* contracted){
        int d=highest->size;
        for (int i=0;i<d;i++){
                double value=0.5*gsl_vector_get(centroid,i)-0.5*gsl_vector_get(highest,i);
                gsl_vector_set(contracted,i,value);
        }
}


void reduction(gsl_matrix* simplex, int lo){
        int d=simplex->size1;
        for (int i=0;i<d+1;i++){
		if (i!=lo){
			for (int j=0;j<d;j++){
				double value=0.5*(gsl_matrix_get(simplex,j,i)+gsl_matrix_get(simplex,j,lo));
				gsl_matrix_set(simplex,j,i,value);
			}
		}
        }
}



double distance(gsl_vector* vec1, gsl_vector* vec2){
	int d=vec1->size;
	double s=0;
	for (int i=0;i<d;i++){
		double vec1i=gsl_vector_get(vec1,i);
		double vec2i=gsl_vector_get(vec2,i);
		s+=(vec2i-vec1i)*(vec2i-vec1i);
	}
return sqrt(s);
}


double size(gsl_matrix* simplex){
	int d=simplex->size1;
	double s=0;
	gsl_vector* vec0=gsl_vector_alloc(d);
	gsl_vector* veci=gsl_vector_alloc(d);

	vec_extract(simplex,vec0,0);

	for(int i=1;i<d+1;i++){
		vec_extract(simplex,veci,i);
		double dist=distance(vec0,veci);
		if (dist>s){
			s=dist;
		}
	}

gsl_vector_free(vec0);
gsl_vector_free(veci);

return s;
}







void simplex_update(gsl_matrix* simplex, gsl_vector* f_values, int* hi, int* lo, gsl_vector* centroid){
int d=simplex->size1;
*hi=0;
*lo=0;
double highest=gsl_vector_get(f_values,0);
double lowest=highest;

for (int i=1;i<d+1;i++){
	double next=gsl_vector_get(f_values,i);
	if (next>highest){
		highest=next;
		*hi=i;
	} else if(next<lowest) {
		lowest=next;
		*lo=i;
	}
}

for (int i=0;i<d;i++){
	double s=0;
	for (int j=0;j<d+1;j++){
		if (j!=*hi){
			s+=gsl_matrix_get(simplex,i,j);
		}
	}
	gsl_vector_set(centroid,i,s/d);
}
}


void simplex_start(double f(gsl_vector*),gsl_matrix* simplex, gsl_vector* f_values, int* hi, int* lo, gsl_vector* centroid){
int d=simplex->size1;
gsl_vector* vec=gsl_vector_calloc(d);
for (int i=0;i<d+1;i++){
        vec_extract(simplex,vec,i);
        gsl_vector_set(f_values,i,f(vec));
}
simplex_update(simplex,f_values,hi,lo,centroid);
gsl_vector_free(vec);
}




int simplex_ringdown(double f(gsl_vector*),gsl_matrix* simplex,double simplex_size_goal){
int hi;
int lo;
int k=0;
int d=simplex->size1;

gsl_vector* centroid=gsl_vector_alloc(d);
gsl_vector* f_values=gsl_vector_alloc(d+1);
gsl_vector* p1=gsl_vector_alloc(d);
gsl_vector* p2=gsl_vector_alloc(d);
gsl_vector* holder=gsl_vector_alloc(d);

simplex_start(f,simplex,f_values,&hi,&lo,centroid);

while(size(simplex)>simplex_size_goal){
        simplex_update(simplex,f_values,&hi,&lo,centroid);
        vec_extract(simplex,holder,hi);
        reflection(holder,centroid,p1);
        double f_re=f(p1);
        if (f_re<gsl_vector_get(f_values,lo)){
                expansion(holder,centroid,p2);
                double f_ex=f(p2);
                if (f_ex<f_re){
                        for (int i=0;i<d;i++){
                                gsl_matrix_set(simplex,i,hi,gsl_vector_get(p2,i));
                                gsl_vector_set(f_values,hi,f_ex);
                        }
                } else {
                        for (int i=0;i<d;i++){
                                gsl_matrix_set(simplex,i,hi,gsl_vector_get(p1,i));
                                gsl_vector_set(f_values,hi,f_re);
                        }
                }
        } else {
                if (f_re<gsl_vector_get(f_values,hi)){
                        for(int i=0;i<d;i++){
                                gsl_matrix_set(simplex,i,hi,gsl_vector_get(p1,i));
                                gsl_vector_set(f_values,hi,f_re);
                        }
                } else {
                        vec_extract(simplex,holder,hi);
                        contraction(holder,centroid,p1);
                        double f_co=f(p1);
                        if (f_co<gsl_vector_get(f_values,hi)){
                                for(int i=0;i<d;i++){
                                        gsl_matrix_set(simplex,i,hi,gsl_vector_get(p1,i));
                                        gsl_vector_set(f_values,hi,f_co);
                                }
                        } else {
                                reduction(simplex,lo);
                                simplex_start(f,simplex,f_values,&hi,&lo,centroid);
                        }
                }
        }
        k++;
}
gsl_vector_free(holder);
gsl_vector_free(centroid);
gsl_vector_free(f_values);
gsl_vector_free(p1);
gsl_vector_free(p2);
return k;
}
















