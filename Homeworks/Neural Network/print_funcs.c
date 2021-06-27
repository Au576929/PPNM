#include <stdio.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void vector_print(gsl_vector* vec){

	for (int i=0; i<vec->size;i++){
		printf("%10f \n",gsl_vector_get(vec,i));
	}


}

