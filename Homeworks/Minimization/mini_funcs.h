void gradient(double (gsl_vector*),gsl_vector*, gsl_vector*);

int minimize(double (gsl_vector*), gsl_vector*,double);

int simplex_ringdown(double f(gsl_vector*),gsl_matrix* simplex,double simplex_size_goal);

void vec_extract(gsl_matrix* simplex, gsl_vector* vec, int j);

void simplex_update(gsl_matrix* simplex, gsl_vector* f_values, int* hi, int* lo, gsl_vector* centroid);

void simplex_start(double f(gsl_vector*),gsl_matrix* simplex, gsl_vector* f_values, int* hi, int* lo, gsl_vector* centroid);

double size(gsl_matrix* simplex);

double distance(gsl_vector* vec1, gsl_vector* vec2);

void reduction(gsl_matrix* simplex, int lo);

void contraction(gsl_vector* highest, gsl_vector* centroid,gsl_vector* contracted);

void expansion(gsl_vector* highest, gsl_vector* centroid,gsl_vector* expanded);

void reflection(gsl_vector* highest, gsl_vector* centroid,gsl_vector* reflected);
