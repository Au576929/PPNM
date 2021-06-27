void plainmc(int dim, double f(int dim, double* x), double* a, double* b, int N, double* res, double* eps);

double corput(int,int);

void quasimc(int dim, double f(int dim, double* x), double* a, double* b, int N, double* res, double* eps);

double stratamc(
	int dim,
	double f(int dim,double*x),
	double*a,double*b,
	double acc,double eps,
	int n_reuse,double mean_reuse);

