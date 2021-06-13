double hidden_neuron_c(double x, gsl_vector* params);

double cost_func_c(double embed_func(double x),double xmin, double xmax, double acc, double eps
                , double x_boundary, double y_boundary, double dy_boundary, gsl_vector* params, int neuron_count);

double NeuralNetwork_c(double Phi(double x, double y, double dy, double ddy),gsl_vector* params, int neuron_count, double xmin, double xmax,
                        double x_boundary, double y_boundary, double dy_boundary);

double interp_c(double x, gsl_vector* params, int neuron_count);

double interp_deriv_c(double x, gsl_vector* params, int neuron_count);

double interp_dderiv_c(double x, gsl_vector* params, int neuron_count);

void NeuralNetworkTrainer_c(double Phi(double x, double y, double dy, double ddy),gsl_vector* params, int neuron_count, double xmin, double xmax,
                        double x_boundary,double y_boundary,double dy_boundary, double acc);
