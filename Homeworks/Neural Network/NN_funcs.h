double act_func(double X);

void hidden_neuron(gsl_vector* Xs, gsl_vector* params, gsl_vector* Res);

double cost_func(gsl_vector* Res, gsl_vector* Ys);

double NeuralNetwork(gsl_vector* Xs, gsl_vector* Ys, gsl_vector* params, int neuron_count);

void NeuralNetworkTrainer(gsl_vector* Xs, gsl_vector* Ys, gsl_vector* params, int neuron_count,double acc);

void interp_print(gsl_vector* Xs, gsl_vector* params, int neuron_count, char* name);

void interp_integ_print(gsl_vector* Xs, gsl_vector* params, int neuron_count, char* name);

void interp_deriv_print(gsl_vector* Xs, gsl_vector* params, int neuron_count, char* name);
