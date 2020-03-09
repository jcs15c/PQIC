#ifndef RBFENT
#define RBFNET

void get_training_data_para( int N, const char input_NN_fn[], const char output_NN_fn[] );
void get_training_data_quad( int N, const char input_NN_fn[], const char output_NN_fn[] );
void get_centers( int J2, const char input_NN_fn[], const char output_NN_fn[], const char centers_fn[] );
void get_weights( const char input_NN_fn[], const char output_NN_fn[], const char centers_fn[], const char weights_fn[] );
void eval_network( double refs[], double coeffs[], const char centers_fn[], const char weights_fn[] );
double RBF( double r, double sigma );
void poly_coefficients( double a[], double c[] );

#endif
