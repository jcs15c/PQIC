using namespace std;

double cost_PLIC( double a[], double ref[] );
void gradient_descent_step_PLIC( double an[], double anp1[], double gamma, double refs[] );
double cost_partial_PLIC( double an[], double refs[], double grad[] );

double cost_PQIC( double a[], double a0[] );
double cost_partial_PQIC( double an[], double refs[], double grad[]);
void gradient_descent_step_PQIC( double an[], double anp1[], double gamma, double refs[] );

double cost_PPIC( double a[], double a0[] ); 
double cost_partial_PPIC( double cn[], double refs[], double grad[] );
void gradient_descent_step_PPIC( double cn[], double cnp1[], double gamma, double refs[] );

void flooding( double a[], double vof, double tol );

void poly_coefficients( double a[], double c[] ); 
void init_bases();
void init_grad_bases( Phi * phi0 );

// Probably dont work well
void cost_gradient_disc( double a[], double a0[], double del_a[] );
double int_partial_PQIC( int i, int j, double a[], double vof, double h, double tol );

// Probably dont work
void cost_gradient( double a[], double a0[], double del_a[] );
