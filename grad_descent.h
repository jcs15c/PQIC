using namespace std;

double cost( double a[], double a0[] );
void cost_gradient( double a[], double a0[], double del_a[] );
void cost_gradient_disc( double a[], double a0[], double del_a[] );
void gradient_descent_step( double an[], double anp1[], double gamma, double refs[] );
void flooding( double a[], double vof, double tol );
double partial_approx( int i, int j, double a[], double vof, double h, double tol );

void init_bases();
void init_grad_bases( Phi * phi0 );
