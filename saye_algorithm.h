#ifndef SAYEALGO
#define SAYEALGO

#include <vector>
using namespace std;

class F_params
{
public:
    F_params();
    F_params( const F_params & fp );
    
    F_params( vector<Psi>& n_psis, vector<int>& n_ss, F_params* n_fp,
              double n_x1, double n_x2, int n_k, int n_q, bool n_S );

    virtual double eval( Point x );

private:
    F_params * fp;
    vector<Psi> psis;
    vector<int> ss;
    double x1, x2;
    int k, q;
    bool S;
};

double gradient( Point x, Phi* phi );

#include "base_functions.h"

void init_quadrature( int q );
double corner_max_phi( const Box& U, Psi psi, Point xc );
double corner_max_grad( const Box& U, Psi psi, double g, int k );
double gauss_tensor_product( F_params * fp, const Box& U, int q );

double F( Point x, F_params * fp, vector<Psi>& psis, vector<int>& ss,
          double x1, double x2, int k, int q );          
double F_surf( Point x, F_params * fp, Psi phi, double x1, double x2, int k );          
double I( vector<Psi> psis, vector<int> ss, const Box& U, F_params * fp,
          bool S, int q, int depth = 0 );

double box_symdiff( int N, double xL, double xU, double yL, double yU, double a[] );
double quad_symdiff( double a[], double b[] );

#endif
