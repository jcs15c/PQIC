#include <stdio.h>
#include <set>
#include <vector>
#include <cmath>
#include "saye_utils.h"
#include "saye_algorithm.h"
using namespace std;

#define MAX_DEPTH 16

double qs[10], ws[10]; 
void init_quadrature( int q )
{
    switch( q )
    {
    case(3):
        qs[0] = 0.11270166537925831148;
        qs[1] = 0.5;
        qs[2] = 0.88729833462074168851;

        ws[0] = 5.0/18.0;
        ws[1] = 4.0/9.0;
        ws[2] = 5.0/18.0;
        break;
    case(5):
        qs[0] = 0.5 * (1 - 0.90617984593866399278);
        qs[1] = 0.5 * (1 - 0.53846931010568309103);
        qs[2] = 0.5 * (1 + 0);
        qs[3] = 0.5 * (1 + 0.53846931010568309103);
        qs[4] = 0.5 * (1 + 0.90617984593866399278);
       
        ws[0] = 0.5 * 0.23692688505618908751;
        ws[1] = 0.5 * 0.47862867049936646804;
        ws[2] = 0.5 * 128.0 / 225.0;
        ws[3] = 0.5 * 0.47862867049936646804;
        ws[4] = 0.5 * 0.23692688505618908751;
        break;
    default:
        printf("DONT HAVE CORRECT QUADRATURE\n");
    }
}


double corner_max_phi( const Box& U, Psi psi, Point xc )
{
    int d = U.getD();
    double the_max = -1;
    double test_max;
    double vals[3];
    Point x;

    for( int i1 = 0; i1 < 2; i1++ )
    {
        vals[0] = U[0][i1];
        for( int i2 = 0; i2 < 2; i2++ )
        {
            vals[1] = U[1][i2];
            if( d == 3 )
            {
                for( int i3 = 0; i3 < 2; i3++ )
                {
                    vals[2] = U[2][i3];
                    x = Point( vals, d );
                    test_max = dabs( psi.eval(x) - psi.eval(xc) );
                    if( test_max > the_max )
                        the_max = test_max;
                }
            }
            else
            {
                x = Point( vals, d );
                test_max = dabs( psi.eval(x) - psi.eval(xc) );
                if( test_max > the_max )
                    the_max = test_max;
            }
        }
    }

    return the_max;
}

double corner_max_grad( const Box& U, Psi psi, double g, int k )
{
    int d = U.getD();
    double the_max = -1;
    double test_max;
    double vals[3];
    Point x;

    for( int i1 = 0; i1 < 2; i1++ )
    {
        vals[0] = U[0][i1];
        for( int i2 = 0; i2 < 2; i2++ )
        {
            vals[1] = U[1][i2];
            if( d == 3 )
            {
                for( int i3 = 0; i3 < 2; i3++ )
                {
                    vals[2] = U[2][i3];
                    x = Point( vals, d );
                    test_max = dabs( psi.eval_k(x, k) - g );
                    if( test_max > the_max )
                        the_max = test_max;
                }
            }
            else
            {
                x = Point( vals, d );
                test_max = dabs( psi.eval_k(x, k) - g );
                if( test_max > the_max )
                    the_max = test_max;
            }
        }
    }

    return the_max;
}

double gauss_tensor_product( F_params * fp, const Box& U, int q )
{
    init_quadrature( q );
    double I = 0;
    int d = U.getD();

    if( d == 1 )
        for( int i1 = 0; i1 < q; i1++ )
        {
            Point x;
            x.saye( U[0][0] + (U[0][1] - U[0][0])*qs[i1], 0 );
            I += ws[i1]*fp->eval( x );
        }

    if( d == 2 )
        for( int i1 = 0; i1 < q; i1++ )
            for( int i2 = 0; i2 < q; i2++ )
            {
                Point x;
                x.saye( U[0][0] + (U[0][1] - U[0][0])*qs[i1], 0 );
                x.saye( U[1][0] + (U[1][1] - U[1][0])*qs[i2], 1 );
                I += ws[i1]*ws[i2]*fp->eval( x );
            }

    if( d == 3 )
        for( int i1 = 0; i1 < q; i1++ )
            for( int i2 = 0; i2 < q; i2++ )
                for( int i3 = 0; i3 < q; i3++ )
                {
                    Point x;
                    x.saye( U[0][0] + (U[0][1] - U[0][0])*qs[i1], 0 );
                    x.saye( U[1][0] + (U[1][1] - U[1][0])*qs[i2], 1 );
                    x.saye( U[2][0] + (U[2][1] - U[2][0])*qs[i3], 2 );
                    I += ws[i1]*ws[i2]*ws[i3]*fp->eval( x );
                }

    return I;
}

F_params::F_params( )
{  }

F_params::F_params( const F_params & nfp )
{   
    fp = nfp.fp; 
    psis = nfp.psis;
    ss = nfp.ss;
    x1 = nfp.x1;
    x2 = nfp.x2;
    k = nfp.k;
    q = nfp.q;
    S = nfp.S;
}

F_params::F_params( vector<Psi>& n_psis, vector<int>& n_ss, 
                    F_params* n_fp, double n_x1, double n_x2, 
                    int n_k, int n_q, bool n_S )
{
    fp = n_fp; 
    psis = n_psis;
    ss = n_ss;
    x1 = n_x1;
    x2 = n_x2;
    k = n_k;
    q = n_q;
    S = n_S;
}
  
double F_params::eval( Point x )
{
    if( S == false )
        return F( x, fp, psis, ss, x1, x2, k, q );
    else
        return F_surf( x, fp, psis[0], x1, x2, k );
}

double F( Point x, F_params * fp, vector<Psi>& psis, vector<int>& ss,
          double x1, double x2, int k, int q )
{
    int bad_sign = 0;

    int d = x.getD();
    int n = psis.size(); 
    init_quadrature( q );

    Point xc;
    double I = 0;
    double L;

    set<double> R{ x1, x2 };
    set<double> R_temp;

    // Get roots for each psi
    for( int i = 0; i < n; i++ )
    {
        R_temp = psis[i].roots( x, x1, x2, k );
        for( auto & elem : R_temp )
            R.insert( elem );
    }

    // Turn into indexable array
    vector<double> R_list( R.begin(), R.end() );
    for( int j = 0; j < R_list.size() - 1; j++ )
    {
        L = R_list[j+1] - R_list[j];
        xc = x;
        xc.saye( 0.5*( R_list[j+1] + R_list[j] ), k );
        
        bad_sign = 0;
        for( int i = 0; i < n; i++ )
        {
            if( ss[i]*psis[i].eval( xc ) < 0 ) 
            {
                bad_sign = 1;
                break;
            }
        }
        if( !bad_sign )
        {
            for( int i = 0; i < q; i++ )
            {
                xc = x;
                xc.saye( R_list[j] + L*qs[i], k );
                I += L * ws[i]*fp->eval( xc );
            }
        }
    }

    return I;
}

double gradient( Point x, Phi* phi )
{
    double grad = 0;

    for( int j = 0; j < x.getD(); j++ )
        grad += phi->eval_k( x, j )*phi->eval_k( x, j );

    return sqrt( grad );
}

double F_surf( Point x, F_params * fp, Psi phi,
          double x1, double x2, int k )
{    
    int d = x.getD();
    set<double> R = phi.roots( x, x1, x2, k );
    double r;
    double grad = 0;

    if( R.size() == 1 )
    {
        r = *R.begin(); 
        x.saye( r, k );
        for( int j = 0; j < x.getD(); j++ )
            grad += phi.eval_k( x, j )*phi.eval_k( x, j );
        
        return fp->eval( x ) * sqrt( grad ) / 
                dabs( phi.eval_k( x, k ) );
    }
    else
        return 0;
}

double I( vector<Psi> psis, vector<int> ss, const Box& U, 
          F_params * fp, bool S, int q, int depth )
{
    int d = U.getD();
    int n = psis.size();
    int k = 0;

    Point tmp, xc;
    
    vector<Psi> psis_t;
    vector<int> ss_t;

    vector<double> g;
    vector<double> delta_k;
    double jacobian;
    // assert s == false etc.

    Box U1, U2, U3;

    F_params fpt;
    
    if( d == 1 )
        return F( tmp, fp, psis, ss, U[0][0], U[0][1], 0, q );

    double xcs[3] = {0.5*(U[0][0] + U[0][1]),
                     0.5*(U[1][0] + U[1][1]),
                     0.5*(U[2][0] + U[2][1])};
    xc = Point( xcs, d );

    double delta;

    for( int i = n-1; i >= 0; i-- )
    {
        delta = corner_max_phi( U, psis[i], xc );
        if( dabs( psis[i].eval( xc ) ) >= delta )
            if( ( ss[i]*psis[i].eval( xc ) > 0 ) ||
                ( ss[i]*psis[i].eval( xc ) >= 0 && S == false ) ) // Originally >=
            {
                n -= 1;
                psis.erase(psis.begin() + i);
                ss.erase(ss.begin() + i);
            }
            else
                return 0;
    }
    if( n == 0 )
        return U.volume() * gauss_tensor_product( fp, U, q );

    for( int j = 0; j < d; j++ )
    {
        if( dabs( psis[0].eval_k( xc, j ) ) > 
            dabs( psis[0].eval_k( xc, k ) ) )
            k = j;
    }

    for( int i = 0; i < n; i++ )
    {
        g.clear();
        delta_k.clear();
        jacobian = 0;
        for( int j = 0; j < d; j++ )
        {
            g.push_back( psis[i].eval_k( xc, j ) );
            delta_k.push_back( corner_max_grad( U, psis[i], g[j], j) );
            jacobian += (g[j] + delta_k[j])*(g[j] + delta_k[j]);
        }
        jacobian /= (g[k] - delta_k[k])*(g[k] - delta_k[k]);

        if( dabs(g[k]) > delta_k[k]  && jacobian < 20 )
        {
            psis_t.push_back( Psi( psis[i], U[k][0], k ) );
            ss_t.push_back( sgn( sign(g[k]), ss[i], S, -1 ) );

            psis_t.push_back( Psi( psis[i], U[k][1], k ) );
            ss_t.push_back( sgn( sign(g[k]), ss[i], S, 1 ) );
        }
        else
        {
            if( depth >= MAX_DEPTH )
                return U.volume() * fp->eval( xc );
             
            U.split( U1, U2 );
    
            double val1 = I( psis, ss, U1, fp, S, q, depth + 1 );
            double val2 = I( psis, ss, U2, fp, S, q, depth + 1 );
            return val1 + val2;
        }
    }

    fpt = F_params( psis, ss, fp, U[k][0], U[k][1], k, q, S );
    
    U.erase( U3, k );
    
    return I( psis_t, ss_t, U3, &fpt, false, q, depth + 1 );
}
