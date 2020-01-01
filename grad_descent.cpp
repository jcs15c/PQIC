#include <stdio.h>
#include <cmath>
#include "saye_utils.h"
#include "saye_algorithm.h"
#include "grad_descent.h"

Mx2 mx2;
My2 my2;
Mxy mxy;
Cx cx;
Cy cy;

Gx4 gx4;
Gx3 gx3;
Gx2 gx2;
Gx  gx;

Gx3y gx3y;
Gx2y gx2y;
Gxy  gxy;
Gy   gy;

Gx2y2 gx2y2;
Gxy2  gxy2;
Gy2   gy2;

Gxy3 gxy3;
Gy3  gy3;

Gy4  gy4;

F_params * bases[5];
F_params * grad_bases[5];
F_params * d_bases[5][5];

// cost function of QLIC 
//  a is coefficients for phi
//  a0 is reference data
double cost( double a[], double a0[] )
{
    init_bases();
    double the_cost = 0;

    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    double xL[2] = {-1, -1};
    double xU[2] = { 1,  1};
    Box U( xL, xU, 2 );
    int q = 5;
    
    for( int i = 0; i < 5; i++ )
    {
        double temp = I( psis, ss, U, bases[i], false, q ) / a0[5];
        //printf("\ttemp: %phi2d phi0( a );f\n", temp);
        the_cost += ( temp - a0[i] )*( temp - a0[i] );
    }

    return the_cost;
}

// Attempted analytical gradient function for QLIC method
//  probably doesn't work
//  del_a is the return value for the gradient
void cost_gradient( double a[], double a0[], double del_a[] )
{
    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss_vol{ -1 };
    vector<int> ss_surf{ 0 };

    double xL[2] = {-1, -1};
    double xU[2] = { 1,  1};
    Box U( xL, xU, 2 );
    int q = 5;

    init_grad_bases( &phi0 );
    init_bases();
    
    double vol_ints[5];
    double surf_ints[5][5];
    double surf_grads[5];

    // Set up integrals
    for( int i = 0; i < 5; i++ )
        vol_ints[i] = I( psis, ss_vol, U, bases[i], false, q ) / a0[5];

    for( int i = 0; i < 5; i++ )
        surf_ints[i][i] = I( psis, ss_surf, U, d_bases[i][i], true, q );

    for( int i = 0; i < 5; i++ )
        for( int j = 0; j < i; j++ )
            surf_ints[j][i] = surf_ints[i][j] =
                    I( psis, ss_surf, U, d_bases[i][j], true, q );
   
    for( int i = 0; i < 5; i++ )
        surf_grads[i] = I( psis, ss_surf, U, grad_bases[i], true, q );
    
    // Calculate gradient
    for( int j = 0; j < 5; j++ )
    {
        del_a[j] = 0;
        for( int i = 0; i < 5; i++ )
            del_a[j] += 2*( vol_ints[i] - a0[i] ) *
                          ( a0[5]*surf_ints[i][j] - vol_ints[i]*surf_grads[j] ) / a0[5] / a0[5];
    }
}

// More correct, less accurate QLIC gradient function
void cost_gradient_disc( double a[], double refs[], double del_a[] )
{
    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss_vol{ -1 };
    vector<int> ss_surf{ 0 };

    double xL[2] = {-1, -1};
    double xU[2] = { 1,  1};
    Box U( xL, xU, 2 );
    int q = 5;

    init_bases();
    
    double vol_ints[5];
    double partials[5][5];


    // Set up integrals
    for( int i = 0; i < 5; i++ )
        vol_ints[i] = I( psis, ss_vol, U, bases[i], false, q ) / refs[5];

    
    // Set up partials
    for( int i = 0; i < 5; i++ )
        for( int j = 0; j < 5; j++ ) 
            partials[i][j] = partial_approx( i, j, a, refs[5], 0.00001, 0.0001 );
    

    // Calculate gradient
    for( int j = 0; j < 5; j++ )
    {
        del_a[j] = 0;
        for( int i = 0; i < 5; i++ )
            del_a[j] += 2*( vol_ints[i] - refs[i] ) * partials[i][j] / refs[5];
    }
}

// Approximate ith term of cost function with respecet to j
double partial_approx( int i, int j, double a[], double vof, double h, double tol )
{
    double partial = 0;

    double original = a[j];    
    vector<int> ss{ -1 };
    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };

    double xL[2] = {-1, -1};
    double xU[2] = { 1,  1};
    Box U( xL, xU, 2 );
    int q = 5;

    init_bases();
    
    a[j] = original + 2*h;
    flooding( a, vof, tol );
    phi0 = Phi2D( a );
    psis[0] = Psi( &phi0 );
    partial += -1*I( psis, ss, U, bases[i], false, q ) / vof;
    
    a[j] = original + h;
    flooding( a, vof, tol );
    phi0 = Phi2D( a );
    psis[0] = Psi( &phi0 );
    partial += 8*I( psis, ss, U, bases[i], false, q ) / vof;
    
    a[j] = original - h;
    flooding( a, vof, tol );
    phi0 = Phi2D( a );
    psis[0] = Psi( &phi0 );
    partial += -8*I( psis, ss, U, bases[i], false, q ) / vof;

    a[j] = original - 2*h;
    flooding( a, vof, tol );
    phi0 = Phi2D( a );
    psis[0] = Psi( &phi0 );
    partial += I( psis, ss, U, bases[i], false, q ) / vof;

    partial = partial / 12 / h;
}

void gradient_descent_step( double an[], double anp1[], double gamma, double refs[] )
{
    double grad_an[5];
    
    cost_gradient_disc( an, refs, grad_an );

    for( int i = 0; i < 5; i++ )
        anp1[i] = an[i] - gamma*grad_an[i];
    

    flooding( anp1, refs[5], 0.0001 );
}

void flooding( double a[], double vof, double tol )
{
    vector<int> ss{ -1 };
    double xL[2] = {-1, -1};
    double xU[2] = { 1,  1};
    Box U( xL, xU, 2 );
    int q = 5;
    Unit unit;

    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };

    double low, high, end;
    double start = I( psis, ss, U, &unit, false, q );

    double test_a5, test_vol;
    // get bounds for bisection method
    if( start < vof ) // Need to decrease a5
    {
        high = a[5];
        low = high - 1;
        while( true )
        {
            a[5] = low;
            phi0 = Phi2D( a );
            psis[0] = Psi( &phi0 );
            
            end = I( psis, ss, U, &unit, false, q );
            
            if( end > vof )
                break;
            else
                low -= 1;
        }
    }
    else // need to increase a5
    {
        low = a[5];
        high = low + 1;
        while( true )
        {
            a[5] = high;
            phi0 = Phi2D( a );
            psis[0] = Psi( &phi0 );

            end = I( psis, ss, U, &unit, false, q );

            if( end < vof )
                break;
            else
                high += 1;
        }
    }

    // Bisection method part
    while( dabs( test_vol - vof ) > tol )
    {
        a[5] = 0.5*(low + high);
        
        phi0 = Phi2D( a );
        psis[0] = Psi( &phi0 );

        test_vol = I( psis, ss, U, &unit, false, q );
        
        if( test_vol > vof )
            low = a[5];
        else
            high = a[5];
    }

    test_vol;
}

void init_bases()
{
    bases[0] = &mx2;
    bases[1] = &my2;
    bases[2] = &mxy;
    bases[3] = &cx;
    bases[4] = &cy;
}

void init_grad_bases( Phi * phi0 )
{   gx4.phi0 = phi0;
    gx3.phi0 = phi0;
    gx2.phi0 = phi0;
    gx3y.phi0 = phi0;
    gx2y.phi0 = phi0;
    gxy.phi0 = phi0;
    gx2y2.phi0 = phi0;
    gxy2.phi0 = phi0;
    gy2.phi0 = phi0;
    gxy3.phi0 = phi0;
    gy3.phi0 = phi0;
    gy4.phi0 = phi0;
    gx.phi0 = phi0;
    gy.phi0 = phi0;

    d_bases[0][0] = &gx4;
    d_bases[0][1] = d_bases[1][0] = d_bases[2][2] = &gx2y2;
    d_bases[0][2] = d_bases[2][0] = &gx3y;
    d_bases[0][3] = d_bases[3][0] = &gx3;
    d_bases[0][4] = d_bases[4][0] = d_bases[3][2] = d_bases[2][3] = &gx2y;
    d_bases[1][1] = &gy4;
    d_bases[1][2] = d_bases[2][1] = &gxy3;
    d_bases[1][3] = d_bases[3][1] = d_bases[4][2] = d_bases[2][4] = &gxy2;
    d_bases[1][4] = d_bases[4][1] = &gy3;
    d_bases[3][3] = &gx2;
    d_bases[3][4] = d_bases[4][3] = &gxy;
    d_bases[4][4] = &gy2;

    grad_bases[0] = &gx2;
    grad_bases[1] = &gy2;
    grad_bases[2] = &gxy;
    grad_bases[3] = &gx;
    grad_bases[4] = &gy;
}
