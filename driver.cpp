#include <stdio.h>
#include "saye_utils.h"
#include "saye_algorithm.h"
#include "grad_descent.h"
#include "level_set_plot.h"

int main()
{
    double xL[2] = {2, 1};
    double xU[2] = {4, 5};
    Box U( xL, xU, 2 );
    
    // Set up initial level set function
    //double coeff[10] = { 1, 4, 9, 0, 0, 0, 0, 0, 0, -1};
    // alpha, h, k, alpha
    //double params[4] = { 1, 0, 1, 0 };
    double coeff_t[6] = {0, 0, 0, -1, -1, 6};
    double refs[6];
    double refs_t[6];
    //get_refs( U, coeff, refs );
    //transform_refs_to_unit( refs, refs_t, U );
    double coeff[6] = {0, 0, 0, -2, -4, 3};
    transform_coeffs_to_nonunit( coeff, coeff_t, U );
    //poly_coefficients( coeff, params ); 
    //for( int i = 1; i < 6; i++ )
    //    coeff[i] = coeff[i] / coeff[0];
    //coeff[0] = 1.0;

    Phi2D phi0( coeff );

    // Initialize lists of functions
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    // Initialize box

    bool S = false;
    int q = 5;
    
    // Get volume fraction
    Unit vof;
    double volume = I( psis, ss, U, &vof, false, q );

    Mx2 mx2;
    double moment_xx = I( psis, ss, U, &mx2, false, q );
    printf( "int{ xx }: %4.10f\n", moment_xx/volume );
    
    My2 my2;
    double moment_yy = I( psis, ss, U, &my2, false, q );
    printf( "int{ yy }: %4.10f\n", moment_yy/volume );
    
    Mxy mxy;
    double moment_xy = I( psis, ss, U, &mxy, false, q );
    printf( "int{ xy }: %4.10f\n", moment_xy/volume );
    
    Cx cx;
    double moment_x = I( psis, ss, U, &cx, false, q );
    printf( "int{ x }: %4.10f\n", moment_x/volume );
    
    Cy cy;
    double moment_y = I( psis, ss, U, &cy, false, q );
    printf( "int{ y }: %4.10f\n", moment_y/volume );

    printf( "Volume: %4.10f\n\n", volume );
    
    return 0;
}

