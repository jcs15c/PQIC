#include <stdio.h>
#include "saye_utils.h"
#include "saye_algorithm.h"

int main()
{
    // Set up initial level set function
    //double coeff[10] = { 1, 4, 9, 0, 0, 0, 0, 0, 0, -1};
    double coeff[6] = { 1, 1, 0, 0, 0, -1};
    Phi2D phi0( coeff );

    // Initialize lists of functions
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    // Initialize box
    double xL[2] = {-1, -1};
    double xU[2] = {1, 1};
    Box U( xL, xU, 2 );

    bool S = false;
    int q = 5;
    
    // Get volume fraction

    Mx2 mx2;
    double moment_xx = I( psis, ss, U, &mx2, false, q );
    printf( "int{ xx }: %4.10f\n", moment_xx );
    
    My2 my2;
    double moment_yy = I( psis, ss, U, &my2, false, q );
    printf( "int{ yy }: %4.10f\n", moment_yy );
    
    Mxy mxy;
    double moment_xy = I( psis, ss, U, &mxy, false, q );
    printf( "int{ xy }: %4.10f\n", moment_xy );
    
    Cx cx;
    double moment_x = I( psis, ss, U, &cx, false, q );
    printf( "int{ x }: %4.10f\n", moment_x );
    
    Cy cy;
    double moment_y = I( psis, ss, U, &cy, false, q );
    printf( "int{ y }: %4.10f\n", moment_y );

    Unit vof;
    double volume = I( psis, ss, U, &vof, false, q );
    printf( "Volume: %4.10f\n\n", volume );
    
    return 0;
}

