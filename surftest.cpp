#include <stdio.h>
#include "SayeUtils.h"
#include "SayeAlgorithm.h"

int main()
{
    // Return variable
    double surface_area = 0;

    // Set up initial level set function
    //double coeff[10] = { 1, 4, 9, 0, 0, 0, 0, 0, 0, -1};
    double coeff[6] = { 1, 4, 0, 0, 0, -1};
    Phi2D phi0( coeff );

    // Initialize lists of functions
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ 0 };
    vector<int> ss1{ -1 };

    // Initialize box
    double xL[2] = {-2, -2};
    double xU[2] = {2, 2};
    Box U( xL, xU, 2 );

    bool S = true;
    int q = 5;
    
    // Get volume fraction
    Unit vof;
    R0R0 r0r0( &phi0 );

    surface_area = I( psis, ss, U, &r0r0, S, q );
    //double volume = I( psis, ss1, U, &vof, false, q );

    printf( "Surface Area: %4.10f\n", surface_area );
    //printf( "Volume: %4.10f\n", volume );

    return 0;
}

