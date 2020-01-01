#include <stdio.h>
#include "saye_utils.h"
#include "saye_algorithm.h"
#include "grad_descent.h"
#include <cmath>

double PI = atan(1) * 4;

//class M2x
int main()
{   
    // Return variable
    double the_cost = 0;

    // Set up initial level set function, adjust vof
    double a0[6] = { 0.0, 0.0, 0.0, -1.0, 0.0, 0.0};
    double refs[3] = { 0.5, 0.5, 2 };
    flooding( a0, refs[2], 0.00001 );

    double an[6];
    double anp1[6];
    
    // calculate initial cost
    double init_cost = cost_PLIC( a0, refs );
        
    double gamma = 0.01;
    double old_cost = init_cost, new_cost;
    
    // Do one step of gradient descent, readjust volume
    gradient_descent_step_PLIC( a0, an, gamma, refs );
    for( int i = 0; i < 1000; i++ )
    {
        //if( (i +1) % 100 == 0 )
        //    gamma = gamma / 2;
        new_cost = cost_PLIC( an, refs );
        printf( "%d\t%3.10f\t%f\n", i, new_cost, gamma );
        gradient_descent_step_PLIC( an, anp1, gamma, refs );
        for( int j = 0; j < 6; j++ )
            an[j] = anp1[j];
        
        flooding( an, refs[2], 0.00001 );

        //Phi2D phi0( an );
        //vector<Psi> psis{ Psi( &phi0 ) };
        //vector<int> ss_vol{ -1 };
        //double xL[2] = {-1, -1};
        //double xU[2] = {1, 1};
        //Box U( xL, xU, 2 );
        //Unit unit;
        //printf("\tvolume: %f\n", I( psis, ss_vol, U, &unit, false, 5 ) );

    }

    // calculate final cost
    double final_cost = cost_PLIC( anp1, refs );

    printf( "Init Cost:  %3.10f\n", init_cost );
    printf( "Final Cost: %3.10f\n", final_cost );

    printf( "Initial Coefficients: " );
    for( int i = 0; i < 6; i++ )
        printf( "% 2.8f, ", a0[i] );
    printf("\b\b\n");

    printf( "Final Coefficients:   " );
    for( int i = 0; i < 6; i++ )
        printf( "% 2.8f, ", anp1[i] );
    printf("\b\b\n");

    return 0;
}

