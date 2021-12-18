#include <stdio.h> 
#include "saye_utils.h"
#include "saye_algorithm.h"
#include "grad_descent.h"
#include "level_set_plot.h"
#include <cmath>

double PI = atan(1) * 4;

//class M2x
int main()
{   
    // Return variable
    double the_cost = 0;

    // Set up initial level set function, adjust vof
    //double a0[6] = { 0.0, 0.0, 0.0, -0.3, 0.1, 0.0};
    //double refs_PLIC[3] = { 0.5, 0.5, 2.0 };
    //flooding( a0, refs_PLIC[2], 0.00001 );

    double p0[4] = {0, 0, 0, 0.1};
    double refs[6] = { 0.8239578526,
                       0.9320568423,
                      -0.1567717007,
                       0.6754312145,
                       0.3945994185,
                       2.9036603218 };

    double a[6];
    poly_coefficients( a, p0 );
    plot_2D( a, refs[3], refs[4], "before.png" );

    double pn[4];
    double pnp1[4];

    // calculate initial cost
    //double init_cost = cost_PLIC( a0, refs );
    double init_cost = cost_PPIC( p0, refs );

    double gamma = 0.1;
    double old_cost = init_cost, new_cost;

    // Do one step of gradient descent, readjust volume
    //gradient_descent_step_PLIC( a0, an, gamma, refs_PLIC );
    gradient_descent_step_PPIC( p0, pn, gamma, refs );
    for( int i = 0; i < 1000; i++ )
    {
        //if( (i +1) % 100 == 0 )
        //    gamma = gamma / 2;

        new_cost = cost_PPIC( pn, refs );

        printf( "%d\t%3.10f\t%f\n", i, new_cost, gamma );

        gradient_descent_step_PPIC( pn, pnp1, gamma, refs );

        for( int j = 0; j < 4; j++ )
            pn[j] = pnp1[j];

        if( old_cost == new_cost ) break;
        old_cost = new_cost; 

        //poly_coefficients( a, pnp1 );
        //Phi2D phi0( a );
        //vector<Psi> psis{ Psi( &phi0 ) };
        //vector<int> ss_vol{ -1 };
        //double xL[2] = {-1, -1};
        //double xU[2] = {1, 1};
        //Box U( xL, xU, 2 );
        //Unit unit;
        //printf("\tvolume: %f\n", I( psis, ss_vol, U, &unit, false, 5 ) );
    }

    // calculate final cost
    //double final_cost = cost_PLIC( anp1, refs );
    double final_cost = cost_PPIC( pnp1, refs );

    poly_coefficients( a, pnp1 );
    plot_2D( a, refs[3], refs[4], "after.png" );

    printf( "Init Cost:  %3.10f\n", init_cost );
    printf( "Final Cost: %3.10f\n", final_cost );

    printf( "Initial Coefficients: " );
    for( int i = 0; i < 4; i++ )
        printf( "% 2.8f, ", p0[i] );
    printf("\b\b\n");

    printf( "Final Coefficients:   " );
    for( int i = 0; i < 4; i++ )
        printf( "% 2.8f, ", pnp1[i] );
    printf("\b\b\n");

    return 0;
}
