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
    double a0[6] = { 1.1, 0.9, 0.0, 0.0, 0, -1};
    double refs[6] = { PI/4, PI/4.0, 0, 0, 0, PI };
    flooding( a0, refs[5], 0.00001 );

    double an[6];
    double anp1[6];
    
    // calculate initial cost
    double init_cost = cost( a0, refs );
        
    double gamma = 0.1;
    double old_cost = init_cost, new_cost;
    
    // Do one step of gradient descent, readjust volume
    gradient_descent_step( a0, an, gamma, refs );
    for( int i = 0; i < 1000; i++ )
    {
        if( (i +1) % 100 == 0 )
            gamma = gamma / 2;
        new_cost = cost( an, refs );
        printf( "%d\t%3.10f\t%f\n", i, new_cost, gamma );
        gradient_descent_step( an, anp1, gamma, refs );
        for( int j = 0; j < 6; j++ )
            an[j] = anp1[j];
        
    }

    // calculate final cost
    double final_cost = cost( anp1, refs );

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

