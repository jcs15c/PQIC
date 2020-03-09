#include <vector>
#include <armadillo>
#include <stdio.h>

#include "saye_utils.h"
#include "saye_algorithm.h"
#include "level_set_plot.h"
#include "rbf_network.h"
#include "grad_descent.h"


int main()
{
    char energy_str[50];
    char filename_str[50];

    double box_size = 3;
    double true_coeffs[6] = {1, 1, 0, -box_size, -box_size, (box_size/2)*(box_size/2)}; 
    double this_symdiff_err;
    
    // Domain
    
    // Guess type (parabolic)
    char guess_type = 'l';

    // maximum number of iterations
    int max_iter = 500;

    // error tolerances
    double error_tols[4] = {1e-7, 1e-7, 1e-7, 1e-4};

    ////////////// test /////////////////
    /*
    double txL[2] = { -1, -1 };
    double txU[2] = {  1,  1 };
    Box tU( txL, txU, 2 );
    double test_coeffs[6] = {0, 0, 0, 0, 0, 0};
    double test_refs[6] = {1.0/3.0, 1.0/3.0, 0.25, -0.5, -0.5, 1};
    parabolic_MOF( test_refs, tU, 'l', error_tols, test_coeffs, max_iter );
    return 0;
    */
    /////////////////////////////////////


    //parabolic_MOF( refs, U, guess_type, error_tols, final_coeffs, max_iter );

    // Do 32 x 32 test
    int k;
    int N = 3;
    
    double final_coeffs[N*N][6];
    double init_guesses[N*N][6];
    double refs[N*N][6];
   
    double xL[2] = { 0,  0};
    double xU[2] = { 1,  1};
    Box U( xL, xU, 2 );
   
    //get_refs( U, true_coeffs, refs[0] );
    //parabolic_MOF( refs[0], U, guess_type, error_tols, final_params[0], max_iter );

    // return 0;
    // Get reference data
    for( int i = 0; i < N; i++ )
    {
        for( int j = 0; j < N; j++ )
        {
            k = N*i + j;
            
            xL[0] = i * box_size / N;
            xL[1] = j * box_size / N;

            xU[0] = (i+1) * box_size / N;
            xU[1] = (j+1) * box_size / N;
    
            U = Box( xL, xU, 2 );
            
            get_refs( U, true_coeffs, refs[k] );
        }
    }

    k = 15;

    xL[0] = 3;
    xL[1] = 3;

    xU[0] = 4;
    xU[1] = 4;
    
    U = Box( xL, xU, 2 );
            
    //parabolic_MOF( refs[k], U, guess_type, error_tols, final_coeffs[k], max_iter );
    this_symdiff_err = quad_symdiff( true_coeffs, final_coeffs[k] ); 
            
    sprintf(filename_str, "rbf_data2/MOF_%d.png", k);
    sprintf(energy_str, "Sym Diff: %.6f", this_symdiff_err);
            
    //plot_2D_reference( true_coeffs, final_coeffs[k], U, energy_str, filename_str );
 
    //return 0;


    // Test on the box
    for( int i = 0; i < N; i++ )
    {
        for( int j = 0; j < N; j++ )
        {  
            k = N*i + j;
    
            xL[0] = i * box_size / N;
            xL[1] = j * box_size / N;

            xU[0] = (i+1) * box_size / N;
            xU[1] = (j+1) * box_size / N;
            
            printf("MOF number %d: [%f, %f] x [%f, %f]\n", k, xL[0], xU[0], xL[1], xU[1]);
            U = Box( xL, xU, 2 );
            
            parabolic_MOF( refs[k], U, guess_type, error_tols, final_coeffs[k], max_iter, init_guesses[k] );
            this_symdiff_err = quad_symdiff( true_coeffs, final_coeffs[k] ); 
            
            sprintf(filename_str, "rbf_data2/MOF_%d.png", k);
            sprintf(energy_str, "Sym Diff: %.6f", this_symdiff_err);
            
            plot_level_sets_3( true_coeffs, final_coeffs[k], init_guesses[k], U, energy_str, filename_str );
        }
    }

    return 0;
}
