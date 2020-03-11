#include <stdio.h>
#include <cmath>
#include "saye_utils.h"
#include "saye_algorithm.h"
#include "grad_descent.h"
#include "rbf_network.h"
#include "level_set_plot.h"
#include <random>

using namespace std;

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

// Find the cost of a linear level set function using only centroid moment data
// Data is provided in parabola format, having an angle and center
// ref is length 5 to flush with quadratic construction
//  ref[5] = centroid x-coordinate
//  ref[4] = centroid y-coordinate
//  ref[3] = volume fraction
double cost_PLIC( double a[], double ref[], Box U )
{
    double the_cost = 0;

    double coeffs[6];
    poly_coefficients( coeffs, a );
    Phi2D phi0( coeffs );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{-1};

    // double xL[2] = {-1, -1};
    // double xU[2] = { 1,  1};
    // Box U( xL, xU, 2 );
    
    double temp;
    // double init_coeff[6];
    temp = ( I( psis, ss, U, &cx, false, 5 ) / ref[5] - ref[3] );
    the_cost += temp*temp;

    temp = ( I( psis, ss, U, &cy, false, 5 ) / ref[5] - ref[4] );
    the_cost += temp*temp;
    
    return the_cost;
}

/*
void gdescent_cost_fdm_PQIC( double an[], double anp1[], double gamma, double refs[], Box U )
{
    double grad_an[5];
    
    cost_fdm_PQIC( an, refs, U, grad_an, 0.00001, 0.0001 );

    for( int i = 0; i < 5; i++ )
        anp1[i] = an[i] - gamma*grad_an[i];
    
    flooding( anp1, refs[5], 0.0001 );
}
*/

/*
void gdescent_integral_fdm_PQIC( double an[], double anp1[], double gamma, double refs[] )
{
    double grad_an[5];
    
    integral_fdm_PQIC( an, refs, grad_an );

    for( int i = 0; i < 5; i++ )
        anp1[i] = an[i] - gamma*grad_an[i];
    
    flooding( anp1, refs[5], 0.0001 );
}
*/ 

/*
void gdescent_direct_PQIC( double an[], double anp1[], double gamma, double refs[] )
{
    double grad_an[5];
    
    direct_PQIC( an, refs, grad_an );

    for( int i = 0; i < 5; i++ )
        anp1[i] = an[i] - gamma*grad_an[i];
    
    flooding( anp1, refs[5], 0.0001 );
}
*/

// Do the linear MOF on the unit circle
void linear_MOF_unit( double refs[], double error_tols[], double final_params[], int max_iter )
{
    //double guess_tol = error_tols[0];
    double gradient_tol = error_tols[1];
    double flood_tol = error_tols[2];
    double steepest_tol = error_tols[0];

    double xL[2] = {0, 0};
    double xU[2] = {1, 1};
    Box U( xL, xU, 2 );

    double gamma = 0.001, gamma_num, gamma_den;

    vector<double> moment_errs;
    double this_moment_err;

    // Store coefficients
    double an[4] = {0, 1 - refs[5], 0, 0};

    double min_theta, min_cost = 100000;
    // Select initial condition kinda ad hoc
    for( an[3] = 0; an[3] < 8*atan(1); an[3] += 0.5 )
    {
        para_flooding( an, refs[5], U, flood_tol );
        if( cost_PLIC( an, refs, U ) < min_cost )
        {
            min_cost = cost_PLIC( an, refs, U );
            min_theta = an[3];
            printf( "Min cost: %f, Min theta: %f\n", min_cost, min_theta );
        }
    }
    
    // Start with the best guess
    an[3] = min_theta;
    para_flooding( an, refs[5], U, flood_tol );

    double anp1[4] = {0, 0, 0, 0};
    double anm1[4] = {0, 0, 0, 0};

    // Store gradients
    double gradn[4] = {0, 0, 0, 0};
    double gradnm1[4] = {0, 0, 0, 0};

    // Get initial error
    this_moment_err = cost_PLIC( an, refs, U );
    moment_errs.push_back( this_moment_err );
    
    printf( "L: 0 %.10f | %f %f %f %f | ---- ----\n", this_moment_err, an[0], an[1], an[2], an[3] );

    /*
    // Calculate first gradient
    cost_fdm_PLIC( an, refs, U, gradn, gradient_tol, flood_tol );  

    // Do the first steepest descent
    anp1[0] = 0;
    anp1[1] = an[1];
    anp1[2] = an[2];
    anp1[3] = an[3] - gamma*gradn[3];
    para_flooding( anp1, refs[5], U, flood_tol );
    
    this_moment_err = cost_PLIC( anp1, refs, U );
    moment_errs.push_back( this_moment_err );
    
    // Refresh the old information
    for( int j = 0; j < 4; j++ )
    {
        anm1[j] = an[j];
        an[j] = anp1[j];
        gradnm1[j] = gradn[j];
    }
    */
    double break_cond;
    double mypi = 4*atan(1);
    for( int i = 1; i < max_iter; i++ )
    {
        // Calculate gradient
        cost_fdm_PLIC( an, refs, U, gradn, gradient_tol, flood_tol );  

        // Calculate gamma
        gamma_num = gamma_den = 0;
        for( int j = 3; j < 4; j++ )
        {
            gamma_num += ( an[j] - anm1[j] ) * ( gradn[j] - gradnm1[j] );
            gamma_den += ( gradn[j] - gradnm1[j] ) * ( gradn[j] - gradnm1[j] );
        }
        gamma = 100.0;//dabs( gamma_num ) / gamma_den;
        
        // Backstepping parts
        //printf("L: %d Goal: %f, %f, %d\n", i, moment_errs[i-1], moment_errs.back(), moment_errs[i-1] == moment_errs.back() );
        int bs = 0;
        for( bs = 0; bs < 100; bs++ )
        {
            // Do the steepest descent
            anp1[0] = 0;
            anp1[1] = an[1];// - gamma*gradn[1];
            anp1[2] = an[2];// - gamma*gradn[2];
            anp1[3] = an[3] - gamma*gradn[3];
            para_flooding( anp1, refs[5], U, flood_tol );

            this_moment_err = cost_PLIC( anp1, refs, U ); 
            moment_errs.push_back( this_moment_err );
   
            //printf( "L: %d %.10f | %f %f %f %f | %f\n", i, this_moment_err, anp1[0], anp1[1], anp1[2], anp1[3], gamma );
   
            if( moment_errs[i] > moment_errs[i-1] )
            {
                moment_errs.pop_back();
                gamma = gamma * 0.9;
            }
            else
                break;
        }
        if( bs == 100 )
        {
            moment_errs.push_back( this_moment_err );
            //printf( "L: Can't backstep enough\n" );
        }
        printf( "L: %d %.10f | %f %f %f %f | %f %f\n", i, this_moment_err, anp1[0], anp1[1], anp1[2], anp1[3], gamma, gradn[3] );

        // Refresh the old information
        for( int j = 0; j < 4; j++ )
        {
            anm1[j] = an[j];
            an[j] = anp1[j];
            gradnm1[j] = gradn[j];
        }
        
        break_cond = 0;
        for( int j = 3; j < 4; j++ )
            break_cond += gradn[j]*gradn[j];
        if( sqrt( break_cond ) < steepest_tol )
            break;
    }

    for( int j = 0; j < 4; j++ )
        final_params[j] = anp1[j];
}

// refs are in U, we want to transform to [0,1], get the optimal coefficients in [0,1], then transform back to U
void parabolic_MOF( double refs[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[], int NN_data[] )
{
    double guess_tol = error_tols[0];
    double gradient_tol = error_tols[1];
    double flood_tol = error_tols[2];
    double steepest_tol = error_tols[3];

    double unit_refs[6];
    double unit_coeffs[6];
    double unit_init_coeffs[6];
    transform_refs_to_unit( refs, unit_refs, U );

    double xL[2] = {0, 0};
    double xU[2] = {1, 1};
    Box UnitU( xL, xU, 2 );

    // IN [0,1]^2 MODE

    double init_params[4] = {0, 1 - unit_refs[5], 0, 0};

    // Dont do any of this if the box is full or empty
    if( unit_refs[5] <= 0 )
    {
        poly_coefficients( unit_coeffs, init_params );
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[5] = 1;
        init_coeffs[5] = 1;
        return;
    }
    if( unit_refs[5] >= 1 )
    {
        poly_coefficients( unit_coeffs, init_params );
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[5] = -1;
        init_coeffs[5] = -1;
        return;
    }
    
    double coeff[6];

    double gamma = 0.001, gamma_num, gamma_den;

    vector<double> moment_errs;
    double this_moment_err;

    // Store coefficients
    double an[4];
    double anp1[4];
    double anm1[4];

    // Store gradients
    double gradn[4];
    double gradnm1[4];

    if( guess_type == 'l' )
    {
        double init_coeffs[6];
        linear_MOF_unit( unit_refs, error_tols, init_params, max_iter ); // get parameters for best line
        //line_to_param( unit_init_coeffs, init_params ); // turns coefficients into parameters?
        //scaling( unit_init_coeffs ); // scales the coefficients???
        poly_coefficients( unit_init_coeffs, init_params ); // gets coefficients from parameters
        init_params[3] = fmod( init_params[3], 8*atan(1) ); // shrinks angle
    }

    if( guess_type == 'p' )
    {
        int n_test_data = NN_data[0];
        int n_hidden_nodes = NN_data[1];

        // Save filenames so you dont need to recompute evrytime
        // char input_fn[50];
        // char output_fn[50]; 
        char centers_fn[50]; 
        char weights_fn[50];

        // sprintf( input_fn, "rbf_data2/input%d.csv", n_test_data);
        // sprintf( output_fn, "rbf_data2/output%d.csv", n_test_data);
        sprintf( centers_fn, "rbf_data2/centers%d_%d.csv", n_test_data, n_hidden_nodes);
        sprintf( weights_fn, "rbf_data2/weights%d_%d.csv", n_test_data, n_hidden_nodes);
       
        /*
        get_training_data_para( n_test_data, input_fn, output_fn ); 
        get_centers( n_hidden_nodes, input_fn, output_fn, centers_fn );
        get_weights( input_fn, output_fn, centers_fn, weights_fn ); 
        */

        eval_network( unit_refs, init_params, centers_fn, weights_fn );

        para_flooding( init_params, unit_refs[5], UnitU, flood_tol );
    }

    //poly_coefficients( coeff, init_params );
    //plot_2D_reference( coeff, coeff, UnitU, "yeet", "rbf_data2/MOF_iteration_0.png"); 
   
    poly_coefficients( unit_init_coeffs, init_params );
    transform_coeffs_to_nonunit( unit_init_coeffs, init_coeffs, U );
    for( int i = 0; i < 4; i++ )
        an[i] = init_params[i];

    this_moment_err = cost_PPIC( an, unit_refs, UnitU );
    moment_errs.push_back( this_moment_err );
    printf( "0 %.10f | %f %f %f %f | %f\n", this_moment_err, an[0], an[1], an[2], an[3], gamma );

    /*
    // Calculate first gradient
    cost_fdm_PPIC( an, unit_refs, UnitU, gradn, gradient_tol, flood_tol );  

    // Do the first steepest descent
    for( int j = 0; j < 4; j++ )
        anp1[j] = an[j] - gamma*gradn[j];
    para_flooding( anp1, unit_refs[5], UnitU, flood_tol );
    
    this_moment_err = cost_PPIC( anp1, unit_refs, UnitU );
    moment_errs.push_back( this_moment_err );
    printf( "1 %.10f | %f %f %f %f | %f\n", this_moment_err, anp1[0], anp1[1], anp1[2], anp1[3], gamma );
    
    // Refresh the old information
    for( int j = 0; j < 4; j++ )
    {
        anm1[j] = an[j];
        an[j] = anp1[j];
        gradnm1[j] = gradn[j];
    }
    */
    
    double break_cond;
    for( int i = 2; i < max_iter; i++ )
    {
    
        /*
        // Calculate gradient
        cost_fdm_PPIC( an, unit_refs, UnitU, gradn, gradient_tol, flood_tol );  

        // Calculate gamma
        gamma_num = gamma_den = 0;
        for( int j = 0; j < 4; j++ )
        {
            gamma_num += ( an[j] - anm1[j] ) * ( gradn[j] - gradnm1[j] );
            gamma_den += ( gradn[j] - gradnm1[j] ) * ( gradn[j] - gradnm1[j] );
        }
        gamma = 1.0;// dabs( gamma_num ) / gamma_den;
        */

        // Backstepping parts
        //printf("%d Goal: %f, %f, %d\n", i, moment_errs[i-1], moment_errs.back(), moment_errs[i-1] == moment_errs.back() );
        int bs = 0;
        gamma = 1.0;
        for( bs = 0; bs < 100; bs++ ) 
        {
            // Do the steepest descent
            for( int j = 0; j < 4; j++ )
                 anp1[j] = an[j];// - gamma*gradn[j];
            anp1[ (i-2) % 4 ] = an[ (i-2) % 4 ] 
            
            para_flooding( anp1, unit_refs[5], UnitU, flood_tol );
        
            this_moment_err = cost_PPIC( anp1, unit_refs, UnitU ); 
            moment_errs.push_back( this_moment_err );
            
            //printf( "%d %.10f | %f %f %f %f | %f\n", i, this_moment_err, an[0], an[1], an[2], an[3], gamma );

            if( moment_errs[i] > moment_errs[i-1] )
            {
                moment_errs.pop_back();
                gamma = gamma * 0.5;
            }
            else
                break;
        }
        if( bs == 100 )
        {
            moment_errs.push_back( this_moment_err );
            //printf("Can't backstep enough\n" );
        }
        printf( "%d %.10f | %f %f %f %f | %f\n", i, this_moment_err, an[0], an[1], an[2], an[3], gamma );
        //printf("\n");

        char energy_str[50];
        char filename_str[50];
        poly_coefficients( coeff, anp1 );

        sprintf(filename_str, "rbf_data2/MOF_iteration_%d.png", i);
        sprintf(energy_str, "Approximate Square Interface" ); //"Moment Error: %.6f", this_moment_err);
            

        //plot_2D_box_reference( coeff, 0, 0.5, 0, 0.5, energy_str, filename_str );
        

        // Refresh the old information
        for( int j = 0; j < 4; j++ )
        {
            anm1[j] = an[j];
            an[j] = anp1[j];
            gradnm1[j] = gradn[j];
        }
        
        //break_cond = 0;
        //for( int j = 0; j < 4; j++ )
        //    break_cond += gradn[j]*gradn[j];
        //if( sqrt( break_cond ) < steepest_tol )
        //    break;
    }
    
    poly_coefficients( coeff, anp1 );
    transform_coeffs_to_nonunit( coeff, final_coeffs, U );

    poly_coefficients( unit_init_coeffs, init_params );
    transform_coeffs_to_nonunit( unit_init_coeffs, init_coeffs, U );
}

// refs are in U, we want to transform to [0,1], get the optimal coefficients in [0,1], then transform back to U
void parabolic_SA_MOF( double refs[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[], int NN_data[] )
{

}



// refs are in U, we want to transform to [0,1], get the optimal coefficients in [0,1], then transform back to U
void quadratic_MOF( double refs[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[] )
{
    double guess_tol = error_tols[0];
    double gradient_tol = error_tols[1];
    double flood_tol = error_tols[2];
    double steepest_tol = error_tols[3];

    double unit_refs[6];
    double unit_coeffs[6];
    double unit_init_coeffs[6];
    transform_refs_to_unit( refs, unit_refs, U );

    double xL[2] = {0, 0};
    double xU[2] = {1, 1};
    Box UnitU( xL, xU, 2 );

    // IN [0,1]^2 MODE

    unit_init_coeffs[0] = unit_init_coeffs[1] = unit_init_coeffs[2] = 0;
    unit_init_coeffs[3] = -1;
    unit_init_coeffs[4] = 0;
    unit_init_coeffs[5] = 1 - unit_refs[5];

    // Dont do any of this if the box is full or empty
    if( unit_refs[5] <= 0 )
    {
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[5] = 1;
        init_coeffs[5] = 1;
        return;
    }
    if( unit_refs[5] >= 1 )
    {
        transform_coeffs_to_nonunit( unit_coeffs, final_coeffs, U );
        final_coeffs[5] = -1;
        init_coeffs[5] = -1;
        return;
    }
    
    double coeff[6];

    double gamma = 0.001, gamma_num, gamma_den;

    vector<double> moment_errs;
    double this_moment_err;

    // Store coefficients
    double an[6];
    double anp1[6];
    double anm1[6];

    // Store gradients
    double gradn[6];
    double gradnm1[6];

    if( guess_type == 'l' )
    {
        double init_params[4] = {0, 1 - unit_refs[5], 0, 0};
        linear_MOF_unit( unit_refs, error_tols, init_params, max_iter );
        poly_coefficients( unit_init_coeffs, init_params );
        scaling( unit_init_coeffs );
    }

    if( guess_type == 'p' )
    {
        int n_test_data = 1000;
        int n_hidden_nodes = 1000;

        // Save filenames so you dont need to recompute evrytime
        char input_fn[50];
        char output_fn[50]; 
        char centers_fn[50]; 
        char weights_fn[50];

        sprintf( input_fn, "rbf_data2/quad_input%d.csv", n_test_data);
        sprintf( output_fn, "rbf_data2/quad_output%d.csv", n_test_data);
        sprintf( centers_fn, "rbf_data2/quad_centers%d_%d.csv", n_test_data, n_hidden_nodes);
        sprintf( weights_fn, "rbf_data2/quad_weights%d_%d.csv", n_test_data, n_hidden_nodes);
        
        
        get_training_data_quad( n_test_data, input_fn, output_fn ); 
        get_centers( n_hidden_nodes, input_fn, output_fn, centers_fn );
        get_weights( input_fn, output_fn, centers_fn, weights_fn ); 
        

        eval_network( unit_refs, unit_init_coeffs, centers_fn, weights_fn );

        flooding( unit_init_coeffs, unit_refs[5], UnitU, flood_tol );
    }

    //poly_coefficients( coeff, init_params );
    //plot_2D_reference( coeff, coeff, UnitU, "yeet", "rbf_data2/MOF_iteration_0.png"); 
   
    for( int i = 0; i < 4; i++ )
        an[i] = unit_init_coeffs[i];

    printf( "0 %f %f | %f %f %f %f %f %f\n", this_moment_err, gamma, an[0], an[1], an[2], an[3], an[4], an[5] );
    this_moment_err = cost_PQIC( an, unit_refs, UnitU );
    moment_errs.push_back( this_moment_err );

    // Calculate first gradient
    cost_fdm_PQIC( an, unit_refs, UnitU, gradn, gradient_tol, flood_tol );  

    // Do the first steepest descent
    for( int j = 0; j < 6; j++ )
        anp1[j] = an[j] - gamma*gradn[j];
    flooding( anp1, unit_refs[5], UnitU, flood_tol );
    
    printf( "1 %f %f | %f %f %f %f %f %f\n", this_moment_err, gamma, anp1[0], anp1[1], anp1[2], anp1[3], anp1[4], anp1[5] );
    
    this_moment_err = cost_PQIC( anp1, unit_refs, UnitU );
    moment_errs.push_back( this_moment_err );
    
    // Refresh the old information
    for( int j = 0; j < 6; j++ )
    {
        anm1[j] = an[j];
        an[j] = anp1[j];
        gradnm1[j] = gradn[j];
    }
    
    double break_cond;
    for( int i = 2; i < max_iter; i++ )
    {

        // Calculate gradient
        cost_fdm_PQIC( an, unit_refs, UnitU, gradn, gradient_tol, flood_tol );  

        // Calculate gamma
        gamma_num = gamma_den = 0;
        for( int j = 0; j < 6; j++ )
        {
            gamma_num += ( an[j] - anm1[j] ) * ( gradn[j] - gradnm1[j] );
            gamma_den += ( gradn[j] - gradnm1[j] ) * ( gradn[j] - gradnm1[j] );
        }
        gamma = dabs( gamma_num ) / gamma_den;
        printf( "%d %f %f | %f %f %f %f %f %f \n", i, this_moment_err, gamma, an[0], an[1], an[2], an[3], an[4], an[5]);

        // Do the steepest descent
        for( int j = 0; j < 6; j++ )
            anp1[j] = an[j] - gamma*gradn[j];
        flooding( anp1, unit_refs[5], UnitU, flood_tol );

        char energy_str[50];
        char filename_str[50];

        sprintf(filename_str, "rbf_data2/MOF_iteration_%d.png", i);
        sprintf(energy_str, "Approximate Square Interface" ); //"Moment Error: %.6f", this_moment_err);
            

        //plot_2D_box_reference( coeff, 0, 0.5, 0, 0.5, energy_str, filename_str );
        
        this_moment_err = cost_PQIC( anp1, unit_refs, UnitU ); 
        moment_errs.push_back( this_moment_err );

        // Refresh the old information
        for( int j = 0; j < 4; j++ )
        {
            anm1[j] = an[j];
            an[j] = anp1[j];
            gradnm1[j] = gradn[j];
        }
        
        break_cond = 0;
        for( int j = 0; j < 6; j++ )
            break_cond += gradn[j]*gradn[j];
        if( sqrt( break_cond ) < steepest_tol )
            break;
    }
    
    transform_coeffs_to_nonunit( anp1, final_coeffs, U );
    transform_coeffs_to_nonunit( unit_init_coeffs, init_coeffs, U );
}



// Perform one gradient descent step for PLIC
//  an: current level set coefficients
//  anp1: next level set coefficients
//  refs: reference data (6 long vector, but only care about indexes 3 4 5
/*
void gdescent_cost_fdm_PLIC( double an[], double anp1[], double gamma, double refs[] )
{
    double grad_an[5];
    
    cost_fdm_PLIC( an, refs, grad_an, 0.00001,0.0001 );

    anp1[0] = anp1[1] = anp1[2] = 0.0;
    anp1[3] = an[3] - gamma*grad_an[3];
    anp1[4] = an[4] - gamma*grad_an[4];
    
    flooding( anp1, refs[5], 0.0001 );
}
*/

// Find cost gradient
void cost_fdm_PLIC( double cn[], double refs[], Box U, double grad[], double h, double flood_tol )
{
    double original;
    double an[6];
    int i = 3;
        
    grad[i] = 0;
    original = cn[i];

    cn[i] = original + 2*h;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol);
    grad[i] += -cost_PQIC( an, refs, U );

    cn[i] = original + h;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol );
    grad[i] += 8*cost_PQIC( an, refs, U );
    
    cn[i] = original - h;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol );
    grad[i] += -8*cost_PQIC( an, refs, U );
    
    cn[i] = original - 2*h;
    poly_coefficients( an, cn );
    flooding( an, refs[5], U, flood_tol );
    grad[i] += cost_PQIC( an, refs, U );

    grad[i] /= 12*h;
    cn[i] = original;
}

// cost function of PQIC 
//  a (R^6) is coefficients for phi
//  a0 (R^6) is reference data (a0[5] is vof)
double cost_PQIC( double a[], double a0[], Box U )
{
    init_bases();
    double the_cost = 0;

    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    //double xL[2] = {-1, -1};
    //double xU[2] = { 1,  1};
    //Box U( xL, xU, 2 );
    int q = 5;
    
    for( int i = 0; i < 5; i++ )
    {
        double temp = I( psis, ss, U, bases[i], false, q ) / a0[5];
        //printf("\ttemp: %phi2d phi0( a );f\n", temp);
        the_cost += ( temp - a0[i] )*( temp - a0[i] );
    }

    return the_cost;
}

double cost_PPIC( double p[], double ref[], Box U )
{
    double a[6];
    poly_coefficients( a, p );

    return cost_PQIC( a, ref, U );
}


void direct_PQIC( double an[], double refs[], double grad[] )
{
    Phi2D phi0( an );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss_vol{ -1 };
    vector<int> ss_surf{ 0 };

    double xL[2] = {-1, -1};
    double xU[2] = { 1,  1};
    Box U( xL, xU, 2 );
    int q = 5;
    
    init_bases();
    init_grad_bases( &phi0 );

    double temp;

    for( int j = 0; j < 5; j++ )
    {
        grad[j] = 0;

        for( int i = 0; i < 5; i++ )
        {
            temp = I( psis, ss_vol, U, bases[i], false, q ) / refs[5] - refs[i];
            temp = temp * I( psis, ss_surf, U, d_bases[i][j], true, q );
        }

        grad[j] *= 2 / refs[5];
    }
}

void integral_fdm_PQIC( double an[], double refs[], Box U, double grad[] )
{
    Phi2D phi0( an );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss_vol{ -1 };
    vector<int> ss_surf{ 0 };

    // double xL[2] = {-1, -1};
    // double xU[2] = { 1,  1};
    // Box U( xL, xU, 2 );
    int q = 5;
    init_bases();
    
    double temp_error, temp_fdm, original;
    double h = 0.00001, tol = 0.0001;

    for( int j = 0; j < 5; j++ )
    {
        grad[j] = 0;

        for( int i = 0; i < 5; i++ )
        {
            temp_error = I( psis, ss_vol, U, bases[i], false, q ) / refs[5] - refs[i];
            temp_fdm = 0;

            original = an[j];

            an[j] = original + 2*h;
            flooding( an, refs[5], U, tol );
            phi0 = Phi2D( an );
            psis[0] = Psi( &phi0 );
            temp_fdm += -1*I( psis, ss_vol, U, bases[i], false, q ) / refs[5];
            
            an[j] = original + h;
            flooding( an, refs[5], U, tol );
            phi0 = Phi2D( an );
            psis[0] = Psi( &phi0 );
            temp_fdm += 8*I( psis, ss_vol, U, bases[i], false, q ) / refs[5];
            
            an[j] = original - h;
            flooding( an, refs[5], U, tol );
            phi0 = Phi2D( an );
            psis[0] = Psi( &phi0 );
            temp_fdm += -8*I( psis, ss_vol, U, bases[i], false, q ) / refs[5];

            an[j] = original - 2*h;
            flooding( an, refs[5], U, tol );
            phi0 = Phi2D( an );
            psis[0] = Psi( &phi0 );
            temp_fdm += I( psis, ss_vol, U, bases[i], false, q ) / refs[5];

            temp_fdm = temp_fdm / 12 / h;

            an[j] = original;

            grad[j] += temp_error * temp_fdm;
        }

        grad[j] *= 2 / refs[5];
    }
}

// Take partial with respect to ith value in an
void cost_fdm_PQIC( double an[], double refs[], Box U, double grad[], double h, double flood_tol )
{
    double original;

    for( int i = 0; i < 5; i++ )
    {
        grad[i] = 0;
        original = an[i];

        an[i] = original + 2*h; 
        flooding( an, refs[5], U, flood_tol );
        grad[i] += -cost_PQIC( an, refs, U );

        an[i] = original + h;
        flooding( an, refs[5], U, flood_tol );
        grad[i] += 8*cost_PQIC( an, refs, U );
        
        an[i] = original - h;
        flooding( an, refs[5], U, flood_tol );
        grad[i] += -8*cost_PQIC( an, refs, U );
        
        an[i] = original + 2*h;
        flooding( an, refs[5], U, flood_tol );
        grad[i] += cost_PQIC( an, refs, U );

        grad[i] /= 12*h;
        an[i] = original;
    }

    grad[5] = 0;
}

void cost_fdm_PPIC( double cn[], double refs[], Box U, double grad[], double h, double flood_tol )
{
    double original;
    double an[6];

    for( int i = 0; i < 4; i++ )
    {
        grad[i] = 0;
        original = cn[i];

        cn[i] = original + 2*h;
        poly_coefficients( an, cn );
        flooding( an, refs[5], U, flood_tol);
        grad[i] += -cost_PQIC( an, refs, U );

        cn[i] = original + h;
        poly_coefficients( an, cn );
        flooding( an, refs[5], U, flood_tol );
        grad[i] += 8*cost_PQIC( an, refs, U );
        
        cn[i] = original - h;
        poly_coefficients( an, cn );
        flooding( an, refs[5], U, flood_tol );
        grad[i] += -8*cost_PQIC( an, refs, U );
        
        cn[i] = original + 2*h;
        poly_coefficients( an, cn );
        flooding( an, refs[5], U, flood_tol );
        grad[i] += cost_PQIC( an, refs, U );

        grad[i] /= 12*h;
        cn[i] = original;
    }
}

void visualize_MOF( double refs[], Box U, double flood_tol, const char* filename )
{
    double dt = 0.01;
    double theta;
    vector<double> costs;
    vector<double> thetas;


    for( theta = 0; theta < 8*atan(1); theta += dt )
    {
        double a[] = {0, 0.5, 0.5, theta};
        para_flooding( a, refs[5], U, flood_tol );
        costs.push_back( cost_PLIC( a, refs, U ) );
        thetas.push_back( theta );
    }

    plot_linear_MOF( costs, thetas, filename );
}

void para_flooding( double c[], double vof, Box U, double tol )
{
    double a[6];

    // Get terms in coefficient terms
    poly_coefficients( a, c );
    double new_const = a[5];
    
    // Do flooding on the coefficients
    flooding( a, vof, U, tol );
    double const_diff = a[5] - new_const;

    // Adjust parameters to match coefficient change
    c[1] += const_diff * cos( c[3] );
    c[2] += const_diff * sin( c[3] );
}

void flooding( double a[], double vof, Box U, double tol )
{
    vector<int> ss{ -1 };
    //double xL[2] = {0, 0};
    //double xU[2] = { 1,  1};
    //U = Box( xL, xU, 2 );
    int q = 5;
    Unit unit;

    Phi2D phi0( a );
    vector<Psi> psis{ Psi( &phi0 ) };

    double low, high, end;
    double start = I( psis, ss, U, &unit, false, q );
    // printf("%f\n", start);

    double test_vol;
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
            
            if( end >= vof )
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

            if( end <= vof )
                break;
            else
                high += 1;
        }
    }

    // Bisection method part
    test_vol = end;
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

}

void get_refs( Box U, double coeffs[], double refs[] )
{
    Phi2D phi0( coeffs );
    vector<Psi> psis{ Psi( &phi0 ) };
    vector<int> ss{ -1 };

    bool S = false;
    int q = 5;

    Unit unit;
    refs[5] = I( psis, ss, U, &unit, S, q );

    if( refs[5] <= 0 || refs[5] >= U.volume() )
    {
        refs[0] = refs[1] = refs[2] = refs[3] = refs[4];
        return;
    }

    Mx2 mx2;
    refs[0] = I( psis, ss, U, &mx2, S, q ) / refs[5]; 
    
    My2 my2;
    refs[1] = I( psis, ss, U, &my2, S, q ) / refs[5]; 
    
    Mxy mxy;
    refs[2] = I( psis, ss, U, &mxy, S, q ) / refs[5]; 
    
    Cx cx;
    refs[3] = I( psis, ss, U, &cx, S, q ) / refs[5]; 
    
    Cy cy;
    refs[4] = I( psis, ss, U, &cy, S, q ) / refs[5]; 
}

void transform_refs_to_unit( double refs[], double refs_t[], Box U )
{
    double a = U[0][0];
    double b = U[0][1];
    double c = U[1][0];
    double d = U[1][1];

    // vof
    refs_t[5] = refs[5] / (b-a) / (d-c);

    // Cx
    refs_t[3] = ( refs[3] - a ) / ( b - a );

    // Cy
    refs_t[4] = ( refs[4] - c ) / ( d - c );

    // Mxx
    refs_t[0] = ( refs[0] - 2*a*refs[3] + a*a ) / (b-a) / (b-a);
    
    // Myy
    refs_t[1] = ( refs[1] - 2*c*refs[4] + c*c ) / (d-c) / (d-c);

    // Mxy
    refs_t[2] = ( refs[2] - c*refs[3] - a*refs[4] + a*c ) / (b-a) / (d-c);
}

// arg[0] is original coefficients on unit, arg[1] is coefficients on nonunit
void transform_coeffs_to_nonunit( double coeffs[], double coeffs_t[], Box U )
{
    double a = U[0][0];
    double b = U[0][1];
    double c = U[1][0];
    double d = U[1][1];

    coeffs_t[0] = coeffs[0] / (b-a) / (b-a);

    coeffs_t[1] = coeffs[1] / (d-c) / (d-c);

    coeffs_t[2] = coeffs[2] / (b-a) / (d-c);

    coeffs_t[3] = -2*coeffs[0]*a / (b-a)/(b-a) - coeffs[2]*c / (b-a)/(d-c) + coeffs[3] / (b-a);

    coeffs_t[4] = -2*coeffs[1]*c / (d-c)/(d-c) - coeffs[2]*a / (b-a)/(d-c) + coeffs[4] / (d-c);

    coeffs_t[5] = coeffs[0]*a*a/(b-a)/(b-a) + coeffs[1]*c*c/(d-c)/(d-c) + 
                  coeffs[2]*a*c/(b-a)/(d-c) - coeffs[3]*a/(b-a) - coeffs[4]*c/(d-c) + coeffs[5];

}

double box_symdiff( int N, double xL, double xU, double yL, double yU, double a[] )
{
    default_random_engine generator;
    uniform_real_distribution<double> distribution(-1.0, 1.0);
    double quad_val;
    int box_val;
    
    int inside_count = 0;

    double x, y;
    for( int i = 0; i < N; i++ )
    {
        x = distribution( generator );
        y = distribution( generator );

        quad_val = a[0]*x*x + a[1]*y*y + a[2]*x*y + a[3]*x + a[4]*y + a[5];
        
        if( xL < x && x < xU && yL < y && y < yU )
            box_val = -1;
        else
            box_val = 1;
    
        if( box_val * quad_val <= 0 )
            inside_count += 1;
    }

    return 4.0 * inside_count / N;
}


double quad_symdiff( double a[], double b[] )
{
    Phi2D phia( a );
    Phi2D phib( b );

    vector<Psi> psis{ Psi( &phia ), Psi( &phib ) };
    vector<int> ssa{ -1, 1 };
    vector<int> ssb{ 1, -1 };

    double xL[2] = {-1, -1};
    double xU[2] = { 1,  1};
    Box U( xL, xU, 2 );

    bool S = false;
    int q = 5;

    Unit vof;

    return I( psis, ssa, U, &vof, S, q ) + I( psis, ssb, U, &vof, S, q );
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
{
    gx4.phi0 = phi0;
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
