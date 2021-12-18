#include <cmath>
#include <iostream>
#include <armadillo>

#include "saye_utils.h"
#include "saye_algorithm.h"
#include "rbf_network.h" 
#include "level_set_plot.h"

using namespace std;
using namespace arma;

// We want this to train the NN on the domain [0,1]^2, but until those transformations
//  are discovered, we instead just take in U as a parameter and do it the slow way
void get_training_data_quad( int N, const char input_NN_fn[], const char output_NN_fn[] )
{
    arma_rng::set_seed_random();

    mat Y = randu<mat>( 6, N );
    mat X = mat( 6, N );

    Mx2 mx2;
    My2 my2;
    Mxy mxy;
    Cx cx;
    Cy cy;
    Unit unit;

    vector<int> ss{-1};
    double xL[2] = { 0, 0};
    double xU[2] = { 1, 1};
    Box U( xL, xU, 2 );
    int q = 5;

    // double param[4];
    double coeff[6];
    double vof;
    

    for( int i = 0; i < N; i++ )
    {
        //printf( "creating data point %d...\n", i );

        coeff[0] = Y(0, i);
        coeff[1] = Y(1, i);
        coeff[2] = Y(2, i);
        coeff[3] = Y(3, i);
        coeff[4] = Y(4, i);
        coeff[5] = Y(5, i);

        Phi2D phi0( coeff );
        vector<Psi> psis{ Psi( &phi0 ) };
    
        X(5, i) = vof = I( psis, ss, U, &unit, false, q );

        if( vof <= 0 || vof >= 1 )
        {
            Y.col(i) = randu<mat>( 6, 1 );
            
            i--;
            continue;
        }
        
        X(0, i) = I( psis, ss, U, &mx2, false, q ) / vof; 
        X(1, i) = I( psis, ss, U, &my2, false, q ) / vof; 
        X(2, i) = I( psis, ss, U, &mxy, false, q ) / vof; 
        X(3, i) = I( psis, ss, U, &cx, false, q ) / vof; 
        X(4, i) = I( psis, ss, U, &cy, false, q ) / vof; 
    
    }

    Y.save( output_NN_fn, csv_ascii );
    X.save( input_NN_fn, csv_ascii );

    //Y.save( "output.txt", arma_ascii );
    //X.save( "input.txt", arma_ascii );
}

// We want this to train the NN on the domain [0,1]^2, but until those transformations
//  are discovered, we instead just take in U as a parameter and do it the slow way
void get_training_data_para( int N, const char input_NN_fn[], const char output_NN_fn[] )
{
    arma_rng::set_seed_random();

    // mat Y = 2 * randu<mat>( 6, N ) - 1;
    mat Y = randn<mat>( 4, N );

    //Y.row(0) = Y.row(0) * 3;
    Y.row(0) = abs( Y.row(0) * 3 );
    Y.row(1) = Y.row(1) / 2 + 0.5;
    Y.row(2) = Y.row(2) / 2 + 0.5;
    Y.row(3) = 2 * datum::pi * randu<mat>(1, N);

    mat X = mat( 6, N );

    Mx2 mx2;
    My2 my2;
    Mxy mxy;
    Cx cx;
    Cy cy;
    Unit unit;

    vector<int> ss{-1};
    double xL[2] = { 0, 0};
    double xU[2] = { 1, 1};
    Box U( xL, xU, 2 );
    int q = 5;

    double param[4];
    double coeff[6];
    double vof;
    
    char filename_str[30];

    for( int i = 0; i < N; i++ )
    {
        //printf( "creating data point %d...\n", i );

        param[0] = Y(0, i);
        param[1] = Y(1, i);
        param[2] = Y(2, i);
        param[3] = Y(3, i);

        poly_coefficients( coeff, param );

        /*
        coeff[0] = Y(0, i);
        coeff[1] = Y(1, i);
        coeff[2] = Y(2, i);
        coeff[3] = Y(3, i);
        coeff[4] = Y(4, i);
        coeff[5] = Y(5, i);
        */

        Phi2D phi0( coeff );
        vector<Psi> psis{ Psi( &phi0 ) };
    
        X(5, i) = vof = I( psis, ss, U, &unit, false, q );

        if( vof <= 0 || vof >= 1 )
        {
            // mat Y = 2 * randu<mat>( 6, N ) - 1;
            Y.col(i) = randn<mat>( 4, 1 );

            //Y(0, i) = Y(0, i) * 3;
            Y(0, i) = abs( Y(0, i) * 3 );
            Y(1, i) = Y(1, i) / 2 + 0.5;
            Y(2, i) = Y(2, i) / 2 + 0.5;
            Y(3, i) = 2 * datum::pi * randu();
            
            i--;
            continue;
        }
        
        X(0, i) = I( psis, ss, U, &mx2, false, q ) / vof; 
        X(1, i) = I( psis, ss, U, &my2, false, q ) / vof; 
        X(2, i) = I( psis, ss, U, &mxy, false, q ) / vof; 
        X(3, i) = I( psis, ss, U, &cx, false, q ) / vof; 
        X(4, i) = I( psis, ss, U, &cy, false, q ) / vof; 
        
        //sprintf( filename_str, "rbf_data2/test_data%d.png", i );
        //plot_2D_reference( coeff, coeff, U, "yeet", filename_str );
    }

    Y.save( output_NN_fn, csv_ascii );
    X.save( input_NN_fn, csv_ascii );

    //Y.save( "output.txt", arma_ascii );
    //X.save( "input.txt", arma_ascii );
}

// J2 is the number of nodes in the NN
void get_centers( int J2, const char input_NN_fn[], const char output_NN_fn[], const char centers_fn[] )
{
    mat X;
    X.load( input_NN_fn );

    int N = X.n_cols;
    
    // select initial generators 
    mat C = X.cols( 0, J2 - 1 );
    ivec bins = randi( N, distr_param(0, J2-1) );

    //cout << X.n_rows << " " << X.n_cols << endl;
    //cout << C.n_rows << " " << C.n_cols << endl;
    
    double energy = 0;
    // print k-means energy
    for( int i = 0; i < N; i++ )
        energy += norm( X.col(i) - C.col( bins(i) ) ) * norm( X.col(i) - C.col( bins(i) ) );
    //printf("Energy 0\t%f\n", energy );
    
    for( int k = 0; k < 3; k++ )
    {
        // Find the closest generator for each point
        for( int i = 0; i < N; i++ )
        {
            double test_dist;
            double min_dist = norm( X.col(i) - C.col(0) ); // Initialize to distance to first cluster
            bins(i) = 0;
            
            for( int j = 1; j < J2; j++ )
            {
                test_dist = norm( X.col(i) - C.col(j) );
                if( test_dist*test_dist < min_dist*min_dist )
                {
                    min_dist = test_dist;
                    bins(i) = j;
                }
            }
        }

        //cout << sum( X.cols( find( bins == 0 ) ), 1 ) << endl;
        //cout << sum( X.cols( find( bins == 0 ) ) ) << endl;
        // Set each cluster to the average of the closest generators
        //cout << X.cols( find( bins == j ) ) << endl;
        //cout << bins.elem( find( bins == j
        for( int j = 0; j < J2; j++ )
        {
            C.col(j) = sum( X.cols( find( bins == j ) ), 1 );
            C.col(j) = C.col(j) / size( X.cols( find( bins == j ) ), 1 );
        }
   
        energy = 0;
        // print k-means energy
        for( int i = 0; i < N; i++ )
            energy += norm( X.col(i) - C.col( bins(i) ) ) * norm( X.col(i) - C.col( bins(i) ) );
        //printf("Energy %d\t%f\n", k+1, energy );
    
    }

    C.save( centers_fn );
}

// J2 is the number of nodes in the NN
void get_weights( const char input_NN_fn[], const char output_NN_fn[], const char centers_fn[], const char weights_fn[] )
{
    mat X, Y, C;

    X.load( input_NN_fn );
    Y.load( output_NN_fn );
    C.load( centers_fn );

    int N = Y.n_cols;
    int J2 = C.n_cols;
    // int J3 = Y.n_rows;
    double sigma = 1;
    
    cout << N << endl;
    cout << J2 << endl;

    //mat W = mat( J2, 6 );
    mat W = mat( J2, J2 );
    mat P = mat( J2, N );

    // set up Phi matrix
    for( int p = 0; p < N; p++ )
        for( int k = 0; k < J2; k++ )
            P(k, p) = RBF( norm( X.col(p) - C.col(k) ), sigma );

    
    mat U, V;
    vec s;

    svd( U, s, V, P );

    mat S = diagmat(s);
    S.resize( J2, N );

    W = inv( P * P.t() ) * P * Y.t();
    //W = U * inv( S*S.t() ) * S * V.t() * Y.t();
    
    W.save( weights_fn );
}

// We want this function to be able to translate refs[] defined on [a,b]x[c,d]
// into a valid shape on [0,1]^2, get the coefficients fron the NN, 
// then translate those back into [a,b]x[c,d]
void eval_network( double refs[], double coeffs[], const char centers_fn[], const char weights_fn[] )
{
    vec x( 6 );
    
    x(0) = refs[0];
    x(1) = refs[1];
    x(2) = refs[2];
    x(3) = refs[3];
    x(4) = refs[4];
    x(5) = refs[5];

    mat C, W;

    C.load( centers_fn );
    W.load( weights_fn );

    int J2 = W.n_rows;
    int J3 = W.n_cols;
    double sigma = 1;

    //for( int i = 0; i < 6; i++ )
    for( int i = 0; i < J3; i++ )
        for( int k = 0; k < J2; k++ )
            coeffs[i] += W(k, i) * RBF( norm( x - C.col(k) ), sigma );
}

double RBF( double r, double sigma )
{
    return exp( -r*r/2/sigma/sigma );
}

/*
void poly_coefficients( double a[], double c[] ) 
{
    double sinth = sin(c[3]);
    double costh = cos(c[3]);
    double cosmsin = c[2]*costh - c[1]*sinth;
    
    a[0] =  c[0]*sinth*sinth/4.0;
    a[1] =  c[0]*costh*costh/4.0;
    a[2] = -c[0]*costh*sinth/2.0;
    a[3] =  c[0]/2.0*sinth*cosmsin - costh;
    a[4] = -c[0]/2.0*costh*cosmsin - sinth;
    a[5] =  c[1]*costh + c[2]*sinth + c[0]/4.0*cosmsin*cosmsin;
}
*/
// CHeck symetric difference error
