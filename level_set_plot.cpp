#include <vector>
#include "saye_utils.h"
#include "matplotlibcpp.h"
#include "level_set_plot.h"
namespace plt = matplotlibcpp;

void plot_2D( double a[], double ref_x, double ref_y, const char * filename )
{
    plt::plot( {-1, 1}, {-1,-1}, "k-" );
    plt::plot( {-1, 1}, { 1, 1}, "k-" );
    plt::plot( {-1,-1}, {-1, 1}, "k-" );
    plt::plot( { 1, 1}, {-1, 1}, "k-" );
    plt::axis( "equal" );

    int N = 1000;
    double tol = 0.001;
    double x, y;

    for( int i = 0; i < N; i++ )
        for( int j = 0; j < N; j++ )
        {
            x = 2.0 * i / ( N - 1 ) - 1; 
            y = 2.0 * j / ( N - 1 ) - 1; 
    
            if( dabs( a[0]*x*x + a[1]*y*y + a[2]*x*y + a[3]*x + a[4]*y + a[5] ) < tol )
                plt::plot( {x}, {y}, "b." );
        }

    plt::plot( {ref_x}, {ref_y}, "r*" );
    plt::save( filename );
    plt::clf();
}

void plot_2D_reference( double a[], double b[], Box U, const char* energy, const char * filename )
{
    plt::plot( {U[0][0], U[0][1]}, {U[1][0], U[1][0]}, "k-" );
    plt::plot( {U[0][0], U[0][1]}, {U[1][1], U[1][1]}, "k-" );
    plt::plot( {U[0][0], U[0][0]}, {U[1][0], U[1][1]}, "k-" );
    plt::plot( {U[0][1], U[0][1]}, {U[1][0], U[1][1]}, "k-" );
    plt::axis( "equal" );

    int N = 500;
    double tol = 0.001;
    double x, y;
    double dx = ( U[0][1] - U[0][0] ) / N;
    double dy = ( U[1][1] - U[1][0] ) / N;

    // Plot edges
    for( x = U[0][0]; x < U[0][1]; x += dx )
        for( y = U[1][0]; y < U[1][1]; y += dy )
        {
            if( dabs( a[0]*x*x + a[1]*y*y + a[2]*x*y + a[3]*x + a[4]*y + a[5] ) < tol )
                plt::plot( {x}, {y}, "b." );

            if( dabs( b[0]*x*x + b[1]*y*y + b[2]*x*y + b[3]*x + b[4]*y + b[5] ) < tol )
                plt::plot( {x}, {y}, "r." );
        }
    
    for( x = U[0][0]; x < U[0][1]; x += 100*dx )
        for( y = U[1][0]; y < U[1][1]; y += 100*dy )
        {
            if( a[0]*x*x + a[1]*y*y + a[2]*x*y + a[3]*x + a[4]*y + a[5] < 0 )
                plt::plot( {x}, {y}, "bo" );

            if( b[0]*x*x + b[1]*y*y + b[2]*x*y + b[3]*x + b[4]*y + b[5] < 0 )
                plt::plot( {x}, {y}, "r." );
        }



    plt::title( energy );
    plt::save( filename );
    plt::clf();

}

void plot_level_sets_3( double a[], double b[], double c[], Box U, const char* plot_title, const char * filename )
{
    plt::plot( {U[0][0], U[0][1]}, {U[1][0], U[1][0]}, "k-" );
    plt::plot( {U[0][0], U[0][1]}, {U[1][1], U[1][1]}, "k-" );
    plt::plot( {U[0][0], U[0][0]}, {U[1][0], U[1][1]}, "k-" );
    plt::plot( {U[0][1], U[0][1]}, {U[1][0], U[1][1]}, "k-" );
    plt::axis( "equal" );

    int N = 500;
    double tol = 0.001;
    double x, y;
    double dx = ( U[0][1] - U[0][0] ) / N;
    double dy = ( U[1][1] - U[1][0] ) / N;

    // Plot edges
    for( x = U[0][0]; x < U[0][1]; x += dx )
        for( y = U[1][0]; y < U[1][1]; y += dy )
        {
            if( dabs( a[0]*x*x + a[1]*y*y + a[2]*x*y + a[3]*x + a[4]*y + a[5] ) < tol )
                plt::plot( {x}, {y}, "b." );

            if( dabs( b[0]*x*x + b[1]*y*y + b[2]*x*y + b[3]*x + b[4]*y + b[5] ) < tol )
                plt::plot( {x}, {y}, "r." );
            
            if( dabs( c[0]*x*x + c[1]*y*y + c[2]*x*y + c[3]*x + c[4]*y + c[5] ) < tol )
                plt::plot( {x}, {y}, "g." );
        }
    
    for( x = U[0][0]; x < U[0][1]; x += 100*dx )
        for( y = U[1][0]; y < U[1][1]; y += 100*dy )
        {
            if( c[0]*x*x + c[1]*y*y + c[2]*x*y + c[3]*x + c[4]*y + c[5] < 0 )
                plt::plot( {x}, {y}, "g*" );
            
            if( a[0]*x*x + a[1]*y*y + a[2]*x*y + a[3]*x + a[4]*y + a[5] < 0 )
                plt::plot( {x}, {y}, "bo" );

            if( b[0]*x*x + b[1]*y*y + b[2]*x*y + b[3]*x + b[4]*y + b[5] < 0 )
                plt::plot( {x}, {y}, "r." );
        }



    plt::title( plot_title );
    plt::save( filename );
    plt::clf();
}



void plot_2D_box_reference( double a[], double xL, double xU, double yL, double yU, 
                            const char* plot_title, const char* filename )
{
    plt::plot( {xL, xU}, {yL, yL}, "r-" );
    plt::plot( {xL, xU}, {yU, yU}, "r-" );
    plt::plot( {xL, xL}, {yL, yU}, "r-" );
    plt::plot( {xU, xU}, {yL, yU}, "r-" );
    
    plt::plot( {0, 1}, {0, 0}, "k-" );
    plt::plot( {0, 1}, {1, 1}, "k-" );
    plt::plot( {0, 0}, {0, 1}, "k-" );
    plt::plot( {1, 1}, {0, 1}, "k-" );

    plt::axis( "equal" );

    int N = 1000;
    double tol = 0.001;
    double x, y;
    double dx = 1.0 / N;
    double dy = 1.0 / N;
    
    // Plot edges
    for( x = 0; x < 1; x += dx )
        for( y = 0; y < 1; y += dy )
        {
            if( dabs( a[0]*x*x + a[1]*y*y + a[2]*x*y + a[3]*x + a[4]*y + a[5] ) < tol )
                plt::plot( {x}, {y}, "b." );
        }

    plt::title( plot_title );
    plt::save( filename );
    plt::clf();
}

//void bisection_plot( 
void plot_histogram( std::vector<double> a, const char* fig_title, const char* filename )
{
    plt::hist( a, 50 );
    plt::title( fig_title );
    plt::save( filename );
    plt::clf();
}
