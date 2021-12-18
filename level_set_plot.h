#ifndef LEVELPLOT
#define LEVELPLOT

#include <vector> 

void plot_2D( double a[], double ref_x, double ref_y, const char * filename );
void plot_2D_reference( double a[], double b[], Box U, const char * energy, const char * filename );
void plot_level_sets_3( double a[], double b[], double c[], Box U, const char * plot_title, const char * filename );
void plot_2D_box_reference( double a[], double xL, double xU, double yL, double yU, 
                            const char* plot_title, const char* filename );
void plot_histogram( std::vector<double> a, const char * fig_title, const char * filename );

#endif
