#ifndef GRADDESC
#define GRADDESC 

#include "saye_utils.h"

double cost_PLIC( double a[], double ref[], Box U ); //
double cost_PQIC( double a[], double a0[], Box U ); //
double cost_PPIC( double p[], double ref[], Box U ); //

//void gdescent_cost_fdm_PQIC( double an[], double anp1[], double gamma, double refs[] );
//void gdescent_integral_fdm_PQIC( double an[], double anp1[], double gamma, double refs[] );
//void gdescent_direct_PQIC( double an[], double anp1[], double gamma, double refs[] );

//void gdescent_cost_fdm_PPIC( double pn[], double pnp1[], double gamma, double refs[] );
//void gdescent_cost_fdm_PLIC( double an[], double anp1[], double gamma, double refs[] );

void cost_fdm_PPIC( double cn[], double refs[], Box U, double grad[], double grad_tol, double flood_tol ); //

void cost_fdm_PLIC( double an[], double refs[], Box U, double grad[], double grad_tol, double flood_tol );

void cost_fdm_PQIC( double an[], double refs[], Box U, double grad[], double grad_tol, double flood_tol ); //
void integral_fdm_PQIC( double an[], double refs[], Box U, double grad[] );
void direct_PQIC( double an[], double refs[], double grad[] );

void flooding( double a[], double vof, Box U, double tol );
void para_flooding( double c[], double vof, Box U, double tol );
void init_bases();
void init_grad_bases( Phi * phi0 );

void linear_MOF_unit( double ref[], double error_tols[], double final_coeffs[], int max_iter );
void parabolic_MOF( double ref[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[] );
void quadratic_MOF( double ref[], Box U, char guess_type, double error_tols[], double final_coeffs[], int max_iter, double init_coeffs[] );

void get_refs( Box U, double coeffs[], double refs[] );
void transform_refs_to_unit( double refs[], double refs_t[], Box U );
void transform_coeffs_to_nonunit( double coeffs[], double coeffs_t[], Box U );

#endif
