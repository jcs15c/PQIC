#ifndef SAYEUTILS
#define SAYEUTILS

#include <set>
using namespace std;
 
class Point
{
public:
	Point();
	Point( double xs[], int n );
	Point( const Point & x );
	void saye( double val, int k );
	double& operator[]( int i );
    Point& operator=( const Point x );

    int getD();
private:
	double vals[3];
	int d;
};

class Box
{
public:
	Box();
    Box( const Box & U );
	Box( double xL[], double xU[], int d );
	const double* operator[]( int dir ) const;
	void split( Box& U1, Box& U2 ) const;
    int getD() const;
    double volume() const;
    void erase( Box& U, int k ) const;

private:
	double x[3][2];
	int d;
};

class Phi
{
public:
    virtual double eval( Point x ) = 0;
    virtual double eval_k( Point x, int k ) = 0;
    virtual set<double> roots( Point x, double x1,
                               double x2, int k ) = 0;
};

class Phi3D : public Phi
{
public:
    Phi3D( );
	Phi3D(double coeff[]);
	Phi3D( const Phi3D & phi );
    double eval( Point x );
	double eval_k( Point x, int k );
	set<double> roots( Point x, double x1, double x2, int k );

private:
	double coeff[10];
};

class Phi2D : public Phi
{
public:
    Phi2D( );
	Phi2D(double coeff[]);
	Phi2D( const Phi2D & phi );
    double eval( Point x );
	double eval_k( Point x, int k );
	set<double> roots( Point x, double x1, double x2, int k );
   
private:
	double coeff[6];
};

class Psi
{
public:
    Psi();
	Psi( const Psi & psi );
	Psi( Phi * phi );
	Psi( Psi psi, double val, int dir );
	
    double eval( Point x );
	double eval_k( Point x, int k );
	set<double> roots( Point x, double x1, double x2, int k );

	int getDir(int i);
	double getVal(int i);
	
    int getD();
    Phi* getPhi();

private:
	double vals[2], dirs[2];
	Phi* phi0;
	int d;
};

set<double> quadratic( double a, double b, double c, double x1, double x2 );
double sgn( int m, int s, bool S, int sigma );
double sign( double d );

void print_arr( double a[], int N );

double dabs( double d );

void poly_coefficients( double a[], double c[] ); 
void line_to_param( double coeff[], double param[] );

void scaling( double a[] );
#endif

