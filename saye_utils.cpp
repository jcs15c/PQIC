#include <set>
#include <cmath>
#include "saye_utils.h"
using namespace std;

Point::Point()
{
	for( int i = 0; i < 3; i++ )
		vals[i] = -1;
	d = 0;
}

Point::Point( double xs[], int n )
{
	for( int i = 0; i < n; i++ )
		vals[i] = xs[i];
	for( int i = n; i < 3; i++ )
		vals[i] = -1;
	d = n;
}

Point::Point( const Point & x )
{
	for( int i = 0; i < 3; i++ )
		vals[i] = x.vals[i];
	d = x.d;
}

void Point::saye( double val, int k )
{
	for( int i = d; i > k; i-- )
		vals[i] = vals[i-1];
	vals[k] = val;
	d++;
}

double& Point::operator[]( int i )
{	return vals[i];	  }

Point& Point::operator=( const Point x )
{
	for( int i = 0; i < 3; i++ )
		vals[i] = x.vals[i];
	d = x.d;
}

int Point::getD()
{   return d;   }

///////////////////////////////////////

Box::Box()
{	d = 0;   }

Box::Box( double nxL[], double nxU[], int nd )
{
	for( int i = 0; i < nd; i++ )
	{
		x[i][0] = nxL[i];
		x[i][1] = nxU[i];
	}

	d = nd;
}

const double* Box::operator[]( int dir ) const
{	return x[dir];   }

void Box::split( Box& U1, Box& U2 ) const
{
	int sd = 0;
	double x1L[3], x1U[3];
	double x2L[3], x2U[3];

	// Get the direction to split along
	for( int i = 0; i < d; i++ )
	{
		if( x[i][1] - x[i][0] > x[sd][1] - x[sd][0] )
			sd = i;
	}

	for( int i = 0; i < sd; i++ )
	{
		x1L[i] = x2L[i] = x[i][0];
		x1U[i] = x2U[i] = x[i][1];
	}
	
	x1L[sd] = x[sd][0];
	x1U[sd] = 0.5*( x[sd][0] + x[sd][1] );

	x2L[sd] = 0.5*( x[sd][0] + x[sd][1] );
	x2U[sd] = x[sd][1];

	for( int i = sd+1; i < d; i++ )
	{
		x1L[i] = x2L[i] = x[i][0];
		x1U[i] = x2U[i] = x[i][1];
	}
	
	U1 = Box( x1L, x1U, d );
	U2 = Box( x2L, x2U, d );
}

int Box::getD() const
{   return d;   }

double Box::volume() const
{
    double vol = 1;
    for( int i = 0; i < d; i++ )
        vol *= x[i][1] - x[i][0];
    return vol;
}

void Box::erase( Box& U, int k ) const
{
	double xL[3], xU[3];
    int nd = 0;

    for( int i = 0; i < k; i++ )
    {
        xL[nd] = x[i][0];
        xU[nd] = x[i][1];
        nd++;
    }

    for( int i = k + 1; i < d; i++ )
    {
        xL[nd] = x[i][0];
        xU[nd] = x[i][1];
        nd++;
    }

    U = Box( xL, xU, d - 1 );
}

///////////////////////////////////////
Phi2D::Phi2D()
{ }

Phi2D::Phi2D( double nCoeff[] )
{
	for( int i = 0; i < 6; i++ )
		coeff[i] = nCoeff[i];
}

Phi2D::Phi2D( const Phi2D & phi )
{
    for( int i = 0; i < 6; i++ )
		coeff[i] = phi.coeff[i];

}

double Phi2D::eval( Point x )
{
    return coeff[0]*x[0]*x[0] + coeff[1]*x[1]*x[1] + 
	       coeff[2]*x[0]*x[1] + coeff[3]*x[0]      + 
           coeff[4]*x[1]      + coeff[5];
}

double Phi2D::eval_k( Point x, int k )
{
	if( k == 0 )
        return 2*coeff[0]*x[0] + coeff[2]*x[1] + 
				 coeff[3];
    if( k == 1 )
        return 2*coeff[1]*x[1] + coeff[2]*x[0] + 
				 coeff[4];
}

set<double> Phi2D::roots( Point x, double x1, double x2, int k )
{
    double a, b, c;

	if(k == 0)
	{
		a = coeff[0];
        b = coeff[2]*x[1] + coeff[3];
        c = coeff[1]*x[1]*x[1] + coeff[4]*x[1] + 
		    coeff[5];
	}
	if(k == 1)
	{
		a = coeff[1];
        b = coeff[2]*x[0] + coeff[4];
        c = coeff[0]*x[0]*x[0] + coeff[3]*x[0] + 
		    coeff[5];
	}

	return quadratic( a, b, c, x1, x2 );
}


///////////////////////////////////////
Phi3D::Phi3D()
{ }

Phi3D::Phi3D( double nCoeff[] )
{
	for( int i = 0; i < 10; i++ )
		coeff[i] = nCoeff[i];
}

Phi3D::Phi3D( const Phi3D & phi )
{
    for( int i = 0; i < 10; i++ )
		coeff[i] = phi.coeff[i];

}

double Phi3D::eval( Point x )
{
    return coeff[0]*x[0]*x[0] + coeff[1]*x[1]*x[1] + 
	       coeff[2]*x[2]*x[2] + coeff[3]*x[0]*x[1] + 
		   coeff[4]*x[1]*x[2] + coeff[5]*x[0]*x[2] + 
           coeff[6]*x[0]      + coeff[7]*x[1]      + 
		   coeff[8]*x[2]      + coeff[9];
}

double Phi3D::eval_k( Point x, int k )
{
	if( k == 0 )
        return 2*coeff[0]*x[0] + coeff[3]*x[1] + 
				 coeff[5]*x[2] + coeff[6];
    if( k == 1 )
        return 2*coeff[1]*x[1] + coeff[3]*x[0] + 
				 coeff[4]*x[2] + coeff[7];
    if( k == 2 )
        return 2*coeff[2]*x[2] + coeff[4]*x[1] + 
				 coeff[5]*x[0] + coeff[8];
}

set<double> Phi3D::roots( Point x, double x1, double x2, int k )
{
    double a, b, c;

	if(k == 0)
	{
		a = coeff[0];
        b = coeff[3]*x[1] + coeff[5]*x[2] + coeff[6];
        c = coeff[1]*x[1]*x[1] + coeff[2]*x[2]*x[2] + 
		    coeff[4]*x[1]*x[2] + coeff[7]*x[1] + 
			coeff[8]*x[2] + coeff[9];
	}
	if(k == 1)
	{
		a = coeff[1];
        b = coeff[3]*x[0] + coeff[4]*x[2] + coeff[7];
        c = coeff[0]*x[0]*x[0] + coeff[2]*x[2]*x[2] + 
		    coeff[5]*x[0]*x[2] + coeff[6]*x[0] + 
			coeff[8]*x[2] + coeff[9];
	}
	if(k == 2)
	{
		a = coeff[2];
        b = coeff[4]*x[1] + coeff[5]*x[0] + coeff[8];
        c = coeff[0]*x[0]*x[0] + coeff[1]*x[1]*x[1] + 
		    coeff[3]*x[0]*x[1] + coeff[6]*x[0] + 
			coeff[7]*x[1] + coeff[9];
	}

	return quadratic( a, b, c, x1, x2 );
}


////////////////////////////////////////////////

Psi::Psi()
{ }

Psi::Psi( const Psi & psi )
{
    for( int i = 0; i < 2; i++ )
    {
        vals[i] = psi.vals[i];
        dirs[i] = psi.dirs[i];
    }
    
    phi0 = psi.phi0;
    d = psi.d;
}

Psi::Psi( Phi * phi )
{
    phi0 = phi;
	for( int i = 0; i < 2; i++ )
		vals[i] = dirs[i] = -1;
	d = 0;
}


Psi::Psi( Psi psi, double val, int dir )
{
	phi0 = psi.phi0;
    for( int i = 0; i < psi.getD(); i++ )
	{
		vals[i] = psi.getVal(i);
		dirs[i] = psi.getDir(i);
	}

	vals[psi.getD()] = val;
	dirs[psi.getD()] = dir;

	d = psi.getD() + 1;
}

double Psi::eval( Point x )
{
	for( int i = d - 1; i >= 0; i-- )
		x.saye( vals[i], dirs[i] );

	return phi0->eval( x );
}

double Psi::eval_k( Point x, int k )
{
	for( int i = d - 1; i >= 0; i-- )
	{
		if( dirs[i] <= k )
			k++;
		x.saye( vals[i], dirs[i] );
	}

	return phi0->eval_k( x, k );
    
}

set<double> Psi::roots( Point x, double x1, double x2, int k )
{
	x.saye( -999, k );
	for( int i = d - 1; i >= 0; i-- )
	{
		if( dirs[i] <= k )
			k++;
		x.saye( vals[i], dirs[i] );
	}

    return phi0->roots( x, x1, x2, k );
}


int Psi::getDir( int i )
{    return dirs[i];    }

double Psi::getVal( int i )
{   return vals[i];     }

int Psi::getD()
{   return d;   }

Phi* Psi::getPhi()
{   return phi0;    }
////////////////////////////////////////////////

set<double> quadratic( double a, double b, double c, double x1, double x2 )
{
	set<double> roots;
	double d = b*b - 4*a*c, r1, r2;

	if( d < 0 || ( a == b && b == 0 ) )
		return roots;
	else if( a == 0 )
	{
		r1 = -c / b;
		r2 = r1;
	}
	else
	{
		r1 = ( -b + sqrt(d) )/2/a;
		r2 = ( -b - sqrt(d) )/2/a;
	}

	if( x1 < r1 && r1 < x2 )
		roots.insert( r1 );
	if( x1 < r2 && r2 < x2 )
		roots.insert( r2 );
	
	return roots;
}

double sgn( int m, int s, bool S, int sigma )
{
    if( m == sigma*s || S == true )
        return sigma*m;
    return 0;
}

double sign( double d )
{
    if( d > 0 )
        return 1;
    else if( d < 0 )
        return -1;
    return 0;
}

double dabs( double d )
{
    if( d < 0 )
        return -d;
    return d;
}

