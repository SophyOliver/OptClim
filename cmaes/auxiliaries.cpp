#include "auxiliaries.hpp"

double normaldistribution( double m, double s )
{
	double r1 = ( double )rand() / RAND_MAX;
	double r2 = ( double )rand() / RAND_MAX;
	double pi = 3.141592653589793;
	return ( s * sqrt( -2 * log( r1 ) ) * sin( 2 * pi * r2 ) + m );
}

void merge( double *v, double *x, int a, int m, int b )
// merge two consecutive segments of a vector 
// and arrange a second vector correspondingly
// 
// v, x     vectors
// a, m, b  component indizes ( counting from 0 )
//
// both subvectors of consecutive components from a up to m 
// and from m+1 up to b of v are supposed to be in ascending order 
// merging brings all components from a up to b in ascending order
// and arranges the components of x correspondingly  
{
	double	*w, *y;
	int	i, j, k;
	w = ( double* ) calloc( m + 1 - a, sizeof( double ) );
	y = ( double* ) calloc( m + 1 - a, sizeof( double ) );
	// copy left segments of v and x to auxially vectors w and y:
	for( i = 0; i < m + 1 - a; i++ )
	{
		*( w + i ) = *( v + a + i );
		*( y + i ) = *( x + a + i );
	}
	i = 0;
	j = m + 1;
	k = a;
	// copy back the next greater elements:
	while( k < j && j <= b )
	{
		if( *( w + i ) <= *( v + j ) )
		{
			*( v + k ) = *( w + i );
			*( x + k ) = *( y + i );
			k++;
			i++;
		} else
		{
			*( v + k ) = *( v + j );
			*( x + k ) = *( x + j );
			k++;
			j++;	
		} 
	}
	// copy back remaining elements of w:
	while ( k < j )
	{
		*( v + k ) = *( w + i );
		*( x + k ) = *( y + i );
		k++;
		i++;
	}
	free( w );
	free( y );
}

void mergeSort( double *v, double *x, int a, int b )
// sort one double array and permute a second double array accordingly
//
// v: array to be sorted in ascending order
// x: array to be permuted like v
// a,b: components a to b (0 <= a <= b <= length(v) - 1 ) are sorted
{
	int m;
	if( a < b )
	{
		m = ( a + b ) / 2;
		mergeSort( v, x, a, m );
		mergeSort( v, x, m + 1, b );
		merge( v, x, a, m, b );
	}	
}

void sortVectors( VectorXld &v, VectorXld &x )
// sort a vector and permute a second vector accordingly
{
	int n = v.size();
	// copy vectors to some double arrays
	double *va, *xa;
	va = ( double* )calloc( n, sizeof( double ) );
	xa = ( double* )calloc( n, sizeof( double ) );
	for ( int i = 0; i < n; i++ )
	{
		va[ i ] = v( i );
		xa[ i ] = x( i );
	}
	// sort double arrays
	mergeSort( va, xa, 0, n - 1 );
	// copy back double arrays to vectors
	for ( int i = 0; i < n; i++ )
	{
		v( i ) = va[ i ];
		x( i ) = xa[ i ];
	}
	free( va );
	free( xa );
}

double median( VectorXld v ) // Computes median of a vector v
{
	int n = v.size();
	VectorXld dummy( n );
	sortVectors( v, dummy );
	if ( n % 2 == 1 ) return( v( n / 2 ) );
	else return( 0.5 * ( v( n / 2 - 1 ) + v( n / 2 ) ) );
}

double percentile( VectorXld v, double p ) // Computes the percentile p from vector v
{
	int n = v.size();
	double r;	
	VectorXld dummy( n );
	sortVectors( v, dummy );
	if ( p <= 100 * ( 0.5 / n ) ) r = v( 0 );
	else if ( p >= 100 * ( ( n - 0.5 ) / n ) ) r = v( n - 1 );
	else
	{
		// find largest index smaller than required percentile
		VectorXld indices = VectorXld::LinSpaced( n, 1, n );
		VectorXld a = 100 * ( indices.array() - 0.5 ) / n;
		int i = ( a.array() >= p ).select( 0, indices ).lpNorm<Eigen::Infinity>();
		// interpolate linearly
		r = v( i ) + ( v( i + 1 ) - v( i ) ) * ( p - a( i ) ) / ( a( i + 1 ) - a( i ) );
	}
	return( r );
}

void xintobounds( VectorXld x, VectorXld lb, VectorXld ub, VectorXld &bx, VectorXld &ix )
// brings a vector x between lower and upper bound vectors lb and ub
// bx is the result, ix indicates the out of bounds components of x  
{
	ix = VectorXld::Zero( ix.size() );
	VectorXld ix2 = VectorXld::Zero( ix.size() );
	ix = ( x.array() < lb.array() ).select( 1, ix );
	bx = ( x.array() < lb.array() ).select( lb, x );
	ix2 = ( bx.array() > ub.array() ).select( 1, ix2 );
	bx = ( bx.array() > ub.array() ).select( ub, bx );
	ix = ix2 - ix;
}

void xIntoUnitCube( VectorXld x, VectorXld &bx, VectorXld &ix )
// map a vector x to unit cube (maping outside components onto closest bound)
// bx is the result, ix indicates the out of bounds components of x  
{
	ix = VectorXld::Zero( ix.size() );
	VectorXld ix2 = VectorXld::Zero( ix.size() );
	ix = ( x.array() < 0 ).select( 1, ix );
	bx = ( x.array() < 0 ).select( 0, x );
	ix2 = ( bx.array() > 1 ).select( 1, ix2 );
	bx = ( bx.array() > 1 ).select( 1, bx );
	ix = ix2 - ix;
}
