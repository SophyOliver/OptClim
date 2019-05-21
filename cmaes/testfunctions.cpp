#include "eigenreq.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>

using namespace std;

long double griewank( VectorXld x )
{
	int n = x.size();
	long double obj = 0;
	long double amp_scale = 0.25;
	long double freq_scale = 63.24555320336759; // sqrt( 4000 )
	for ( int j = 0; j < n; j++ )
	{
		obj = obj + x( j ) * x( j );
	}
	long double prod = 1;
	for ( int j = 0; j < n; j++ )
	{
		prod = prod * cos( freq_scale / ( j + 1 ) * x( j ) );
	}
	obj = obj + amp_scale * ( 1 - prod );
	return( obj );
}

long double rosenbrock( VectorXld x )
{
	int n = x.size();
	long double obj = 0;
	long double q1, q2, q3, d;
	for ( int j = 0; j < n - 1; j++ )
	{
		q1 = x( j ) * x( j );
		q2 = q1 + 1 - 2 * x( j );
		d = x( j + 1 ) - q1;
		q3 = d * d;
		obj = obj + q2 + 100 * q3;
	}
	return( obj );
}



long double sphere( VectorXld x )
{
	return( ( x.array() * x.array() ).sum() );
}



long double ssphere( VectorXld x )
{
	return( pow( ( x.array() * x.array() ).sum(), ( long double )0.5 ) );
}
  


long double schwefel( VectorXld x )
{
	int n = x.size();
	long double r = 0;
	for ( int i = 1; i <= n; i++ )
	{
		r = r + ( x.head( i ).sum() ) * ( x.head( i ).sum() );
	}
	return( r );
}



long double cigar( VectorXld x )
{
	int n = x.size();
	return( x( 0 ) * x( 0 ) + 1e6 * ( x.tail( n - 1 ).array() * x.tail( n - 1 ).array() ).sum() );
}  



long double cigtab( VectorXld x )
{
	int n = x.size();
	return( x( 0 ) * x( 0 ) + 1e8 * x( n - 1 ) * x( n - 1 ) + 1e4 * ( x.segment( 1, n - 2 ).array() * x.segment( 1, n - 2 ).array() ).sum() );
}  



long double tablet( VectorXld x )
{
	int n = x.size();
	return ( 1e6 * x( 0 ) * x( 0 ) + ( x.tail( n - 1 ).array() * x.tail( n - 1 ).array() ). sum() );
}



long double elli( VectorXld x )
{
	int n = x.size();
	VectorXld v = VectorXld::LinSpaced( n, 0, 1 );
	for ( int i = 0; i < n; i++ ) v( i ) = pow( ( long double )1e6, v( i ) );
	return( v.transpose() * ( x.array() * x.array() ).matrix() );
}



long double elli100( VectorXld x )
{
	int n = x.size();
	VectorXld v = VectorXld::LinSpaced( n, 0, 1 );
	for ( int i = 0; i < n; i++ ) v( i ) = pow( ( long double )1e4, v( i ) );
	return( v.transpose() * ( x.array() * x.array() ).matrix() );
}



long double plane( VectorXld x )
{
	return( x( 0 ) );
}  



long double twoaxes( VectorXld x )
{
	int n = x.size();
	int k = floor( n / 2 );
	VectorXld v1 = x.head( k ).array() * x.head( k ).array();
	VectorXld v2 = x.tail( n - k ).array() * x.tail( n - k ).array();
	return( v1.sum() + 1e6 * v2.sum() );
}



long double parabR( VectorXld x )
{
	int n = x.size();
	VectorXld v = x.tail( n - 1 ).array() * x.tail( n - 1 ).array();
	return( -x( 0 ) + 100 * v.sum() );
}



long double sharpR( VectorXld x )
{
	int n = x.size();
	return( -x( 0 ) + 100 * x.tail( n - 1 ).norm() );
}
 


long double diffpow( VectorXld x )
{
	int n = x.size();
	VectorXld v = VectorXld::LinSpaced( n, 0, 10 );
	v = 2 + v.array();
	for ( int i = 0; i < n; i++ ) v( i ) = pow( abs( x( i ) ), v( i ) );
	return( v.sum() );
} 
  


long double rastrigin10( VectorXld x )
{
	int n = x.size();
	long double pi = 3.141592653589793238462643383279503;
	VectorXld v = VectorXld::LinSpaced( n, 0, 1 );
	for ( int i = 0; i < n; i++ ) v( i ) = pow( ( long double )10, v( i ) );
	v = v.array() * x.array();
	VectorXld w( n );
	for ( int i = 0; i < n; i++ ) w( i ) = 10 * cos( 2 * pi * v( i ) );
	return( 10 * n + ( ( v.array() * v.array() ).matrix() - w ).sum() );
}



int main( int argc, char *argv[] )
{
    // command line arguments: nPar nIter nSession fName
	assert( argc == 4 );
//	int nParams = atoi( argv[ 1 ] );
//	string fName = argv[ 2 ];
//	int nIter = atoi( argv[ 3 ] );
//	int nSession = atoi( argv[ 4 ] );
//	int kind = atoi( argv[ 5 ] );
	string instanceFileName = argv[ 1 ];
	string inFileName = argv[ 2 ];
	string outFileName = argv[ 3 ];
	string fName;
	int nParams;
	ifstream instanceFile( instanceFileName.c_str() );
//	inFile >> N;
//	inFile >> lambda;
//	inFile >> functionName;
	string line;
	stringstream ss;
	getline( instanceFile, line );
	getline( instanceFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> fName );
	getline( instanceFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> nParams );
	instanceFile.close();
//	cout << "nParams" << endl << nParams << endl << endl;
	VectorXld x( nParams );
	long double v;
//	stringstream ss;
//	ss << "parameters_" << nIter << "_" << nSession << "_" << kind << ".txt";
//	ifstream inFile( ss.str().c_str() );
	ifstream inFile( inFileName.c_str() );
	int i = 0;
	while ( inFile >> v && i < nParams )
	{
		x( i ) = v;
		i++;
	}
	inFile.close();
	assert ( i == nParams );
	long double obj;
	if ( fName == "rastrigin10" ) obj = rastrigin10( x );
	else if ( fName == "griewank" ) obj = griewank( x );
	else if ( fName == "rosenbrock" ) obj = rosenbrock( x );
	else if ( fName == "sphere" ) obj = sphere( x );
	else if ( fName == "ssphere" ) obj = ssphere( x );
	else if ( fName == "schwefel" ) obj = schwefel( x );
	else if ( fName == "cigar" ) obj = cigar( x );
	else if ( fName == "cigtab" ) obj = cigtab( x );
	else if ( fName == "tablet" ) obj = tablet( x );
	else if ( fName == "elli" ) obj = elli( x );
	else if ( fName == "elli100" ) obj = elli100( x );
	else if ( fName == "plane" ) obj = plane( x );
	else if ( fName == "twoaxes" ) obj = twoaxes( x );
	else if ( fName == "parabR" ) obj = parabR( x );
	else if ( fName == "sharpR" ) obj = sharpR( x );
	else if ( fName == "diffpow" ) obj = diffpow( x );
	else obj = rosenbrock( x );
//	ss.str( "" ); ss.clear();
//	ss << "fitness_" << nIter << "_" << nSession << "_" << kind << ".txt";
//	cout << "testfunction evaluating parameterss_" << nIter << "_" << nSession << " is waiting, now" << endl;
//	usleep( ( ( rand() % 5 ) + 1 ) * 1000000 ); // wait some seconds to simulate slow process
	usleep( 1 * 1000000 ); // wait 10 seconds to simulate slow process
//	cout << "testfunction evaluating parameters_" << nIter << "_" << nSession << "_" << kind << " is ready and will write fitness_" << nIter << "_" << nSession << "_" << kind << ", now" << endl;
	cout << "testfunction evaluating " << inFileName << " is ready and will write " << outFileName << ", now" << endl;
	ofstream outFile( outFileName.c_str() );
	outFile.precision( 43 );
	outFile << obj << "\n";
	outFile.close();
	return( 0 );
}
