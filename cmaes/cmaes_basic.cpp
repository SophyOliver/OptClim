// Implementation of the Covariance Matrix Adaption Evolution Strategy (CMAES)
// This implements the basic strategy using C++ and the Eigen package for
// matrix calculations
// We follow [1] for the pure CMAES, linked on "The CMA Evolution Strategy"
// Website [2]

#include <iostream>
#include <fstream>
#include "eigenreq.hpp"
#include "auxiliaries.hpp"

using namespace std;

int writeSamples( int nI, int lambda, int N, MatrixXld arx )
{
	for ( int j = 0; j < lambda; j++ )
	{
		stringstream fileName;
		fileName << "parameters_" << ( nI + 1 ) << "_" << ( j + 1 ) << ".txt";
		ofstream f( fileName.str().c_str() );
		f.precision( 43 );
		for ( int i = 0; i < N; i++ )
		{
			f << arx( i, j ) << "\n";
		}
		f.close();
	}
}

int readResults( int nI, int lambda, VectorXld &arfitness )
{
	for ( int i = 0; i < lambda; i++ )
	{
		string line;
		stringstream ss;
		stringstream fileName;
		fileName << "fitness_" << nI << "_" << ( i + 1 ) << ".txt";
		ifstream f( fileName.str().c_str() );
		getline( f, line );
		ss << line;
		assert( ss >> arfitness( i ) );
		f.close();
	}
}

int writeAlgVars( int nI, int lambda, int N, int counteval, int eigeneval, long double sigma, VectorXld pc, VectorXld ps, MatrixXld C, MatrixXld B, MatrixXld D, MatrixXld arz, MatrixXld arx )
{
	stringstream fileName;
	fileName << "algVars_" << ( nI + 1 ) << ".txt";
	ofstream f( fileName.str().c_str() );
	f.precision( 43 );
	// 1.) write counteval
	f << counteval << "\n";
	// 2.) write eigeneval
	f << eigeneval << "\n";
	// 3.) write sigma
	f << sigma << "\n";
	// 4.) write pc
	for ( int i = 0; i < N; i++ )
	{
		f << pc( i ) << "\n";
	}
	// 5.) write ps
	for ( int i = 0; i < N; i++ )
	{
		f << ps( i ) << "\n";
	}
	// 6.) write C
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << C( i, j ) << " ";
		}
		f << "\n";
	}
	// 7.) write B
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << B( i, j ) << " ";
		}
		f << "\n";
	}
	// 8.) write D
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << D( i, j ) << " ";
		}
		f << "\n";
	}
	// 9.) write arz
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < lambda; j++ )
		{
			f << arz( i, j ) << " ";
		}
		f << "\n";
	}
	// 10.) write arx
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < lambda; j++ )
		{
			f << arx( i, j ) << " ";
		}
		f << "\n";
	}
	f.close();
}

int readAlgVars( int nI, int lambda, int N, int *counteval, int *eigeneval, long double *sigma, VectorXld &pc, VectorXld &ps, MatrixXld &C, MatrixXld &B, MatrixXld &D, MatrixXld &arz, MatrixXld &arx )
{
	string line;
	stringstream ss;
	stringstream fileName;
	fileName << "algVars_" << nI << ".txt";
//	printf( "reading from %s\n", fileName.str().c_str() );
	ifstream f( fileName.str().c_str() );
	// 1.) read counteval
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
//	printf( "counteval as line is %s\n", ss.str().c_str() );
	assert( ss >> *counteval );
	// 2.) read eigeneval
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
//	printf( "eigeneval as line is %s\n", ss.str().c_str() );
	assert( ss >> *eigeneval );
	// 3.) read sigma
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
//	printf( "sigma as line is %s\n", ss.str().c_str() );
	assert( ss >> *sigma );
	// 4.) read pc
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
//		printf( "line %i of pc as line is %s\n", i + 1, ss.str().c_str() );
		assert( ss >> pc( i ) );
	}
	// 5.) read ps
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
//		printf( "line %i of ps as line is %s\n", i + 1, ss.str().c_str() );
		assert( ss >> ps( i ) );
	}
	// 6.) read C
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line.c_str();
//		printf( "line %i of C as line is %s\n", i + 1, ss.str().c_str() );
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> C( i, j ) );
		}
	}
	// 7.) read B
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
//		printf( "line %i of B as line is %s\n", i + 1, ss.str().c_str() );
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> B( i, j ) );
		}
	}
	// 8.) read D
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
//		printf( "line %i of D as line is %s\n", i + 1, ss.str().c_str() );
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> D( i, j ) );
		}
	}
	// 9.) read arz
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
//		printf( "line %i of arz as line is %s\n", i + 1, ss.str().c_str() );
		for ( int j = 0; j < lambda; j++ )
		{
			assert( ss >> arz( i, j ) );
		}
	}
	// 10.) read arx
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
//		printf( "line %i of arx as line is %s\n", i + 1, ss.str().c_str() );
		for ( int j = 0; j < lambda; j++ )
		{
			assert( ss >> arx( i, j ) );
		}
	}
	f.close(); 
}



int main( int argc, char* argv[] )
{
	// command line arguments: problem dimension, objective function name, current iteration 

	// ----------------- process command line arguments -------------------
	int N, nI;
	string functionName;
	if ( argc != 4 ) 
	{
		printf( "Program cmaes_basic requires 3 command line arguments:\nproblem size (integer), function name (string), current iteration (int)\n\n" );
	}
	else
	{
		N = atoi( argv[ 1 ] );
		functionName = argv[ 2 ];
		nI = atoi( argv[ 3 ] );
	}
	assert( argc == 4 );

	// derived operational parameters: Selection
	int lambda = 4 + floor( 3 * log( N ) ); // population size, offspring number
	int mu = floor( ( double )lambda / 2 ); // number of parents/points for recombination
	VectorXld weights = VectorXld::LinSpaced( mu, 1, mu );
	weights = log( mu + 0.5 ) - weights.array().log();// muXone recombination weights
	weights = weights / weights.sum(); // normalize recombination weights array
	double mueff = weights.sum() * weights.sum() / ( weights.array() * weights.array() ).sum(); // variance-effective size of mu

	// derived operational parameters: Adaptation
	double cc = ( 4 + mueff / N ) / ( N + 4 + 2 * mueff / N ); // time constant for cumulation for C
	double cs = ( mueff + 2 ) / ( N + mueff + 5 ); // t-const for cumulation for sigma control
	double c1 = 2 / ( ( N + 1.3 ) * ( N + 1.3 ) + mueff ); // learning rate for rank-one update of C
	double cmu = 2 * ( mueff - 2 + 1 / mueff ) / ( ( N + 2 ) * ( N + 2 ) + 2 * mueff / 2 ); // and for rank-mu update
	double damps = 1 + cs;
	double damps_plus = 2 * sqrt( ( mueff - 1 ) / ( N + 1 ) ) - 1;
	if ( damps_plus > 0 ) damps = damps + damps_plus; // damping for sigma

	// further operational parameter
	double chiN = sqrt( N ) * ( 1 - ( double )1 / ( 4 * N ) + ( double )1 / ( 21 * N * N ) ); // expectation of || N( 0, I ) || == norm( randn( N, 1 ) )

	// declare dynamic (internal) strategy parameters
	VectorXld pc( N ); // evolution path for C
	VectorXld ps( N ); // evolution path for sigma
	MatrixXld B( N, N );  // B defines coordinate system
	MatrixXld D( N, N );  // diagonal matrix D defines scaling
	MatrixXld C( N, N );  // covariance matrix
	int counteval; // count total number of function evaluations
	int eigeneval; // B and D updated at counteval == 0
	MatrixXld arz( N, lambda );
	MatrixXld arx( N, lambda );
	MatrixXld arz_sub( N, mu );
	MatrixXld arx_sub( N, mu );
	VectorXld arfitness( lambda );
	VectorXld xmean( N );
	VectorXld zmean( N );
	long double sigma; // coordinate wise standard deviation (step-size)

	if ( nI == 0 )
	{
		// initialize dynamic strategy parameters
		xmean = VectorXld::Random( N );
		xmean = ( xmean.array() + 1 ) / 2;
		sigma = 0.5;
		pc = VectorXld::Zero( N );
		ps = VectorXld::Zero( N );
		B = MatrixXld::Identity( N, N );
		D = MatrixXld::Identity( N, N );
		C = B * D; C = C * C.transpose();
		counteval = 0;
		eigeneval = 0;
//		cout << "xmean" << nI << endl; cout << xmean << endl << endl;
//		cout << "sigma" << nI << endl; cout << sigma << endl << endl;
//		cout << "C" << nI << endl; cout << C << endl << endl;
//		cout << "B" << nI << endl; cout << B << endl << endl;
//		cout << "D" << nI << endl; cout << D << endl << endl;
	}
	else // ( nI > 0 )
	{
		// read current dynamic strategy parameters
		// and evaluation results of previous iteration
	 	readAlgVars( nI, lambda, N, &counteval, &eigeneval, &sigma, pc, ps, C, B, D, arz, arx );
		readResults( nI, lambda, arfitness );
	
		// Sort by fitness and compute weighted mean into new xmean
		VectorXld arindex = VectorXld::LinSpaced( lambda, 0, lambda - 1 );
		sortVectors( arfitness, arindex ); // minimization
//		cout << "arfitness" << nI << endl; cout << arfitness << endl << endl;
//		cout << "arindex" << nI << endl; cout << arindex << endl << endl;
		for ( int k = 0; k < mu; k++ )
		{
			arx_sub.col( k ) = arx.col( arindex( k ) );
			arz_sub.col( k ) = arz.col( arindex( k ) );
		}
		xmean = arx_sub * weights; // recombination, [1] (38),(39)
		zmean = arz_sub * weights; //== D^ -1 * B’ * ( xmean - xold ) / sigma
//		cout << "xmean" << nI << endl; cout << xmean << endl << endl;
//		cout << "zmean" << nI << endl; cout << zmean << endl << endl;
		
		// Cumulation: Update evolution paths
		ps = ( 1 - cs ) * ps + ( sqrt( cs * ( 2 - cs ) * mueff ) ) * ( B * zmean ); // [1] (40)
		int hsig = 0;
		if ( ps.norm() / sqrt( 1 - pow( 1 - cs, 2 * counteval / lambda ) ) / chiN < 1.4 + 2 / ( N + 1 ) ) hsig = 1;
		pc = ( 1 - cc ) * pc + hsig * sqrt( cc * ( 2 - cc ) * mueff ) * ( B * D * zmean ); // [1] (42)
//		cout << "ps" << nI << endl; cout << ps << endl << endl;
//		cout << "hsig" << nI << endl; cout << hsig << endl << endl;
//		cout << "pc" << nI << endl; cout << pc << endl << endl;

		// Adapt covariance matrix C
		MatrixXld A = B * D * arz_sub;
		MatrixXld M = weights.asDiagonal();
		C =	  ( 1 - c1 - cmu ) * C // [1] (43)
			+ c1 * ( pc * pc.transpose() // plus rank one update
			+ ( 1 - hsig ) * cc * ( 2 - cc ) * C ) // minor correction
			+ cmu * A * M * A.transpose(); // plus rank mu update
	
		// Adapt step-size sigma
		sigma = sigma * exp( ( cs / damps ) * ( ps.norm() / chiN - 1 ) ); // [1] (41)
//		cout << "sigma" << nI << endl; cout << sigma << endl << endl;

		// Update B and D from C
		if ( counteval - eigeneval > lambda / ( c1 + cmu ) / N / 10 )// to achieve O( N ^ 2 )
		{
			eigeneval = counteval;
			for ( int i = 0; i < N - 1; i++ ) for ( int j = i + 1; j < N; j++ ) C( j, i ) = C( i, j ); // enforce symmetry
//			cout << "C" << nI << endl; cout << C << endl << endl;
			SelfAdjointEigenSolver<MatrixXld> es( C );
			D = es.eigenvalues().asDiagonal();
			B = es.eigenvectors(); // B==normalized eigenvectors
			D = D.array().sqrt(); // D contains standard deviations now
		}
//		cout << "B" << nI << endl; cout << B << endl << endl;
//		cout << "D" << nI << endl; cout << D << endl << endl;
		// Escape flat fitness, or better terminate?
		if ( arfitness( 1 ) == arfitness( ceil( 0.7 * lambda ) ) )
		{
			sigma = sigma * exp( 0.2 + cs / damps );
			cout << "warning: flat fitness, consider reformulating the objective" << endl;
			cout << "arfitness( 1 ) =" << arfitness( 1 ) << "and arfitness( ceil( 0.7 * lambda ) ) =" << arfitness( ceil( 0.7 * lambda ) ) << endl;
		}
		cout << counteval << ": " << arfitness( 1 ) << endl;
	}// else ( nI > 0 )

	// Generate lambda offspring
	srand( counteval );
	printf( "calling srand( %i )\n", counteval );
	for ( int k = 0; k < lambda; k++ )
	{
		for ( int i = 0; i < N; i++ ) arz( i, k ) = normaldistribution( 0.0, 1.0 );
		arx.col( k ) = xmean.array() + sigma * ( B * D * arz.col( k ) ).array(); // add mutation, [1] (37)
		counteval = counteval + 1;
	}
//	cout << "arz" << nI << endl; cout << arz << endl << endl;
//	cout << "arx" << nI << endl; cout << arx << endl << endl;

	writeSamples( nI, lambda, N, arx );
	writeAlgVars( nI, lambda, N, counteval, eigeneval, sigma, pc, ps, C, B, D, arz, arx );
}

/* -------------------------------- Literature --------------------------------

[1] N. Hansen (2011). The CMA Evolution Strategy: A Tutorial

[2] CMAES website: https://www.lri.fr/~hansen/cmaesintro.html

-----------------------------------------------------------------------------*/
