// Implementation of the Covariance Matrix Adaption Evolution Strategy (CMAES)
// This implements the basic strategy with box constraints using C++ and the
// Eigen package for matrix calculations
// We follow [1] for the pure CMAES extended by a "boundary handling" as
// decribed in [2].
// Both papers are linked on "The CMA Evolution Strategy" Website [3]
// This version allows to provide an optional "attractor" to resist
// in the population until the mean radius of the Gaussian bell 
// becomes small or a better individual has been sampled

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include "eigenreq.hpp"
#include "auxiliaries.hpp"
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int writeSamples( int nI, int lambda, int N, MatrixXld arx, MatrixXld arx_feas, VectorXld isColumnOutOfBounds )
{
	for ( int j = 0; j < lambda; j++ )
	{
		stringstream fileName;
		fileName << "parameters_input_" << ( nI + 1 ) << "_" << ( j + 1 ) << ".txt";
		ofstream f( fileName.str().c_str() );
		f.precision( 43 );
		for ( int i = 0; i < N; i++ )
		{
			if ( isColumnOutOfBounds( j ) ) f << arx_feas( i, j ) << "\n";
			else f << arx( i, j ) << "\n";
		}
		f.close();
	}
}

int writeScaledSamples( int nI, int lambda, int N, VectorXld lb, VectorXld ub, MatrixXld arx, MatrixXld arx_feas, VectorXld isColumnOutOfBounds )
{
	for ( int j = 0; j < lambda; j++ )
	{
		stringstream fileName;
		fileName << "parameters_input_" << ( nI + 1 ) << "_" << ( j + 1 ) << ".txt";
		ofstream f( fileName.str().c_str() );
		f.precision( 43 );
		for ( int i = 0; i < N; i++ )
		{
			if ( isColumnOutOfBounds( j ) ) f << lb( i ) + ( ub( i ) - lb( i ) ) * arx_feas( i, j ) << "\n";
			else f << lb( i ) + ( ub( i ) - lb( i ) ) * arx( i, j ) << "\n";
		}
		f.close();
	}
}

int readResults( int nI, int lambda, string inDirName, VectorXld isColumnOutOfBounds, VectorXld penalty, VectorXld &arfitness )
{
	cout << "isColumnOutOfBounds" << nI << endl << isColumnOutOfBounds << endl << endl;
	
	double bestCost = 1e20;
	int besti;
			
	for ( int i = 0; i < lambda; i++ )
	{
		char prevDir[3];
		string line;
		stringstream ss;
		stringstream fileName;
		
/*		if ( isColumnOutOfBounds( i ) )
		{
			fileName << "fitness_feas_" << nI << "_" << ( i + 1 ) << ".txt";
		}
		else*/
		
		{
			fileName << inDirName << "/misfit_output_" << nI << "_" << ( i + 1 ) << ".txt";
		}
		printf( "read results from file %s\n", fileName.str().c_str() );
		
		ifstream f( fileName.str().c_str() );
		
		double sumSquare = 0;
		
		while (getline( f, line ))
		{
			if ( (line.find( "." ) == string::npos) && (line.size() > 1) )
			{
				sumSquare += 1e12;
			}
			else if ( (line.find( "." ) == string::npos) && (line.size() <= 1) )
			{
				// Do nothing
			}
			else
			{
				double thisVal = stod(line);
				sumSquare += thisVal * thisVal;
			}
		}
		f.close();
		ss << sumSquare;
		assert( ss >> sumSquare );
		
		cout << sumSquare << endl;
		
		arfitness( i ) = sumSquare;
		
		//fC1 << arfitness( i ) << "\n";
	
		if ( isColumnOutOfBounds( i ) )
		{
			arfitness( i ) = arfitness( i ) + penalty( i );
			cout << "Out of bounds: " << arfitness( i ) << endl;
		}
		
		// Update best cost
		if (arfitness( i ) < bestCost)
		{
			bestCost = arfitness( i );
			besti = i+1;
		}
		
	}
	
	cout << arfitness << endl;
	
	ofstream fBest(inDirName + "/best_cmaes.txt");
	fBest << nI << "\n" << besti << "\n";
	fBest << arfitness;
	fBest.close();

}

int writeAlgVars( int nI, int lambda, int N, int counteval, int eigeneval, bool useAttractor, long double sigma, VectorXld pc, VectorXld ps, MatrixXld C, MatrixXld B, MatrixXld D, MatrixXld arz, MatrixXld arx, bool boundWeightsInitialized, int valHistN, int valHistSize, VectorXld valHist, VectorXld gamma, VectorXld isColumnOutOfBounds, VectorXld penalty, VectorXld xmean )
{
	stringstream fileName;
	fileName << "algVars_" << ( nI + 1 ) << ".txt";
	ofstream f( fileName.str().c_str() );
	f.precision( 43 );
	// 1.) write counteval
	f << "counteval" << nI << "\n"; 
	f << counteval << "\n";
	// 2.) write eigeneval
	f << "eigeneval" << nI << "\n"; 
	f << eigeneval << "\n";
	// 3.) write useAttractor
	f << "useAttractor" << nI << "\n";
	f << ( int )useAttractor << "\n";
	// 4.) write sigma
	f << "sigma" << nI << "\n"; 
	f << sigma << "\n";
	// 5.) write pc
	f << "pc" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		f << pc( i ) << "\n";
	}
	// 6.) write ps
	f << "ps" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		f << ps( i ) << "\n";
	}
	// 7.) write C
	f << "C" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << C( i, j ) << " ";
		}
		f << "\n";
	}
	// 8.) write B
	f << "B" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << B( i, j ) << " ";
		}
		f << "\n";
	}
	// 9.) write D
	f << "D" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << D( i, j ) << " ";
		}
		f << "\n";
	}
	// 10.) write arz
	f << "arz" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < lambda; j++ )
		{
			f << arz( i, j ) << " ";
		}
		f << "\n";
	}
	// 11.) write arx
	f << "arx" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < lambda; j++ )
		{
			f << arx( i, j ) << " ";
		}
		f << "\n";
	}
	// 12.) write boundWeightsInitialized
	f << "boundWeightsInitialized" << nI << "\n"; 
	f << boundWeightsInitialized << "\n";
	// 13.) write valHistSize
	f << "valHistSize" << nI << "\n"; 
	f << valHistSize << "\n";
	// 14.) write valHist 
	f << "valHist" << nI << "\n"; 
	for ( int i = 0; i < valHistN; i++ )
	{
		f << valHist( i ) << "\n";
	}
	// 15.) write gamma 
	f << "gamma" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		f << gamma( i ) << "\n";
	}
	// 16.) write isColumnOutOfBounds
	f << "isColumnOutOfBounds" << nI << "\n"; 
	for ( int i = 0; i < lambda; i++ )
	{
		f << isColumnOutOfBounds( i ) << "\n";
	}
	// 17.) write penalty
	f << "penatly" << nI << "\n"; 
	for ( int i = 0; i < lambda; i++ )
	{
		f << penalty( i ) << "\n";
	}
	// 18.) write xmean
	f << "xmean" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		f << xmean( i ) << "\n";
	}
	f.close();
}

int readAlgVars( int nI, int lambda, int N, string inDirName, int *counteval, int *eigeneval, bool *useAttractor, long double *sigma, VectorXld &pc, VectorXld &ps, MatrixXld &C, MatrixXld &B, MatrixXld &D, MatrixXld &arz, MatrixXld &arx, bool *boundWeightsInitialized, int valHistN, int *valHistSize, VectorXld &valHist, VectorXld &gamma, VectorXld &isColumnOutOfBounds, VectorXld &penalty, VectorXld &xmean_old )
{
	char prevDir[4];
	string bashLine;
	string bashOut;
	string line;
	stringstream ss;
	stringstream fileName;
	FILE *fp;
	char path[230];
	
  	fileName << inDirName << "/algVars_" << nI << ".txt";
	
	ifstream f( fileName.str().c_str() );
	// 1.) read counteval
	getline( f, line );
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> *counteval );
	// 2.) read eigeneval
	getline( f, line );
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> *eigeneval );
	// 3.) read usaAttractor
	getline( f, line );
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
	int bol; assert( ss >> bol ); *useAttractor = ( bool )bol;
	// 4.) read sigma
	getline( f, line );
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> *sigma );
	// 5.) read pc
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> pc( i ) );
	}
	// 6.) read ps
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> ps( i ) );
	}
	// 7.) read C
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line.c_str();
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> C( i, j ) );
		}
	}
	// 8.) read B
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> B( i, j ) );
		}
	}
	// 9.) read D
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> D( i, j ) );
		}
	}
	// 10.) read arz
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < lambda; j++ )
		{
			assert( ss >> arz( i, j ) );
		}
	}
	// 11.) read arx
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < lambda; j++ )
		{
			assert( ss >> arx( i, j ) );
		}
	}
	// 12.) read boundWeightsInitialized
	getline( f, line );
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> *boundWeightsInitialized );
	// 13.) read valHistSize
	getline( f, line );
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> *valHistSize );
	// 14.) read valHist 
	getline( f, line );
	for ( int i = 0; i < valHistN; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> valHist( i ) );
	}
	// 15.) read gamma 
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> gamma( i ) );
	}
	// 16.) read isColumnOutOfBounds 
	getline( f, line );
	for ( int i = 0; i < lambda; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> isColumnOutOfBounds( i ) );
	}
	// 17.) read penalty 
	getline( f, line );
	for ( int i = 0; i < lambda; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> penalty( i ) );
	}
	// 18.) read xmean_old 
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> xmean_old( i ) );
	}
	f.close(); 
}



int main( int argc, char* argv[] )
{

	int N, lambda, nI; // problem size, population size  and iteration number
	string functionName; // fitness function
	int srandOffset; // vsa: for different random number sequences
	bool allowAttractor; // vsa: allow attractor in population flag
	double attractorFitness; // vsa: fitness of optional attractor
 
	if ( argc < 4 ) 
	{
		printf( "program requires 3 command line arguments:\ninstance filename (string) and current iteration (int)\nthe associated file is supposed to contain:\n - problem size,\n - fitness function name,\n - lower and upper bounds,\n - random number seed, \n - attractor flag\nan example file is nIter0.dat\n and folder name of algvars_X.txt \n" );
	}
	assert( argc == 4 );
//	assert( argc == 2 );
	nI = atoi( argv[ 1 ] );
	string inFileName = argv[ 2 ];
	string inDirName = argv[3];
	string line;
	stringstream ss;
	ifstream inFile( inFileName.c_str() );
	
//	inFile >> N;
//	inFile >> lambda;
//	inFile >> functionName;
	getline( inFile, line );
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> functionName );
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> N );
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> lambda );
	getline( inFile, line );
	getline( inFile, line ); // skip iteration entry, here
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> srandOffset );
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> allowAttractor );
	getline( inFile, line );
	VectorXld lb( N );
	VectorXld ub( N );
	for ( int i = 0; i < N; i++ )
	{
		getline( inFile, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss  >> lb( i ) >> ub( i ) );
	}
	getline( inFile, line );
	VectorXld a = VectorXld::Zero( N );
	for ( int i = 0; i < N; i++ )
	{
		getline( inFile, line );
		if( allowAttractor )
		{
			ss.str( "" ); ss.clear(); ss << line;
			assert( ss  >> a( i ) );
			// normalized attractor:
			a( i ) = ( a( i ) - lb( i ) ) / ( ub( i ) - lb( i ) );
		}
	}
	getline( inFile, line );
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> attractorFitness );
	inFile.close();

	cout << "functionName = " << functionName << endl;
	cout << "N = " << N << endl;
	cout << "lambda = " << lambda << endl;
	cout << "lb:" << endl << lb << endl << endl;
	cout << "ub:" << endl << ub << endl << endl;
	cout << "srandOffset = " << srandOffset << endl;
	cout << "allowAttractor = " << allowAttractor << endl << endl;
	cout << "a:" << endl << a << endl << endl;
	cout << "f(a) = " << attractorFitness << endl;

	// derived operational parameters: Selection
//	int lambda = 4 + floor( 3 * log( N ) ); // population size, offspring number (vsa: original, now set in nIter.txt)
	int mu = floor( ( long double )lambda / 2 ); // number of parents/points for recombination
	VectorXld weights = VectorXld::LinSpaced( mu, 1, mu );
	weights = log( mu + 0.5 ) - weights.array().log();// muXone recombination weights
	weights = weights / weights.sum(); // normalize recombination weights array
	long double mueff = weights.sum() * weights.sum() / ( weights.array() * weights.array() ).sum(); // variance-effective size of mu

	// derived operational parameters: Adaptation
	long double cc = ( 4 + mueff / N ) / ( N + 4 + 2 * mueff / N ); // time constant for cumulation for C
	long double cs = ( mueff + 2 ) / ( N + mueff + 5 ); // t-const for cumulation for sigma control
//	long double c1 = 2 / ( ( N + 1.3 ) * ( N + 1.3 ) + mueff ); // learning rate for rank-one update of C (vsa: original)
//	long double cmu = 2 * ( mueff - 2 + 1 / mueff ) / ( ( N + 2 ) * ( N + 2 ) + 2 * mueff / 2 ); // and for rank-mu update (vsa: original)
	long double attenuation = 1; // vsa: attenuation factor for updating C
	long double c1 = 2 / attenuation / ( ( N + 1.3 ) * ( N + 1.3 ) + mueff ); // learning rate for rank-one update of C (vsa: attenuated)
	long double cmu = 2 / attenuation * ( mueff - 2 + 1 / mueff ) / ( ( N + 2 ) * ( N + 2 ) + 2 * mueff / 2 ); // and for rank-mu update (vsa: attenuated)
	long double damps = 1 + cs;
	long double damps_plus = 2 * sqrt( ( mueff - 1 ) / ( N + 1 ) ) - 1;
	if ( damps_plus > 0 ) damps = damps + damps_plus; // damping for sigma
	int valHistN = 20 + 3 * N / lambda; // maximum size of value history 

	// further operational parameter
	long double chiN = sqrt( N ) * ( 1 - ( long double )1 / ( 4 * N ) + ( long double )1 / ( 21 * N * N ) ); // expectation of || N( 0, I ) || == norm( randn( N, 1 ) )

	// declare dynamic (internal) strategy parameters
	bool useAttractor; // vsa: "use attractor in population flag"
	VectorXld pc( N ); // evolution path for C
	VectorXld ps( N ); // evolution path for sigma
	MatrixXld B( N, N );  // B defines coordinate system
	MatrixXld D( N, N );  // diagonal matrix D defines scaling
	MatrixXld C( N, N );  // covariance matrix
	int counteval; // count total number of function evaluations
	int eigeneval; // B and D updated at counteval == 0
	int terminate = 0; // flag, if a termination criterion is satisfied
	MatrixXld arz( N, lambda );
	MatrixXld arx( N, lambda );
	VectorXld isColumnOutOfBounds( lambda ); // indicates arx columns that are out of provided boundaries
	MatrixXld arx_feas( N, lambda );
	MatrixXld arz_sub( N, mu );
	MatrixXld arx_sub( N, mu );
	VectorXld arfitness( lambda );
	VectorXld penalty( lambda ); // penalties for samples
	VectorXld xmean( N );
	VectorXld xmean_old( N );
	VectorXld zmean( N );
	long double sigma; // coordinate wise standard deviation (step-size)
	VectorXld gamma( N ); // boundary weight vector
	VectorXld valHist( valHistN ); // init value history
	int valHistSize; // init size of value history
	bool boundWeightsInitialized; // indicator if boundary weights are already initialized

	if ( nI == 0 )
	{
		// initialize dynamic strategy parameters
		useAttractor = allowAttractor; // vsa: additional feature

//		// vsa: original (random point in unit cube)
//		xmean = VectorXld::Random( N );
//		xmean = ( xmean.array() + 1 ) / 2;

//		// vsa: center of unit cube
		xmean = VectorXld::Zero( N );
		xmean = xmean.array() + 0.5;

		// vsa: (only for one run) fixing initial mean for parameter 2 and parameter 4
//		xmean( 1 ) = 6.0 / 44.0; // => 10 in [4,48] 
//		xmean( 3 ) = 1.0 / 3.9; // => 1 in [0.1,4]

		sigma = 0.5;
		pc = VectorXld::Zero( N );
		ps = VectorXld::Zero( N );
		B = MatrixXld::Identity( N, N );
		D = MatrixXld::Identity( N, N );
		C = B * D; C = C * C.transpose();
		cout << "xmean" << nI << endl << xmean << endl << endl;
		cout << "sigma" << nI << endl << sigma << endl << endl;
		cout << "C" << nI << endl << C << endl << endl;
		cout << "D" << nI << endl << D << endl << endl;
		cout << "B" << nI << endl << B << endl << endl;
		counteval = 0;
		eigeneval = 0;
		isColumnOutOfBounds = VectorXld::Zero( lambda );
		gamma = VectorXld::Zero( N );
		valHist = VectorXld::Zero( valHistN );
		valHistSize = 0;
		boundWeightsInitialized = false;
		printf( "initialized\n" );
	}
	else // ( nI > 0 )
	{
		// read current dynamic strategy parameters and fitness values of
		// previous iteration samples
		readAlgVars( nI, lambda, N, inDirName, &counteval, &eigeneval, &useAttractor, &sigma, pc, ps, C, B, D, arz, arx, &boundWeightsInitialized, valHistN, &valHistSize, valHist, gamma, isColumnOutOfBounds, penalty, xmean_old );
		readResults( nI, lambda, inDirName, isColumnOutOfBounds, penalty, arfitness );

		// calculate already scaled deltas of [2] (11):
		long double val = ( percentile( arfitness, 75 ) - percentile( arfitness, 25 ) ) / N / C.diagonal().mean() / ( sigma * sigma ); 
		if ( valHistSize < valHistN )
		{
			valHist( valHistSize ) = val;
			valHistSize++;
		}
		else
		{
			valHist.head( valHistN - 1 ) = valHist.tail( valHistN - 1 );
			valHist( valHistN - 1 ) = val;
		}
//		cout << "valHist" << endl << valHist << endl << endl;

		// Sort by fitness and compute weighted mean into xmean
		VectorXld arindex = VectorXld::LinSpaced( lambda, 0, lambda - 1 );
		sortVectors( arfitness, arindex ); // minimization
		cout << "arfitness" << nI << endl << arfitness << endl << endl;
		cout << "arindex" << nI << endl << arindex << endl << endl;
		for ( int k = 0; k < mu; k++ )
		{
			arx_sub.col( k ) = arx.col( arindex( k ) );
			arz_sub.col( k ) = arz.col( arindex( k ) );
		}
		xmean_old = xmean;
		xmean = arx_sub * weights; // recombination [1] (38),(39)
		zmean = arz_sub * weights; // = D^ -1 * B’ * ( xmean - xmean_old ) / sigma
		cout << "xmean" << nI << endl << xmean << endl << endl;
		cout << "zmean" << nI << endl << zmean << endl << endl;
		
		// Cumulation: Update evolution paths
		ps = ( 1 - cs ) * ps + ( sqrt( cs * ( 2 - cs ) * mueff ) ) * ( B * zmean ); // [1] (40)
		int hsig = 0;
		if ( ps.norm() / sqrt( 1 - pow( 1 - cs, 2 * counteval / lambda ) ) / chiN < 1.4 + 2 / ( N + 1 ) ) hsig = 1;
		pc = ( 1 - cc ) * pc + hsig * sqrt( cc * ( 2 - cc ) * mueff ) * ( B * D * zmean ); // [1] (42)
		cout << "ps" << nI << endl << ps << endl << endl;
		cout << "hsig" << nI << endl << hsig << endl << endl;
		cout << "pc" << nI << endl << pc << endl << endl;

		// Adapt covariance matrix C
		MatrixXld A = B * D * arz_sub;
		MatrixXld M = weights.asDiagonal();
		C =	  ( 1 - c1 - cmu ) * C // [1] (43)
			+ c1 * ( pc * pc.transpose() // plus rank one update
			+ ( 1 - hsig ) * cc * ( 2 - cc ) * C ) // minor correction
			+ cmu * A * M * A.transpose(); // plus rank mu update
	
		// Adapt step-size sigma
		sigma = sigma * exp( ( cs / damps ) * ( ps.norm() / chiN - 1 ) ); // [1] (41)
		cout << "sigma" << nI << endl << sigma << endl << endl;

		// Update B and D from C
		if ( counteval - eigeneval > lambda / ( c1 + cmu ) / N / 10 ) // to achieve O( N ^ 2 )
		{
			eigeneval = counteval;
			for ( int i = 0; i < N - 1; i++ ) for ( int j = i + 1; j < N; j++ ) C( j, i ) = C( i, j ); // enforce symmetry
			SelfAdjointEigenSolver<MatrixXld> es( C );
			D = es.eigenvalues().asDiagonal();
			B = es.eigenvectors(); // B==normalized eigenvectors
			D = D.array().sqrt(); // D contains standard deviations now
			cout << "C" << nI << endl << C << endl << endl;
			cout << "B" << nI << endl << B << endl << endl;
			cout << "diag( D )'" << nI << endl << D.diagonal().transpose() << endl << endl;
		}
		// Escape flat fitness, or better terminate?
//		long double eps = 1e-12; // vsa: original
		long double eps = 1e-5; // vsa: adaption 
		if ( ( 1 - eps ) * arfitness( ceil( 0.7 * lambda ) ) < arfitness( 0 ) && arfitness( 0 ) < ( 1 + eps ) * arfitness( ceil( 0.7 * lambda ) ) )
		{
			terminate = 1;
//			sigma = sigma * exp( 0.2 + cs / damps );
//			cout << "warning: flat fitness, consider reformulating the objective" << endl;
//			cout << "arfitness( 1 ) =" << arfitness( 1 ) << "and arfitness( ceil( 0.7 * lambda ) ) =" << arfitness( ceil( 0.7 * lambda ) ) << endl;
		}
		cout << counteval << ": " << arfitness( 0 ) << endl;
	} // else ( nI > 0 )

	// mean radius of the current Gaussian bell:
	double r = sigma * pow( D.diagonal().prod(), ( double )1 / N );
	cout << "r" << nI << endl << r << endl << endl;
	
	// Update useAttractor
//	if ( nI > 0 && useAttractor && ( arfitness( 0 ) < attractorFitness || r < 0.1 ) ) useAttractor = false;
	if ( nI > 0 && useAttractor && ( arfitness( 0 ) < attractorFitness ) ) useAttractor = false;
	cout << "useAttractor" << nI << endl << useAttractor << endl << endl;

	// Generate lambda offspring
	srand( counteval + ( lambda + 1 ) * srandOffset ); // vsa: different random number sequence
	printf( "calling srand( %i )\n", counteval + srandOffset );
	for ( int k = 0; k < lambda; k++ )
	{
		for ( int i = 0; i < N; i++ ) arz( i, k ) = normaldistribution( 0.0, 1.0 );
		arx.col( k ) = xmean.array() + sigma * ( B * D * arz.col( k ) ).array(); // add mutation, [1] (37)
		if ( useAttractor && k == 0 ) // vsa: use attractor
		{
			for ( int i = 0; i < N; i++ )
			{
				arx( i, k ) = a( i ); 
			}
			arz.col( k ) = D.inverse() * B.transpose() * ( ( 1 / sigma ) * ( xmean - a ) );
		}
		counteval = counteval + 1;
	}
	cout << "arz" << nI << endl; cout << arz << endl << endl;
	cout << "arx" << nI << endl; cout << arx << endl << endl;

	// Boudary handling
	printf( "boundary handling ...\n" );
	isColumnOutOfBounds = VectorXld::Zero( lambda );
	VectorXld tx( N );
	VectorXld ti( N );
	xIntoUnitCube( xmean, tx, ti );
	if ( ti.any() ) // distribution mean out of bounds
	{
		cout << "ti" << endl << ti << endl << endl;	
		// Set initial weights if necessary:
		if ( !boundWeightsInitialized || nI == 2 )
		{	
			gamma = 2.0002 * median( valHist.head( valHistSize ) ) * VectorXld::Constant( N, 1 ); // gamma / N in [2] (11)
			cout << "gamma = " << gamma << endl;
			boundWeightsInitialized = true;
		}
		// Increase weights:
		tx = xmean - tx;
		VectorXld idx = VectorXld::Zero(N );
		idx = (    ti.array() != 0
			&& tx.array().abs() > 3 * max( ( long double )1, ( long double )sqrt( valHistN ) / mueff ) * sigma * C.diagonal().array().sqrt()
			&& tx.array() * ( xmean - xmean_old ).array() > 0
		      ).select( 1, idx );
		cout << "idx" << endl << idx << endl << endl;
		if ( idx.any() )  // increase
		{
			gamma = ( idx ).select( pow( min( ( long double )1, mueff / 10 / N ), ( long double )1.2 ) * gamma, gamma ); 
		}
		// calculate corrected (feasible) offspring, if necessary
		arx_feas = MatrixXld::Zero( N, lambda );
		isColumnOutOfBounds = VectorXld::Zero( lambda );
		for ( int k = 0; k < lambda; k++ )
		{
			VectorXld dummy( N );
			VectorXld v( N );
			xIntoUnitCube( arx.col( k ), v, dummy );
			arx_feas.col( k ) = v;
			if ( dummy.any() )
			{
				isColumnOutOfBounds( k ) = 1;
				counteval = counteval + 1;
			}
		}
		VectorXld xi( N );// scaling vector
		xi = ( 0.9 * ( C.diagonal().array().log() - C.diagonal().array().log().mean() ) ).exp(); 
		penalty = ( gamma.array() / xi.array() ).matrix().transpose() * ( ( arx_feas - arx ).array() * ( arx_feas - arx ).array() ).matrix(); // vsa: original
//		penalty = 2 * ( gamma.array() / xi.array() ).matrix().transpose() * ( ( arx_feas - arx ).array() * ( arx_feas - arx ).array() ).matrix(); // vsa: multiplying original with constant factor to increase boundary precision
//		cout << "arfitness_feas" << endl << arfitness_feas << endl << endl;
//		cout << "gamma" << endl << gamma << endl << endl;
		cout << "xi" << endl << xi << endl << endl;
		cout << "arx_feas" << endl << arx_feas << endl << endl;
		cout << "penalty" << endl << penalty << endl << endl;
	} // mean out of bounds case

	writeScaledSamples( nI, lambda, N, lb, ub, arx, arx_feas, isColumnOutOfBounds );
	// writeSamples( nI, lambda, N, arx, arx_feas, isColumnOutOfBounds );
	writeAlgVars( nI, lambda, N, counteval, eigeneval, useAttractor, sigma, pc, ps, C, B, D, arz, arx, boundWeightsInitialized, valHistN, valHistSize, valHist, gamma, isColumnOutOfBounds, penalty, xmean );
	// update nIter.txt
	ofstream f;
	f.open( "nIter.txt", ofstream::out | ofstream::app );
	f << nI + 1 << " " << arfitness( 0 ) << " " << terminate << "\n";
	f.close();
}

/* -------------------------------- Literature --------------------------------

[1] N. Hansen (2011). The CMA Evolution Strategy: A Tutorial

[2] N. Hansen, A.S.P. Niederberger, L. Guzzella and P. Koumoutsakos (2009). A Method for Handling Uncertainty in Evolutionary Optimization with an Application to Feedback Control of Combustion. IEEE Transactions on Evolutionary Computation, 13(1), pp. 180-197

[3] CMAES website: https://www.lri.fr/~hansen/cmaesintro.html

-----------------------------------------------------------------------------*/
