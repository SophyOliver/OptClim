
#include <stdio.h>
#include <math.h>
#define SINGLE_MISFIT 0 // if 1 produces single misfit for BOBYQA (no longer needed), if 0 produces multiple misfits for DFOLS (and now BOBYQA!)

int main()
{

	// OBSERVED TRACERS

	// Read in obs data
    float obs1[10], obs2[10], obs3[10], x[10];

    FILE *obs1File, *obs2File, *obs3File, *xFile;
    obs1File = fopen("obs1.txt", "r");
    obs2File = fopen("obs2.txt", "r");
    obs3File = fopen("obs3.txt", "r");
    xFile = fopen("x.txt", "r");

    int i;
    for ( i=0; i<10; i++ ) {
        fscanf(obs1File, "%f", &obs1[i]);
        fscanf(obs2File, "%f", &obs2[i]);
        fscanf(obs3File, "%f", &obs3[i]);
        fscanf(xFile, "%f", &x[i]);
    }

	fclose(obs1File);
	fclose(obs2File);
	fclose(obs3File);
	fclose(xFile);
	
    // MODELLED TRACERS

    // Read input trial model parameters
    float test_params[2] ;

    FILE *paramFile;
    paramFile = fopen("parameters_input.txt", "r");

	printf("Test Parameters = ");
    for ( i=0; i<2; i++ ) {
        fscanf(paramFile, "%f", &test_params[i]);
        printf("%f ", test_params[i]);
    }
    printf("\n\n");
    
    fclose(paramFile);

    double tparam1 = test_params[0];
    double tparam2 = test_params[1];

    // Create modelled tracers
    double mod1[10], mod2[10], mod3[10];

    printf("Model Outputs \nMod1 = ");
    for ( i=0; i<10; i++ ) {
    	mod1[i] = tparam2 * x[i] + tparam1;
    	printf("%f ", mod1[i]);
    }
    printf("\n");

    printf("Mod2 = ");
    for ( i=0; i<10; i++ ) {
    	mod2[i] = (tparam1/10) * sin(x[i]) * tparam2;
    	printf("%f ", mod2[i]);
    }
    printf("\n");

    printf("Mod3 = ");
    for ( i=0; i<10; i++ ) {
    	mod3[i] = (tparam1/15) * cos(x[i]) * tparam2;
    	printf("%f ", mod3[i]);
    }
    printf("\n\n");

    // MISFIT

    double misfit1 = 0;
    double misfit2 = 0;
    double misfit3 = 0;

    // Sum the squared differences
    for ( i=0; i<10; i++ ) {
    	misfit1 = misfit1 + ( (obs1[i] - mod1[i]) * (obs1[i] - mod1[i]) );
    	misfit2 = misfit2 + ( (obs2[i] - mod2[i]) * (obs2[i] - mod2[i]) );
    	misfit3 = misfit3 + ( (obs3[i] - mod3[i]) * (obs3[i] - mod3[i]) );
    }

    // Mean of the squared differences
    double misfit[3];
    misfit[0] = misfit1/10;
    misfit[1] = misfit2/10;
    misfit[2] = misfit3/10;

    // Square root of the mean squared differences
    misfit[0] = sqrt(misfit[0]);
    misfit[1] = sqrt(misfit[1]);
    misfit[2] = sqrt(misfit[2]);

    printf("RMSE misfits = %f %f %f\n", misfit[0], misfit[1], misfit[2]);

    FILE *f = fopen("misfit_output.txt", "w");

	if (SINGLE_MISFIT) {
		double singleMis;
		singleMis = (misfit[0]+misfit[1]+misfit[2])/3;
		fprintf(f, "%f\n", singleMis);
		printf("Single RMSE misfit = %f\n", singleMis);
	}
	else {
    	for ( i=0; i<3; i++ ) {
        	fprintf(f, "%f\n", misfit[i]);
    	}
	}
	
    fclose(f);
    
    printf("\n");

}
