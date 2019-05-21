#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{

	// OBSERVED TRACERS
	
	// Noise to add (maximum random noise as percentage of mean)
	double addNoise = 0;

	// Actual best parameters
	double param1 = 175.00;
	double param2 = 2.50;
	
	printf("\nActual Best Parameters = %f %f\n\n", param1, param2);

    // Create Observed tracers
    double obs1[10], obs2[10], obs3[10];
    double sum1=0, sum2=0, sum3=0;
    double x[10];
    double xcount;

    int i;

	// Obs 1
    xcount = 0.1;
    printf("Observations without noise \nObs1 = ");
    for ( i=0; i<10; i++ ) {
        x[i] = xcount;
    	obs1[i] = param2 * x[i] + param1;
        xcount = xcount + 0.1;
        sum1 = sum1 + obs1[i];
        printf("%f ", obs1[i]);
    }
    printf("\n");

	// Obs 2
	printf("Obs2 = ");
    for ( i=0; i<10; i++ ) {
    	obs2[i] = (param1/10) * sin(x[i]) * param2;
    	sum2 = sum2 + obs2[i];
    	printf("%f ", obs2[i]);
    }
    printf("\n");

	//Obs 3
	printf("Obs3 = ");
    for ( i=0; i<10; i++ ) {
    	obs3[i] = (param1/15) * cos(x[i]) * param2;
    	sum3 = sum3 + obs3[i];
    	printf("%f ", obs3[i]);
    }
    printf("\n");
    
    // Add noise
    if (addNoise == 0) {
    	printf("No noise added \n");
    }
    else {
    	srand(16666);
    	double percentNoise = addNoise*100;
    	double rand1, rand2, rand3, obsbefore, noise;
    	
    	printf("Percentage of noise to add = %f \nNse1 = ", percentNoise);
    	for ( i=0; i<10; i++ ) {
    		rand1 = rand() / (RAND_MAX + 1.0);
    		obsbefore = obs1[i];
    		obs1[i] = obs1[i] + ((rand1 * addNoise) * (sum1/10));
    		noise = obs1[i] - obsbefore;
    		printf("%f ", noise);
    	}
    	printf("\n");
    
    	printf("Nse2 = ");
    	for ( i=0; i<10; i++ ) {
    		rand2 = rand() / (RAND_MAX + 1.0);
    		obsbefore = obs2[i];
    		obs2[i] = obs2[i] + ((rand2 * addNoise) * (sum2/10));
    		noise = obs2[i] - obsbefore;
    		printf("%f ", noise);
    	}
    	printf("\n");
    
    	printf("Nse3 = ");
    	for ( i=0; i<10; i++ ) {
    		rand3 = rand() / (RAND_MAX + 1.0);
    		obsbefore = obs3[i];
    		obs3[i] = obs3[i] + ((rand3 * addNoise) * (sum3/10));
    		noise = obs3[i] - obsbefore;
    		printf("%f ", noise);
    	}
    	printf("\n");
    }
    
    // Write observed tracers
    
    // Obs 1
    FILE *f1 = fopen("obs1.txt", "w");
    printf("Observations with noise \nObs1 = ");
    for ( i=0; i<10; i++ ) {
        fprintf(f1, "%f\n", obs1[i]);
        printf("%f ", obs1[i]);
    }
    fclose(f1);
    printf("\n");
    
    // Obs 2
    FILE *f2 = fopen("obs2.txt", "w");
    printf("Obs2 = ");
    for ( i=0; i<10; i++ ) {
        fprintf(f2, "%f\n", obs2[i]);
        printf("%f ", obs2[i]);
    }
    fclose(f2);
    printf("\n");
    
    // Obs 3
    FILE *f3 = fopen("obs3.txt", "w");
    printf("Obs3 = ");
    for ( i=0; i<10; i++ ) {
        fprintf(f3, "%f\n", obs3[i]);
        printf("%f ", obs3[i]);
    }
    fclose(f3);
    printf("\n");
    
    // x
    FILE *f4 = fopen("x.txt", "w");
    for ( i=0; i<10; i++ ) {
        fprintf(f4, "%f\n", x[i]);
    }
    fclose(f4);
    
}