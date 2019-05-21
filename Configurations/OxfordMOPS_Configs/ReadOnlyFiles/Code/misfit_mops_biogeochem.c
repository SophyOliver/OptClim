/* $Header: /Users/ikriest/CVS/mops/tmm_misfit.c,v 1.1 2015/11/17 14:18:51 ikriest Exp $ */
/* $Name: mops-2_0 $*/
/* Edited by Sophy Oliver 2018 */

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_main.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_timer.h"

#include "mops_biogeochem_misfit_data.h"

#define TR v[0]

#define numPO4regions nRegions[0]
#define numNO3regions nRegions[1]
#define numO2regions nRegions[2]


#undef __FUNCT__
#define __FUNCT__ "iniMisfit"
PetscErrorCode iniMisfit(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v)
{
  PetscErrorCode ierr;
  PetscBool flg;
  PetscScalar zero = 0.0;
  PetscViewer fd;
  char *maskFile[MAXNUMTRACERS], *obsFile[MAXNUMTRACERS];
  PetscInt nTracersObs = 3;
  
/* Add your code here */

   averageCost = PETSC_FALSE;
   
/* IK: added for misfit function - start */
 
/* Type of cost function; so far, only annual averages to PO4, O2 and NO3 available */

  ierr = PetscOptionsGetString(PETSC_NULL,"-misfit_file",misfitFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
  strcpy(misfitFile,"");
  sprintf(misfitFile,"%s","misfit.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"misfit will be written to %s\n",misfitFile);CHKERRQ(ierr);
  ierr = PetscFOpen(PETSC_COMM_WORLD,misfitFile,"w",&misfitf);CHKERRQ(ierr);  

  ierr = PetscOptionsHasName(PETSC_NULL,"-average_cost",&averageCost);CHKERRQ(ierr);
	if (averageCost) {

    ierr = iniStepTimer("cost_", Iter0, &costTimer);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Cost for misfit function will be computed starting at (and including) time step: %d\n", costTimer.startTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Cost for misfit function will be computed over %d time steps\n", costTimer.numTimeSteps);CHKERRQ(ierr);	

/* Read volume of each box for scaling according to region size */

	ierr = VecDuplicate(TR,&boxVol);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"volb.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(boxVol,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Read individual grid point volume file for calculating weighting\n");CHKERRQ(ierr);	

/* Read number of regions for each tracer type */
	nRegions = malloc(nTracersObs*sizeof(PetscInt *));
    ierr = PetscOptionsGetIntArray(PETSC_NULL,"-num_mask_regions",nRegions,&nTracersObs,&flg);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of mask regions with the -num_mask_regions option");
	for (iobs=0; iobs<nTracersObs; iobs++) {
	   ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of regions for observation %d = %d\n",iobs,nRegions[iobs]); CHKERRQ(ierr);	
    }

/* Read vectors of masks for individual tracer types and gridboxes to determine regions */
	for (iobs=0; iobs<nTracersObs; iobs++) {
	  maskFile[iobs] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = nTracersObs;
	ierr = PetscOptionsGetStringArray(PETSC_NULL,"-mask_files",maskFile,&maxValsToRead,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"MUST specify mask files with the -mask_files option");
    if (maxValsToRead != nTracersObs) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of mask file names specified");
    }    

	mask = malloc(nTracersObs*sizeof(Vec **));
	for (iobs=0; iobs<nTracersObs; iobs++) {   
	  ierr = VecDuplicateVecs(TR,nRegions[iobs],&mask[iobs]);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer obs %d: reading mask from file %s\n", iobs,maskFile[iobs]);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,maskFile[iobs],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  for (ireg=0; ireg<nRegions[iobs]; ireg++) {
		ierr = VecLoad(mask[iobs][ireg],fd);CHKERRQ(ierr);
	  }
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	}

/* Create the two weightings (one for every grid point: fracMaskVol, and one for every region: fracTotMaskVol) */
	fracMaskVol = malloc(nTracersObs*sizeof(Vec **));
	totMaskVol = malloc(nTracersObs*sizeof(PetscScalar **));
	fracTotMaskVol = malloc(nTracersObs*sizeof(PetscScalar **));
	
	for (iobs=0; iobs<nTracersObs; iobs++) {   
	  ierr = VecDuplicateVecs(TR,nRegions[iobs],&fracMaskVol[iobs]);CHKERRQ(ierr);
      totMaskVol[iobs]=malloc(nRegions[iobs]*sizeof(PetscScalar *));	  
      fracTotMaskVol[iobs]=malloc(nRegions[iobs]*sizeof(PetscScalar *));	  
      
	  ierr = VecMDot(boxVol,nRegions[iobs],mask[iobs],totMaskVol[iobs]);CHKERRQ(ierr);
	  PetscScalar tmpvol = 0.0;
	  for (ireg=0; ireg<nRegions[iobs]; ireg++) {
	  	ierr = VecPointwiseMult(fracMaskVol[iobs][ireg], boxVol, mask[iobs][ireg]);CHKERRQ(ierr);
		ierr = VecScale(fracMaskVol[iobs][ireg],1.0/totMaskVol[iobs][ireg]);CHKERRQ(ierr);
		tmpvol = tmpvol + totMaskVol[iobs][ireg];
	  }

	  for (ireg=0; ireg<nRegions[iobs]; ireg++) {
	    fracTotMaskVol[iobs][ireg] = totMaskVol[iobs][ireg]/tmpvol;
	  }
	  
	}

/* Read the observations */ 

	for (iobs=0; iobs<nTracersObs; iobs++) {
	  obsFile[iobs] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = nTracersObs;
	ierr = PetscOptionsGetStringArray(PETSC_NULL,"-obs_files",obsFile,&maxValsToRead,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"MUST specify observational data files with the -obs_files option");
    if (maxValsToRead != nTracersObs) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Number observational files given = %d\n", maxValsToRead);CHKERRQ(ierr);	   
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of observational data file names specified");
    }
    
    obsbgc = malloc(nTracersObs*sizeof(Vec *));
	for (iobs=0; iobs<nTracersObs; iobs++) {   
	  ierr = VecDuplicate(TR,&obsbgc[iobs]);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer obs %d: reading observations from file %s\n", iobs,obsFile[iobs]);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,obsFile[iobs],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = VecLoad(obsbgc[iobs],fd);CHKERRQ(ierr);
	}
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

/* Initialize the cost function vectors */

	modbgc = malloc(nTracersObs*sizeof(Vec *));
	avmodbgc = malloc(nTracersObs*sizeof(Vec *));
	costbgc = malloc(nTracersObs*sizeof(Vec *));
	for (iobs=0; iobs<nTracersObs; iobs++) {
		ierr = VecDuplicate(TR,&modbgc[iobs]);CHKERRQ(ierr);
		ierr = VecSet(modbgc[iobs],zero);CHKERRQ(ierr);
		ierr = VecDuplicate(TR,&avmodbgc[iobs]);CHKERRQ(ierr);
		ierr = VecSet(avmodbgc[iobs],zero);CHKERRQ(ierr);
		ierr = VecDuplicate(TR,&costbgc[iobs]);CHKERRQ(ierr);
		ierr = VecSet(costbgc[iobs],zero);CHKERRQ(ierr);
	}
	Gavembgc = malloc(nTracersObs*sizeof(PetscScalar **));
	costbgcw1 = malloc(nTracersObs*sizeof(PetscScalar **));
	obsbgcw1 = malloc(nTracersObs*sizeof(PetscScalar **));
	costbgcw2 = malloc(nTracersObs*sizeof(PetscScalar **));
	obsbgcw2 = malloc(nTracersObs*sizeof(PetscScalar **));
	costnorm = malloc(nTracersObs*sizeof(PetscScalar **));

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"cbgc1.petsc",FILE_MODE_WRITE,&fdcbgcavg[0]);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"cbgc2.petsc",FILE_MODE_WRITE,&fdcbgcavg[1]);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"cbgc3.petsc",FILE_MODE_WRITE,&fdcbgcavg[2]);CHKERRQ(ierr);
  
  /* variables specific to externam_forcing_mops_biogeochem.c */
  ierr = VecGetArray(modbgc[0],&localmbgc1);CHKERRQ(ierr);
  ierr = VecGetArray(modbgc[1],&localmbgc2);CHKERRQ(ierr);
  ierr = VecGetArray(modbgc[2],&localmbgc3);CHKERRQ(ierr);

//    costCount=0;

        } else {

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Only average cost function available, currently.\n");CHKERRQ(ierr);	        
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Specify -average_cost, and start time step and number of time steps to average.\n");CHKERRQ(ierr);	        
	ierr = PetscPrintf(PETSC_COMM_WORLD,"No misfit will be computed. \n");CHKERRQ(ierr);	        
        
        }
 /* IK: added for misfit function - end */
 
  return 0;
}

/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "calcMisfit"
PetscErrorCode calcMisfit(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)

{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;
  PetscScalar one = 1.0;  
  PetscInt nTracersObs = 3;
  
  if (Iter0+iLoop>=costTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */	

/* Add your code here */

   if (averageCost) {   /* only option, currently */

/* IK : Add all tracer snapshots, and increase counter by one */ 
 
    if (costTimer.count<=costTimer.numTimeSteps) {
    	for (iobs=0; iobs<nTracersObs; iobs++) {
			ierr = VecAXPY(avmodbgc[iobs],one,modbgc[iobs]);CHKERRQ(ierr);
		}
		costTimer.count++; // costCount+1;
	}
   }
 }
 
 return 0;
 }

/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/


#undef __FUNCT__
#define __FUNCT__ "writeMisfit"
PetscErrorCode writeMisfit(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;
  PetscScalar minusone = -1.0, globalMisfit = 0;
  PetscInt nTracersObs = 3;  
  
  if (Iter0+iLoop>=costTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */	

/* Add your code here */

  if (averageCost){   /* only option, currently */

/* IK : added for misfit function : start */

	  if (costTimer.count==costTimer.numTimeSteps) { /* time to compute cost function and write to file */
	  
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Computing cost function of time average over %d steps at time %10.5f, step %d\n", costTimer.count, tc, Iter0+iLoop);CHKERRQ(ierr);                      

/* IK : average model over a year; this will be stored again in mbgc1avg, ... */

		for (iobs=0; iobs<nTracersObs; iobs++) {
		
			/* average model over a year */
			ierr = VecScale(avmodbgc[iobs],1.0/costTimer.count);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d : Computed annual average concentrations \n", iobs);CHKERRQ(ierr);
			
			/* write vectors of cost functions for a posteriori analysis */
			ierr = VecView(avmodbgc[iobs],fdcbgcavg[iobs]);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d : Wrote model annual average concentrations to cbgc%d.petsc\n", iobs, iobs);CHKERRQ(ierr);  
			
			/* get regional average model tracers to display later, into Gavembgc */
			Gavembgc[iobs]=malloc(nRegions[iobs]*sizeof(PetscScalar *));	  
			ierr = VecMDot(avmodbgc[iobs], nRegions[iobs],fracMaskVol[iobs],Gavembgc[iobs]);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d : Computed regional mean model tracers \n", iobs);CHKERRQ(ierr);
			
			/* copy average model concentrations for next calculations */
			ierr = VecCopy(avmodbgc[iobs],costbgc[iobs]);CHKERRQ(ierr);
			
			/* subtract observations from average model concentrations */
			ierr = VecAXPY(costbgc[iobs],minusone,obsbgc[iobs]);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d : Subtracted model-observations  \n", iobs);CHKERRQ(ierr);
			
			/* square the difference between observations and average model concentrations */
			ierr = VecPow(costbgc[iobs],2.0);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d : Squared the difference  \n", iobs);CHKERRQ(ierr);
			
			/* sum all local misfits and observations, and scale with grid point fractional volume at the same time */
			costbgcw1[iobs]=malloc(nRegions[iobs]*sizeof(PetscScalar *));	  
			obsbgcw1[iobs]=malloc(nRegions[iobs]*sizeof(PetscScalar *));	  
			ierr = VecMDot(costbgc[iobs],nRegions[iobs],fracMaskVol[iobs],costbgcw1[iobs]);CHKERRQ(ierr);
			ierr = VecMDot(obsbgc[iobs],nRegions[iobs],fracMaskVol[iobs],obsbgcw1[iobs]);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d : Summed local misfits and scaled with grid point fractional volume  \n", iobs);CHKERRQ(ierr);
			
			/* scale regional summed misfits and observations by misfit region size weighting */
			costbgcw2[iobs]=malloc(nRegions[iobs]*sizeof(PetscScalar *));	  
			obsbgcw2[iobs]=malloc(nRegions[iobs]*sizeof(PetscScalar *));
			costnorm[iobs]=malloc(nRegions[iobs]*sizeof(PetscScalar *));	 
			
			for (ireg=0; ireg<nRegions[iobs]; ireg++) {
				costbgcw2[iobs][ireg] = costbgcw1[iobs][ireg]*fracTotMaskVol[iobs][ireg];
				obsbgcw2[iobs][ireg] = obsbgcw1[iobs][ireg]*fracTotMaskVol[iobs][ireg];
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d region %d : Scaled regional summed misfits by misfit region size weighting  \n", iobs, ireg);CHKERRQ(ierr);
				
				/* Square root costs */
				costbgcw2[iobs][ireg] = sqrt(costbgcw2[iobs][ireg]);
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d region %d : Square rooted the misfits  \n", iobs, ireg);CHKERRQ(ierr); 
				
				/* Normalise costs using observations */
				costnorm[iobs][ireg] = costbgcw2[iobs][ireg]/obsbgcw2[iobs][ireg];
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Observational tracer %d region %d : Normalised the misfits using observations\n", iobs, ireg);CHKERRQ(ierr);
				 
			}
		}

/* some diagnostic output for checking */

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Check translation\n");CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing cost functions at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);     
		
		for (iobs=0; iobs<nTracersObs; iobs++) {
			
			ierr = PetscPrintf(PETSC_COMM_WORLD,"\nObservational tracer %d : regional total volume fractions (should add up to 1) are:\n", iobs);CHKERRQ(ierr);
			for (ireg=0; ireg<nRegions[iobs]; ireg++) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f\n", fracTotMaskVol[iobs][ireg]);CHKERRQ(ierr);
			}
			ierr = PetscPrintf(PETSC_COMM_WORLD,"\nObservational tracer %d : model averages before region size weighting are:\n", iobs);CHKERRQ(ierr);
			for (ireg=0; ireg<nRegions[iobs]; ireg++) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f\n", Gavembgc[iobs][ireg]);CHKERRQ(ierr);
			}
			ierr = PetscPrintf(PETSC_COMM_WORLD,"\nObservational tracer %d : observation averages before region size weighting are:\n", iobs);CHKERRQ(ierr);
			for (ireg=0; ireg<nRegions[iobs]; ireg++) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f\n", obsbgcw1[iobs][ireg]);CHKERRQ(ierr);
			}
			ierr = PetscPrintf(PETSC_COMM_WORLD,"\nObservational tracer %d : regional misfits after first weighting are:\n", iobs);CHKERRQ(ierr);
			for (ireg=0; ireg<nRegions[iobs]; ireg++) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f\n", costbgcw1[iobs][ireg]);CHKERRQ(ierr);
			}
			ierr = PetscPrintf(PETSC_COMM_WORLD,"\nObservational tracer %d : square rooted non-normalised regional misfits after second weighting are:\n", iobs);CHKERRQ(ierr);
			for (ireg=0; ireg<nRegions[iobs]; ireg++) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f\n", costbgcw2[iobs][ireg]);CHKERRQ(ierr);
			}
			ierr = PetscPrintf(PETSC_COMM_WORLD,"\nObservational tracer %d : normalised final regional misfits are:\n", iobs);CHKERRQ(ierr);
			for (ireg=0; ireg<nRegions[iobs]; ireg++) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f\n", costnorm[iobs][ireg]);CHKERRQ(ierr);
			}
		
/* Write costs to file */

			ierr = PetscPrintf(PETSC_COMM_WORLD,"\nObservational tracer %d : writing to misfit file: %s\n", iobs, misfitFile);CHKERRQ(ierr);
			for (ireg=0; ireg<nRegions[iobs]; ireg++) {
				ierr = PetscFPrintf(PETSC_COMM_WORLD,misfitf,"%.15f\n ", costnorm[iobs][ireg]);CHKERRQ(ierr);
				globalMisfit = globalMisfit + (costnorm[iobs][ireg] * costnorm[iobs][ireg]);
			}
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGlobal Misfit (Sum of all squared misfits): %f\n", globalMisfit);CHKERRQ(ierr);

        ierr = updateStepTimer("cost_", Iter0+iLoop, &costTimer);CHKERRQ(ierr);
//                 costCount = 0;

      } /* costCount==costNumTimeSteps */ 
      
   } else { /* no other option that average cost, currently */ 
   
     ierr = PetscPrintf(PETSC_COMM_WORLD,"No cost function available \n");CHKERRQ(ierr);                      
   }

} /* Iter0+iLoop>=costTimer.startTimeStep */

/* IK : added for misfit function : end */

  return 0;
}



/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "finalizeMisfit"
PetscErrorCode finalizeMisfit(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscErrorCode ierr;
  PetscInt nTracersObs = 3;

/* Add your code here */

/* IK : added for misfit function */

	for (iobs=0; iobs<nTracersObs; iobs++) {
		ierr = VecDestroy(&modbgc[iobs]);CHKERRQ(ierr);
		ierr = VecDestroy(&avmodbgc[iobs]);CHKERRQ(ierr);
		ierr = VecDestroy(&obsbgc[iobs]);CHKERRQ(ierr);
		ierr = VecDestroy(&costbgc[iobs]);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&fdcbgcavg[iobs]);CHKERRQ(ierr);	
	}
	ierr = PetscFClose(PETSC_COMM_WORLD,misfitf);CHKERRQ(ierr);

/* IK : added for misfit function */

  return 0;
}

