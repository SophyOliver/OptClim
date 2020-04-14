#!/usr/bin/env python
"""
Run DFO-LS or gauss-Newton optimisation or BOBYQA using MOPS on a local mac computer.
This relies on algorithm being deterministic. So repeatedly runs generating model simulations
as necessary. Only makes sense as function evaluatiosn (running the model) very expensive. So
rereading them is really cheep!

Should be relatively easy to change to other model, optimisation algorithm and  other cluster/super-computers
general approach is that there is a class for specific model inheriting from a base class.
Main work needed is to define parameters and how they map to namelists (or what ever is needed to change model)
Algorithm runs and gets None if the model has not run. This error then gets trapped and
required models get generated.
xxxxSubmit handles the actual submission of the jobs, post processing to the cluster/super-computer
and submits itself.

Command line args:
jsonFile: path to jsonFile defining configuration.
"""
import argparse # parse command line arguments
import collections # want ordered dict
import exceptions
import os
import shutil
import stat
import sys
import numpy as np
import pandas as pd
import subprocess
import re
import copy
from OptClimSO import ModelSimulation, Optimise, StudyConfig
import Submit, MOPS
import logging
import dfols # optimisation ...
import pybobyqa
from contextlib import contextmanager

global MODELRUN
# information on Models -- used to store information to decide what runs been done and what need to be done.
# see optFunction for its actual use.
MODELRUN={'modelRuns':collections.OrderedDict(),'modelsToMake':[]}  # create it with empty dicts & lists

@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)

def storeEvalInfo(config, errCov):
	"""
	Store information from model evaluations  and put them in the config file
	:param config: The configuration which gets updated from information in MODELRUN
	:param errCov: The error covariance
	:return: cost, simulated observations and parameters in that order.
	"""


	plist=[] ;   olist=[] ;   c=[]
	dirNames=[] ; dirs=[]
	invCov=np.linalg.inv(errCov)
	for k,m in MODELRUN['modelRuns'].iteritems():

		name = os.path.basename(m.dirPath)
		dirNames.append(name)
		dirs.append(m.dirPath)
		p = m.getParams(params=MODELRUN['paramNames'])
		o = pd.Series(m.getObs(),name=name)*MODELRUN['scales']
		plist.append(p)
		olist.append(o)
		resid=o-MODELRUN['target']
		c.append(resid.dot(invCov).dot(resid))

	index=pd.Index(dirNames,name='Directory') # make index

	#cost=config.cost(np.sqrt(pd.Series(c,index=index)/nObs))
	cost = config.cost(pd.Series(c,index=index))
	params = config.parameters(pd.DataFrame(plist,index=index))
	simObs=config.simObs(simObs=pd.DataFrame(olist,index=index))
	# work out best directory.
	paramNames=MODELRUN['paramNames']

	best=config.optimumParams(paramNames=paramNames) # get best parameters
	if best is not None:
		# got a best case. Work out its key
		k = genKey(paramValues=best.values,paramNames=paramNames, fixParams=MODELRUN['fixParams'], refDir=config.referenceConfig())
		dirPath = MODELRUN['modelRuns'].get(k).dirPath # and (full path)
		bestEval = os.path.basename(dirPath)
		config.setv('bestEval',bestEval) # best evaluation

	else:
		bestEval=None
	return (cost, simObs, params, dirs, bestEval)

def strSeries(Series):
	"""
	Utility function to convert a Pandas Series to a string for pretty printing
	"""
	import re
	return re.sub( '\s+', ' ', Series.to_string(float_format='%6.3e').replace('\n',' ')).strip()


def genKey(paramNames, paramValues, fixParams, refDir):
	"""
	Generate key from paramNames, paramValuees, fixed parameters and reference dir.
	This should be unique (to some rounding on float parameters)
	:param paramNames -- names of parameters
	:param paramValues -- values for parameters -- done this way because optimisation does not know what parameter names are.
	:param fixParams -- an orderedDict of the fixed parameters index by param name
	:param refDir -- name of reference directory. 
	:return: a tuple as an index.
	"""

	key=[]
	fpFmt='%.34g' # format for floating point numbers.
	# deal with variable parameters -- produced by optimisation so have names and values.
	for k,v in zip(paramNames,paramValues):
		key.append(k)
		if isinstance(v,float):
			key.append((fpFmt%v)) # float point number so use rounded to 4 sig fig version.
		else: # just append the value.
			key.append(repr(v)) # use the object repr method.
	# now do the same with the fixed parameters.
	for k,v in fixParams.iteritems():
		key.append(k)
		if isinstance(v, float):
			key.append((fpFmt % v))  # float point number so use rounded to 4 sig fig version.
		else:  # just append the value.
			key.append(repr(v))  # use the object repr method.
	key.extend(['refDir',refDir])
	key=tuple(key) # convert to tupple
	return  key





# function for optimisation.
def optFunctionDFOLS(params):
	"""
	Function used for optimisation. Returns values from cache if already got it.
	If not got it return array of Nan  and set value to None

	:param: params -- a numpy array with the parameter values. TODO: convince myself that ordering is fixed.
	"""

	global MODELRUN # info on models that we know about and other things needed.

	paramNames= MODELRUN['paramNames'] # save some typing.
	if len(params) != len(paramNames):
		raise Exception("No of parameters not consistent with no of varying parameters")

	key=genKey(paramValues=params, paramNames=paramNames, fixParams=MODELRUN['fixParams'], refDir=MODELRUN['refDir'])
	#print key
	modelConfig = MODELRUN['modelRuns'].get(key, None)  # get the model config or None.
	#print modelConfig
	if modelConfig is None:  # set up models to make.
		MODELRUN['modelsToMake'].append(params) # add params to  list of models to make when we get round to it.
		#print "Need to make ", key
		#return np.repeat(np.nan,len(MODELRUN['obsNames'])) # return array of nans.
		return None # return nothing which should trigger a FP error.

	#print key
	#print "WE GOT TO NON-NONE MODELCONFIG **************************"
	obs=modelConfig.getObs() # get obs from the modelConfig
	#print "obs=modelConfig.getObs() = ", obs
	# need to scale obs then project onto eigenvalues.
	obs=pd.Series(obs) # TODO make getObs return pandas series.
	#print "obs=pd.Series(obs) = ", obs
	#obs=obs*MODELRUN['scales'] # scale it
	#print "obs=obs*MODELRUN['scales'] = ", obs
	#resid=MODELRUN['target'] - obs # make it a residual
	#print "resid=MODELRUN['target'] - obs = ", resid
	resid=obs # in this case we have no scalings and the obs values are already residuals relative to a target
	resid = resid.loc[MODELRUN['obsNames']] # order it consistently
	#print "resid = resid.loc[MODELRUN['obsNames']] = ", resid

	if np.any(pd.isnull(resid)):
		raise Exception("Got missing values  obs")

	resid=resid.values # downgrade to numpy
	#print "resid=resid.values = ", resid
	resid=MODELRUN['cov_diag_trans'].dot(resid) # transform into orthogonal basis of covar matrix
	#print "resid=MODELRUN['cov_diag_trans'].dot(resid) = ", resid
	# and divide by root obs for comparability with earlier code.
	#resid=resid/np.sqrt(len(obsNames))
	#for idx, r in enumerate(resid):
	#	if r < 0:
	#		resid[idx] = resid[idx] * -1
			
	#print "Input to DFOLS = ", resid
	#print "Number of misfits to DFOLS = ", len(resid)
	
	sumSquare = 0
	for r in resid:
		sumSquare = sumSquare + (r*r) # for printing purposes only, as DFO-LS does this internally itself
	#print "Sum-squared cost of all DFOLS inputs values = ", sumSquare


	return  resid # return the observations
	
def optFunctionBOBYQA(params):
	"""
	Function used for optimisation. Returns values from cache if already got it.
	If not got it return array of Nan  and set value to None

	:param: params -- a numpy array with the parameter values. TODO: convince myself that ordering is fixed.
	"""

	global MODELRUN # info on models that we know about and other things needed.

	paramNames= MODELRUN['paramNames'] # save some typing.
	if len(params) != len(paramNames):
		raise Exception("No of parameters not consistent with no of varying paraameters")

	key=genKey(paramValues=params, paramNames=paramNames, fixParams=MODELRUN['fixParams'], refDir=MODELRUN['refDir'])
	modelConfig = MODELRUN['modelRuns'].get(key, None)  # get the model config or None.
	if modelConfig is None:  # set up models to make.
		MODELRUN['modelsToMake'].append(params) # add params to  list of models to make when we get round to it.
		#print "Need to make ", key
		#return np.repeat(np.nan,len(MODELRUN['obsNames'])) # return array of nans.
		return None # return nothing which should trigger a FP error.

	#print key
	#print "WE GOT TO NON-NONE MODELCONFIG **************************"
	obs=modelConfig.getObs() # get obs from the modelConfig
	#print "obs=modelConfig.getObs() = ", obs
	# need to scale obs then project onto eigenvalues.
	obs=pd.Series(obs) # TODO make getObs return pandas series.
	#print "obs=pd.Series(obs) = ", obs
	#obs=obs*MODELRUN['scales'] # scale it
	#print "obs=obs*MODELRUN['scales'] = ", obs
	#resid=MODELRUN['target'] - obs # make it a residual
	#print "resid=MODELRUN['target'] - obs = ", resid
	resid=obs # in this case we have no scalings and the obs values are already residuals relative to a target
	resid = resid.loc[MODELRUN['obsNames']] # order it consistently
	#print "resid = resid.loc[MODELRUN['obsNames']] = ", resid

	if np.any(pd.isnull(resid)):
		raise Exception("Got missing values  obs")

	resid=resid.values # downgrade to numpy
	#print "resid=resid.values = ", resid
	resid=MODELRUN['cov_diag_trans'].dot(resid) # transform into orthogonal basis of covar matrix
	#print "resid=MODELRUN['cov_diag_trans'].dot(resid) = ", resid
	# and divide by root nobs for comparability with earlier code.
	#resid=resid/np.sqrt(len(obsNames))
	#for idx, r in enumerate(resid):
	#	if r < 0:
	#		resid[idx] = resid[idx] * -1
			
	#print "resid = resid * -1 = ", resid
	
	sumSquare = 0
	for r in resid:
		sumSquare = sumSquare + (r*r)
	print "Input to BOBYQA = ", sumSquare


	return  sumSquare # return the observation

def nextName(base='aa',start =0 ):
	"""
	generate the next name -- this is a generator function.
	:param base -- base characters for name
	:param prevRuns -- list of previois runs
	:return:
	"""
	# initialisation
	lenBase=len(base)
	if lenBase>10:  # truncate name length to 10 if too long
		base = base[:10]
		lenBase = len(base)
	maxLen=lenBase+3
	strLen=maxLen-lenBase
	#TODO allow letters a-z as well as numnbers 0-9. That gives 36^3 values ~ 46,000 values for 2 letter base
	# and use lower/upper case to, That gives 62^3 cases or 238,328 cases for 2 letter base or 3844 cases for 3 letter base
	# or ~1000 values for 3 letter base.
	maxValue=10**strLen-1 # can't have value larger than this.
	num = start+1
	fmtStr='%%%d.%dd'%(maxLen,maxLen)
	while True:  # keep going
		str=fmtStr%num
		result=base+str[len(base):]
		result=result.encode('ascii') # convert to ascii
		yield result
		num += 1
		if num > maxValue: raise ValueError("num to large")

## main script
		
finished = "False" # output to bash loop so it knows to stop

# set up command line args
#modelFn=ModelSimulation # function for model cretion/read.
modelFn=MOPS.MOPS # function for model creation/read.
#modelFn=ModelSimulation.ModelSimulation
parser = argparse.ArgumentParser(description="Run study")
parser.add_argument("-dir",help="path to root directory where model runs will be created")
parser.add_argument("jsonFile",help="json file that defines the study")
parser.add_argument("-r","--restart",action='store_true',help="Restart the optimisation by deleting all files except json file")
parser.add_argument("-v","--verbose",action='store_true',help="Provide Verbose information")
parser.add_argument("-n","--noresubmit",action='store_false',help="If set then do not resubmit this script. Good for testing")
parser.add_argument("-d","--dryrun",action='store_true',help="if set do not submit any jobs but do create directories. Good for testing")
parser.add_argument("--nonew",action='store_true',help="If set fail if generate a new model.")
# TODO merge dryrun and noresubmit together.
parser.add_argument("-t","--test",action='store_true',help='If set run test optimisation rather than submitting models.')
parser.add_argument("-m","--monitor",action='store_true',help='Producing monitoring plot after running')
args=parser.parse_args()
verbose=args.verbose
resubmit=args.noresubmit
dryRun=args.dryrun
testRun=args.test
jsonFile=os.path.expanduser(os.path.expandvars(args.jsonFile))
monitor=args.monitor
genNew=not args.nonew

config= StudyConfig.readConfig(filename=jsonFile,ordered=True) # parse the jsonFile.

# setup MODELRUN

obsNames=config.obsNames(add_constraint=False)
MODELRUN['obsNames']= obsNames
MODELRUN['target']=config.targets(obsNames=obsNames,scale=True)
MODELRUN['scales']=config.scales() # get the scales.
constraintName = config.constraintName()
if constraintName is not None:
	MODELRUN['scales'][constraintName]=1.0

	
# TODO consider changing all above so MODELRUN is an object with some sensible methods.
refDir=config.referenceConfig() # get the reference configuration.
MODELRUN['refDir']=refDir
baseRunID=config.baseRunID # get the baseRunID.
runscriptName=config.runscriptName() # get the name of the runscript that runs the 'black box' model
rootDir=os.path.dirname(jsonFile) # get the studyDir
fixParams = config.fixedParams()  # get the fixed parameters and their values.
MODELRUN['fixParams']=fixParams # store fixed parameters for use by optFucntion.

runTime=config.runTime() # run time. None if not present.
runCode=config.runCode() # code to run jobs. None if not present.
restartCMD=[os.path.realpath(sys.argv[0]), args.jsonFile] # generate restart cmd
if args.dir: 
	rootDir=os.path.expanduser(os.path.expandvars(args.dir))# directory defined so overwrite rootDir
	restartCMD.extend(['-dir',args.dir])
if args.monitor:
	restartCMD.extend(['--monitor'])
if verbose:
	print "Running from config %s named %s"%(jsonFile,config.name())
	restartCMD.extend(['--verbose'])
	logging.basicConfig(level=logging.INFO, format='%(message)s')  # detailed info: show every function evaluation


if not resubmit: restartCMD=None # nothing to resubmit.
if args.restart and os.path.isdir(rootDir): # starting anew
	# go and clean all directories and empty _obsLookup by removing everything EXCEPT args.jsonFile
	# algorithm -- iterate over all files in rootDir
	#  if file is a file and not jsonFile delete it,   if file is a dir delete it
	for p in os.listdir(rootDir):
		fp=os.path.join(rootDir,p)
		if os.path.isfile(fp) and os.path.basename(p) !=os.path.basename(jsonFile):
			if verbose: print "Deleting %s"%fp
			os.chmod(fp,stat.S_IWRITE)
			os.remove(fp) # remove the file
		elif os.path.isdir(fp):
			if verbose: print "Deleting %s and contents"%fp
			shutil.rmtree(fp, onerror=ModelSimulation.errorRemoveReadonly)

if not os.path.isdir(rootDir):
	os.mkdir(rootDir) # create dir

optimise = config.optimise() # get optimisation info
algorithmName = optimise['algorithm'].upper()

modelDirs=next(os.walk(rootDir))[1] # all sub-directories
lend = len(modelDirs)
varParamNames=config.paramNames() # extract the parameter names if have them
MODELRUN['paramNames']=varParamNames # store the param names.
startRunCount=0 # the count for the start run
for d in modelDirs:
	try:
		dir=os.path.join(rootDir,d)
		m= modelFn(dir) # read in the model. # what to do if model crash but want to continue anyhow???
		# generate model with key simulatedObservations set to values wanted.
		params=m.getParams(params=varParamNames, verbose=verbose)	
		key=genKey(params.keys(),params.values(),fixParams,refDir) # Generate key. Note order is set my varParamNames.
		if ((m.getObs() is None or None in m.getObs().values())) or (algorithmName != 'CMAES' and lend != 1 and (m.getParams() is None or None in m.getParams().values())) :
			print "*** Model in %s has no observations so ignoring, or if cmaes then will update later ****"%dir
		else:
			MODELRUN['modelRuns'][key]=m # store the model in modelRuns indexed by its key.

	except IOError:
		pass
nameGen = nextName(base=baseRunID,start = len(MODELRUN['modelRuns']))  # start the name generator.
cov = config.Covariances(trace=verbose, scale=True)  # get covariances. Scaling done for compatability with optFunction.

errCov=cov['CovTotal']

# now add in mu.
nObs=len(obsNames)
mu=optimise.get('mu')
if  mu is not None:
	constraintName=config.constraintName()
	errCov.loc[constraintName,:]=0.0
	errCov.loc[:,constraintName]=0.0
	errCov.loc[constraintName,constraintName]=1.0/(2*mu)
	errCov=errCov.loc[obsNames,obsNames]


# compute eigenvector and eigenvalues of covariances so we can transform residual into diagonal space.
evalue,evect= np.linalg.eigh(errCov)
transMatrix=(np.diag(evalue**(-0.5)).dot(evect.T))# what we need to do to transform to
MODELRUN['cov_diag_trans']=transMatrix
# transform intCov into eigenvector space of total error.
start = config.beginParam(paramNames=varParamNames)


##############################################################
# Main block of code
#  Now we actually run the optimisation code.
##############################################################

np.seterr('raise') # allow fp errors to raise an exception -- this is necessary for optimisation to work.
# algorithm is that if a fp error is triggered then there should be new models to submit.
# common stuff across all algorithms


np.random.seed(123456) # init RNG

try:
	if  algorithmName == 'DFOLS':
		if verbose:
			print("running DFOLS *********")
			
		prange = (config.paramRanges(paramNames=varParamNames).loc['minParam', :].values,
					config.paramRanges(paramNames=varParamNames).loc['maxParam', :].values)
		userParams = {'logging.save_diagnostic_info': True, 'logging.save_xk': True, "noise.quit_on_noise_level": False,
					'init.run_in_parallel': True, 'general.check_objfun_for_overflow': False}
		# merge in values from namedSetting into userParams
		namedSettings = optimise.get('DFOLS_namedSettings', {})
		for k in namedSettings.keys():  # loop over keys
			if not re.search('_comment\s*$', k):  # not a comment
				userParams[k] = namedSettings[k]
		solution= dfols.solve(optFunctionDFOLS,start.values,objfun_has_noise=False, bounds=prange,
									maxfun=optimise.get('maxfun',100),
									rhobeg=optimise.get('rhobeg',1e-1),
									rhoend=optimise.get('rhoend',1e-3),
									user_params=userParams, scaling_within_bounds=True)
		if solution.flag == solution.EXIT_LINALG_ERROR: # linear algebra error
			raise np.linalg.LinAlgError # re-raise the linear algebra error which will trigger doing more runs..
		elif solution.flag not in (solution.EXIT_SUCCESS,solution.EXIT_MAXFUN_WARNING):
			print("dfols failed with flag %i error : %s"%(solution.flag,solution.msg))
			raise Exception("Problem with dfols")
		## code here will be run when DFOLS has completed. It mostly is to put stuff in the final JSON file
		## so can easily be looked at for subsequent analysis.
		## some of it could be done even if DFOLS did not complete.
		

		# need to wrap best soln.
		best=pd.Series(solution.x,index=varParamNames)
		config.optimumParams(**(best.to_dict())) # write the optimum params
		solution.diagnostic_info.index=range(0,solution.diagnostic_info.shape[0]) # this does not include "burn in" evals
		info=config.DFOLSinfo(solution.diagnostic_info)
		solution.diagnostic_info.to_csv(os.path.join(rootDir,'diagnostics.txt'))
		print("DFOLS completed: Solution status: %s"%(solution.msg))
		
	elif  algorithmName == 'BOBYQA':
		prange = (config.paramRanges(paramNames=varParamNames).loc['minParam', :].values,
					config.paramRanges(paramNames=varParamNames).loc['maxParam', :].values)
		userParams = {'logging.save_diagnostic_info': True, 'logging.save_xk': True, "noise.quit_on_noise_level": False,
					'init.run_in_parallel': True, 'general.check_objfun_for_overflow': False}
		# merge in values from namedSetting into userParams
		namedSettings = optimise.get('BOBYQA_namedSettings', {})
		for k in namedSettings.keys():  # loop over keys
			if not re.search('_comment\s*$', k):  # not a comment
				userParams[k] = namedSettings[k]
		solution= pybobyqa.solve(optFunctionBOBYQA,start.values,objfun_has_noise=False, bounds=prange,
									maxfun=optimise.get('maxfun',100),
									rhobeg=optimise.get('rhobeg',1e-1),
									rhoend=optimise.get('rhoend',1e-3),
									user_params=userParams, scaling_within_bounds=True)
		if solution.flag == solution.EXIT_LINALG_ERROR: # linear algebra error
			raise np.linalg.LinAlgError # re-raise the linear algebra error which will trigger doing more runs..
		elif solution.flag not in (solution.EXIT_SUCCESS,solution.EXIT_MAXFUN_WARNING):
			print("bobyqa failed with flag %i error : %s"%(solution.flag,solution.msg))
			raise Exception("Problem with bobyqa")
		## code here will be run when PYBOBYQA has completed. It mostly is to put stuff in the final JSON file
		## so can easily be looked at for subsequent analysis.
		## some of it could be done even if DFOLS did not complete.
		print("BOBYQA completed: Solution status: %s"%(solution.msg))

		# need to wrap best soln.
		best=pd.Series(solution.x,index=varParamNames)
		config.optimumParams(**(best.to_dict())) # write the optimum params
		solution.diagnostic_info.index=range(0,solution.diagnostic_info.shape[0]) # this does not include "burn in" evals
		info=config.BOBYQAinfo(solution.diagnostic_info)
		solution.diagnostic_info.to_csv(os.path.join(rootDir,'diagnostics.txt'))
		
	elif  algorithmName == 'CMAES':
		
		with cd(rootDir):
		
			if not os.path.isfile('nIter.txt'):
				# Write nIter.txt using json file
				
				# Get json param bounds
				prange = (config.paramRanges(paramNames=varParamNames).loc['minParam', :].values,
					config.paramRanges(paramNames=varParamNames).loc['maxParam', :].values)
					
				# Get CMAES settings
				namedSettings = optimise.get('CMAES_namedSettings', {})
					
				with open('nIter.txt', "w") as nIterFile:
				
					nIterFile.write("operational parameters: model name, #parameters, #sessions, #iterations, random number seed (integer, 0=default), attractor flag (0=default=""use no attractor""):\n")
					nIterFile.write("%s\n"%config.name()) # experiment name
					nIterFile.write("%d\n"%len(prange[0])) # number of params
					nIterFile.write("%d\n"%namedSettings['populationSize']) # population size
					nIterFile.write("%d\n"%optimise.get('maxfun',100)) # max cmaes iterations
					nIterFile.write("%d\n"%namedSettings['randomNumberSeed']) # random number seed
					nIterFile.write("%d\n"%namedSettings['attractorFlag']) # attractor flag
					nIterFile.write("lower and upper bounds (box constraints) for parameters:\n")
					
					# Write json param bounds
					for idx, lower in enumerate(prange[0]):
						nIterFile.write("%.15f"%lower)
						nIterFile.write(" %.15f\n"%prange[1][idx])
					
					# Write json optional attractors
					nIterFile.write("optional attractor:\n")
					for attractor in namedSettings['attractor']:
						nIterFile.write("%.15f\n"%attractor)
					
					# Iteration log
					nIterFile.write("iteration log: iteration, objective value, termination flag:\n")
					nIterFile.write("0 0 0\n")
			
			# Read nIter.txt last line info
			maxIter = optimise.get('maxfun',100);
			with open('nIter.txt', 'r') as f:
				lines = f.read().splitlines()
				last_line = lines[-1]
			nIter = [int(s) for s in last_line.split() if s.isdigit()][0]
			terminate = int(last_line[-1])
			
			# Check if need to run more models then run CMAES
			if nIter <= maxIter:
				if terminate < 1:
					prevDir = baseRunID + str(nIter).zfill(3)
					cmaesString = '$OPTCLIMTOP/cmaes/./cmaes ' + str(nIter) + ' nIter.txt ' + prevDir + " > serial_" + str(nIter) + ".log"
					subprocess.check_output(cmaesString, shell=True)
					
					# Code to read in cmaes-specific file if it exists. Directs us to the best param file of the cmaes population to update config with
					if nIter > 0:
						m = modelFn(prevDir)
						cmaesFile = os.path.join(prevDir, 'best_cmaes.txt')
						if os.path.isfile(cmaesFile):
							with open(cmaesFile) as f:
								nI = f.readline().strip()
								besti = f.readline().strip()
							parFile=os.path.join(prevDir, 'parameters_input_' + str(nI) + '_' + str(besti) + '.txt')
							if os.path.isfile(parFile):
								pardf = pd.read_csv(parFile, header=None, index_col=False)
								pardf = pardf.values.flatten()
								pardf = pardf.tolist()
								params=pd.Series(pardf,index=MODELRUN['paramNames']) # parameter values we want wrapped into  Series.
								params=params.to_dict(into=collections.OrderedDict)
								m.setParams(params)
								key=genKey(params.keys(),params.values(),fixParams,refDir) # Generate key. Note order is set my varParamNames.
								MODELRUN['modelRuns'][key]=m # store the model in modelRuns indexed by its key.
							else:
								print("File %s not found"%parFile)
						else:
							print("File %s not found"%cmaesFile)
				else:
					solutionMsg = "CMAES termination flag = 1 (True)"
			if nIter < maxIter:
				nameGen = nextName(base=baseRunID,start = len(MODELRUN['modelRuns']))  # start the name generator.
				MODELRUN['modelsToMake'].append('RUN_CMAES')
				raise np.linalg.LinAlgError # re-raise the linear algebra error which will trigger doing more runs..
			else:
				solutionMsg = "The maximum number of iterations has been reached"	
					
			# If got to here then all done
			
			print("CMAES completed: Solution status: %s"%(solutionMsg))

			# Need to wrap best soln.
			# best=pd.Series(solution.x,index=varParamNames)
			# config.optimumParams(**(best.to_dict())) # write the optimum params
			
	elif algorithmName == 'GAUSSNEWTON':
		optimise['sigma']=False # wrapped optimisation into cost function.
		# this means need to transform internal var matrix
		optimise['deterministicPerturb']=True # deterministic perturbations.
		nObs=len(MODELRUN['target'])
		intCov=cov['CovIntVar']
		# now add in mu.
		nObs=len(obsNames)
		mu=optimise.get('mu')
		if  mu is not None:
			# TODO figure out what to do with internal var -- think we ignore mu so we set it to 0 but that means it gets
			# enormous wt in the inversion for noise distances.. which is not what we want. Could give it a very large values
			# then it doesn't figure in the distance comp...Note that intvar is only used in *termination* criteria. I think for
			# the moment I will got with the large value = 10^9 * max value in array.
			constraintName=config.constraintName()
			intCov.loc[constraintName,:]=0.0
			intCov.loc[:,constraintName]=0.0
			intCov.loc[constraintName,constraintName]=1e9*np.abs(np.max(intCov.values))
			intCov=intCov.loc[obsNames,obsNames] # last stage may not be necessary but tries to get data in right order.
		intCov=transMatrix.dot(intCov).dot(transMatrix.T) # note we now "downgrade" from pandas dataframe to numpy array.

		best, status, info = Optimise.gaussNewton(optFunction, start.values,
									config.paramRanges(paramNames=varParamNames).values.T,
									config.steps(paramNames=varParamNames).values,
									np.zeros(nObs), optimise, cov = np.diag(np.ones(nObs)), cov_iv=intCov, trace=verbose)
		jacobian = config.GNjacobian(info['jacobian'])
		hessian = config.GNhessian(info['hessian'])
		params = config.GNparams(info['bestParams'])
		cost = config.GNcost(info['err_constraint'])
		alpha = config.GNalpha(info['alpha'])
		best = pd.Series(best, index=config.paramNames(), name=config.name())  # wrap best result as pandas series
		print "All done with status: %s " % (status)
		# iterate through linesearch printing best param values and error.
		for it in params.Iteration:
			idx = int(it)
			print 'iter: %d %5.2f %5.2f %s' % (
			it, cost[idx], alpha[idx], strSeries(params.sel(Iteration=it).to_pandas()))
		best=config.optimumParams(**(best.to_dict()))  # write the optimum params
	else:
		raise Exception("Don't know what to do with Algorithm: %s "%algorithmName)


	cost, simObs, params, dirs, bestEval = storeEvalInfo(config,errCov=errCov)
	# print out some useful information
	#start=config.beginParam() ; start.name='start'
	#best=config.optimumParams() ; best.name='best'
	#print("start to best params are: ")
	#print(pd.DataFrame([start,best])) 
		
	finished = "True" # output to bash loop so it knows to stop

#except (FloatingPointError,TypeError,np.linalg.LinAlgError, ValueError): #  error whch triggers need to run more models.
except (np.linalg.LinAlgError, TypeError): #  error whch triggers need to run more models.
	print("we got to 532")
	np.seterr('warn')
	
	# If CMAES bypass Submit.py for now and submit it right here - not very user-friendly t the moment though!
	# Assumes the runscript is called runscript_cmaes
	# update config info and save to temp file.
	cost, simObs, params,dirs, bestEval = storeEvalInfo(config, errCov=errCov)

	submitModels=[] # models to be submitted.
	for p in MODELRUN['modelsToMake']: # loop over models to create
			
		name=nameGen.next()
		createDir=os.path.join(rootDir,name)
		obsNames=config.obsNames()
		constraintName=config.constraintName()
		if constraintName is not None: obsNames.append(constraintName) # add constraint.
			
		#if type(p) is str and p == "RUN_CMAES":
			#pDict = dict.fromkeys(MODELRUN['paramNames'])
		#else:
		param=pd.Series(p,index=MODELRUN['paramNames']) # parameter values we want wrapped into  Series.
		print("Param = ", param)
		pDict=param.to_dict(into=collections.OrderedDict) # convert the dataframe to a dict.
				
		pDict['RUNID']=name
		pDict.update(fixParams) # add any fixed values in.
		#print("pDict = ", pDict)
		
		# HadCM3 model creation uses pythons facility to deal with keyword arguments
		# ModelSimulation behaves differently which makes life more tricky.
		# TODO unify the approach.
		if genNew: # asked to make new runs?
			model= modelFn(createDir, runscriptName=runscriptName, obsNames=obsNames, create=True,
							name=name,
							ppExePath=config.postProcessScript(), refDirPath=refDir,
							ppOutputFile=config.postProcessOutput(), parameters = pDict)
			if verbose:
				print("createDir is %s path is %s"%(createDir,os.getcwd()))
			submitModels.append(model) # add to list of models to submit
			modelDirs = modelDirs.append(pd.Series([createDir], index=[name])) # Could do this outside the loop.
		else:
			print "Params are ",pDict
			print' all keys are :'
			for k in MODELRUN['modelRuns'].keys():
				print k

			raise Exception("Generating model when set -nonew")

	if len(submitModels) > 0 and not dryRun: # got some models to submit and dryRun not set?
		if testRun:
			status = Submit.submit(submitModels, config, rootDir, verbose=verbose, resubmit=False) # fake submit the models.
		else: # really submit!
		
			# As OptClimSO has been changed to run everything as one large job, this arcSubmit is actually the same as testRun submit above.
			status=Submit.arcSubmit(submitModels,config,rootDir,verbose=verbose,
									resubmit=restartCMD, runCode=runCode) # submit the models

		# models submitted check all worked.. Note there is a race condition in that the models get submitted before the dict
		# gets updated.. But it is not much  of a race.
		if status != True:
			raise  exceptions.RuntimeError("Submission failed")
		# update the directories in dictFile
		config.directories(modelDirs)
	elif dryRun: # dryRun case.
		print "Dry run so no models submitted. Would have submitted %d"%(len(submitModels))
	else:
		print "No models to submit"
		raise # re-raise the exception as we've triggered an error but have nothing to do...
#all end of run stuff
# work out what final json file is called reagrdless of why run finished.
rootDiagFiles, ext = os.path.splitext(os.path.basename(jsonFile))
finalJsonFile = os.path.join(rootDir, rootDiagFiles + "_final" + ext)
monitorFile= os.path.join(rootDir,rootDiagFiles+"_monitor.png")
config.save(filename=finalJsonFile) # save the (updated) configuration file.
#print finalJsonFile
# optionally produce monitoring picture.
if args.monitor and len(config.cost())!=0:
	import matplotlib.pyplot as plt
	fig=plt.figure('monitor',figsize=[11.7,8.3])
	fig.clear()
	ax=fig.add_subplot(121)
	config.cost().plot(ax=ax,marker='s')
	ax.set_title("Cost (Sum of Squares)")
	ax.margins(0.1)

	ax=fig.add_subplot(122)
	config.parameters(normalise=True).plot(ax=ax,marker='s')
	ax.set_title('Normalised Parameter Values')
	ax.set_ylim(-0.1,1.1)
	ax.margins(0.1)

	fig.tight_layout()
	fig.show()
	fig.savefig(monitorFile) # save the figure
	
print finished
