#!/usr/bin/env python
"""
Run DFO-LS or gauss-Newton optimisation using HadCM3 on Eddie.
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
import re
import copy
from OptClimVn2 import ModelSimulation, Optimise, StudyConfig,  Submit, HadCM3
import logging
import dfols # optimisation ...

global MODELRUN
# information on Models -- used to store information to decide what runs been done and what need to be done.
# see optFunction for its actual use.
MODELRUN={'modelRuns':collections.OrderedDict(),'modelsToMake':[]}  # create it with empty dicts & lists

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

    cost=config.cost(np.sqrt(pd.Series(c,index=index)/nObs))
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
    fpFmt='%.4g' # format for floating point numbers.
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
def optFunction(params):
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
    obs=modelConfig.getObs() # get obs from the modelConfig
    # need to scale obs then project onto eigenvalues.
    obs=pd.Series(obs) # TODO make getObs return pandas series.
    obs=obs*MODELRUN['scales'] # scale it
    resid=MODELRUN['target'] - obs # make it a residual
    resid = resid.loc[MODELRUN['obsNames']] # order it consistently

    if np.any(pd.isnull(resid)):
        raise Exception("Got missing values  obs")

    resid=resid.values # downgrade to numpy
    resid=MODELRUN['cov_diag_trans'].dot(resid) # transform into orthogonal basis of covar matrix
    # and divide by root nobs for comparability with earlier code.
    #resid=resid/np.sqrt(len(obsNames))


    return  resid # return the observations

def nextName(base='aa',start =0 ):
    """
    generate the next name -- this is a generator function.
    :param base -- base characters for name
    :param prevRuns -- list of previois runs
    :return:
    """
    # initialisation
    maxLen=5
    lenBase=len(base)
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
# set up command line args
modelFn=ModelSimulation # function for model cretion/read.
modelFn=HadCM3.HadCM3 # function for model creation/read.
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

obsNames=config.obsNames(add_constraint=True)
MODELRUN['obsNames']= obsNames
MODELRUN['target']=config.targets(obsNames=obsNames,scale=True)
MODELRUN['scales']=config.scales() # get the scales.
MODELRUN['scales'][config.constraintName()]=1.0
# TODO consider changing all above so MODELRUN is an object with some sensible methods.
refDir=config.referenceConfig() # get the reference configuration.
MODELRUN['refDir']=refDir
baseRunID=config.baseRunID # get the baseRunID.
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



modelDirs=next(os.walk(rootDir))[1] # all sub-directories
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
        if m.getObs() is None or None in m.getObs().values():
                print "*** Model in %s has no observations so ignoring ****"%dir
        else:
            MODELRUN['modelRuns'][key]=m # store the model in modelRuns indexed by its key.

   except IOError:
        pass
nameGen = nextName(base=baseRunID,start = len(MODELRUN['modelRuns']))  # start the name generator.
optimise = config.optimise() # get optimisation info
cov = config.Covariances(trace=verbose, scale=True)  # get covariances. Scaling done for compatability with optFunction.

errCov=cov['CovTotal']
intCov=cov['CovIntVar']
# now add in mu.
nObs=len(obsNames)
mu=optimise.get('mu')
if  mu is not None:
    constraintName=config.constraintName()
    errCov.loc[constraintName,:]=0.0
    errCov.loc[:,constraintName]=0.0
    errCov.loc[constraintName,constraintName]=1.0/(2*mu)
    errCov=errCov.loc[obsNames,obsNames]
    # TODO figure out what to do with internal var -- think we ignore mu so we set it to 0 but that means it gets
    # enormous wt in the inversion for noise distances.. which is not what we want. Could give it a very large values
    # then it doesn't figure in the distance comp...Note that intvar is only used in *termination* criteria. I think for
    # the moment I will got with the large value = 10^9 * max value in array.
    constraintName=config.constraintName()
    intCov.loc[constraintName,:]=0.0
    intCov.loc[:,constraintName]=0.0
    intCov.loc[constraintName,constraintName]=1e9*np.abs(np.max(intCov.values))
    intCov=intCov.loc[obsNames,obsNames] # last stage may not be necessary but tries to get data in right order.


# compute eigenvector and eigenvalues of covariances so we can transform residual into diagonal space.
evalue,evect= np.linalg.eigh(errCov)
transMatrix=(np.diag(evalue**(-0.5)).dot(evect.T))# what we need to do to transform to
MODELRUN['cov_diag_trans']=transMatrix
# transform intCov into eigenvector space of total error.
intCov=transMatrix.dot(intCov).dot(transMatrix.T) # note we now "downgrade" from pandas dataframe to numpy array.
start = config.beginParam(paramNames=varParamNames)


##############################################################
# Main block of code
#  Now we actually run the optimisation code.
##############################################################

np.seterr('raise') # allow fp errors to raise an exception -- this is necessary for optimisation to work.
# algorithm is that if a fp error is triggered then there should be new models to submit.
# common stuff across all algorithms


np.random.seed(123456) # init RNG

algorithmName = optimise['algorithm'].upper()
try:
    if  algorithmName == 'DFOLS':
        prange = (config.paramRanges(paramNames=varParamNames).loc['minParam', :].values,
                  config.paramRanges(paramNames=varParamNames).loc['maxParam', :].values)
        userParams = {'logging.save_diagnostic_info': True, 'logging.save_xk': True, "noise.quit_on_noise_level": True,
                      'init.run_in_parallel': True, 'general.check_objfun_for_overflow': False}
        # merge in values from namedSetting into userParams
        namedSettings = optimise.get('namedSettings', {})
        for k in namedSettings.keys():  # loop over keys
            if not re.search('_comment\s*$', k):  # not a comment
                userParams[k] = namedSettings[k]
        solution= dfols.solve(optFunction,start.values,objfun_has_noise=True, bounds=prange,
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
        print("DFOLS completed: Solution status: %s"%(solution.msg))

        # need to wrap best soln.
        best=pd.Series(solution.x,index=varParamNames)
        config.optimumParams(**(best.to_dict())) # write the optimum params
        solution.diagnostic_info.index=range(0,solution.diagnostic_info.shape[0]) # this does not include "burn in" evals
        info=config.DFOLSinfo(solution.diagnostic_info)
    elif algorithmName == 'GAUSSNEWTON':
        optimise['sigma']=False # wrapped optimisation into cost function.
                # this means need to transform internal var matrix
        optimise['deterministicPerturb']=True # deterministic perturbations.
        nObs=len(MODELRUN['target'])
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
        raise Exception("Don't know what to do with Algorithm: %s")


    cost, simObs, params, dirs, bestEval = storeEvalInfo(config,errCov=errCov)
    # print out some useful information
    start=config.beginParam() ; start.name='start'
    best=config.optimumParams() ; best.name='best'
    print("start to best params are: ")
    print(pd.DataFrame([start,best]))

#except (FloatingPointError,TypeError,np.linalg.LinAlgError, ValueError): #  error whch triggers need to run more models.
except (np.linalg.LinAlgError, TypeError,ValueError): #  error whch triggers need to run more models.
    np.seterr('warn')
    # update config info and save to temp file.
    cost, simObs, params,dirs, bestEval = storeEvalInfo(config, errCov=errCov)

    submitModels=[] # models to be submitted.
    for p in MODELRUN['modelsToMake']: # loop over models to create
            param=pd.Series(p,index=MODELRUN['paramNames']) # parameter values we want wrapped into  Series.
            name=nameGen.next()
            createDir=os.path.join(rootDir,name)

            obsNames=config.obsNames()
            constraintName=config.constraintName()
            if constraintName is not None: obsNames.append(constraintName) # add constraint.
            pDict=param.to_dict() # convert the dataframe to a dict.
            pDict['RUNID']=name
            pDict.update(fixParams) # add any fixed values in.
            # HadCM3 model creation uses pythons facility to deal with keyword arguments
            # ModelSimulation behaves differently which makes life more tricky.
            # TODO unify the approach.
            if genNew: # asked to make new runs?
                model= modelFn(createDir,obsNames=obsNames, create=True,
                               name=name, runTime=runTime, runCode=runCode,
                               ppExePath=config.postProcessScript(), refDirPath=refDir,
                               ppOutputFile=config.postProcessOutput(),**pDict)
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
            status=Submit.eddieSubmit(submitModels,config,rootDir,verbose=verbose,
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
# optionally produce monitoring picture.
if args.monitor:
    import matplotlib.pyplot as plt
    fig=plt.figure('monitor',figsize=[11.7,8.3])
    fig.clear()
    ax=fig.add_subplot(121)
    config.cost().plot(ax=ax,marker='s')
    ax.set_title("Cost")
    ax.margins(0.1)

    ax=fig.add_subplot(122)
    config.parameters(normalise=True).plot(ax=ax,marker='s')
    ax.set_title('Normalised Paramter Values')
    ax.set_ylim(-0.1,1.1)
    ax.margins(0.1)

    fig.tight_layout()
    fig.show()
    fig.savefig(monitorFile) # save the figure
