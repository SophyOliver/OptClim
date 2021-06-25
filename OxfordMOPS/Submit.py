"""
EDITED FOR MOPS ONLY: Functions for model submission. Idea is one function per machine/cluster/os.
Currently defined are:
	1) submit -- which doesn't actually submit anything but runs fakemodel on each model to generate observations
	2) arcSubmit - which submits models (and ancillary scripts) to ARC.
"""
import os
import time
import netCDF4
import subprocess
import sys
import fileinput
import pandas as pd
import numpy as np
from contextlib import contextmanager

@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)

def submitFake(model_list, config, rootDir, verbose=False, resubmit=None, runCode=None, runTime = None):
	"""
	Submit the model -- depends on the machine on which the model is to be run.
	This variant does nothing.. except generates the observations depending on the parameters.

	:param: model_list -- a list of model configurations to submit
	:param: config -- configuration  for studies
	:param: rootDir -- root directory. Not used but available to be parallel to eddieSubmit 
	:param: verbose -- default False; if True provide more info
	:param resubmit -- default None. If not None what gets submitted for next run. For dummy submission this doesn't do much!
	:param runCode -- default None. If not None the project to be used for jobs. For dummy submission this doesn't do much!
	:return: status of submission
	"""

	outputPath=config.postProcessOutput()
	# iterate over list of models to submit
	for model in model_list:
		# put some dummy data in the ouput file

		param=model.getParams() # get the parameters
		modelDir=model.dirPath
		outPath=os.path.join(modelDir,outputPath)
		obs=fakeModel(param,config)
		obs.name=model.name()
		# generate text file...

		with open(outPath, "w") as text_file:
			for k,v in obs.items():
				text_file.write("%f\n"%v)
				
				
def submit(model_list, config, rootDir, verbose=False, resubmit=None, runCode=None, runTime = None):
	"""
	Submit the model -- depends on the machine on which the model is to be run.
	This variant does nothing.. except generates the observations depending on the parameters.

	:param: model_list -- a list of model configurations to submit
	:param: config -- configuration  for studies
	:param: rootDir -- root directory. Not used but available to be parallel to eddieSubmit 
	:param: verbose -- default False; if True provide more info
	:param resubmit -- default None. If not None what gets submitted for next run. For dummy submission this doesn't do much!
	:param runCode -- default None. If not None the project to be used for jobs. For dummy submission this doesn't do much!
	:return: status of submission
	"""

	outputPath=config.postProcessOutput()
	
	# iterate over list of models to submit
	for model in model_list:
		# put some dummy data in the ouput file
		modelSubmitName=model.submit()
		if verbose: print "Submitting ",modelSubmitName
		with cd(model.dirPath):
			subprocess.check_output(modelSubmitName, shell=True) # submit the script


# end of submit
	return True # submission worked!

def fakeModel(paramV, studyCfg, obsNames=None, trace=False):
	"""
	Fake model values --
	Assume that simulated values are a linear sum of effects from each parameter and that each parameter affects two
	types of parameters. With each region behaving similarly.
	:param paramV: parameter vector to run at
	:param studyCfg: study configuration -- which contains much useful information
	:param obsNames: (optional) the observations to make. If not provided studyCfg will be queried for them
	:param trace (optional) If True (default is False) more information will be printed out.
	:return: simulated observations as a pandas object
	"""

	if obsNames is None :
		use_obsNames=studyCfg.obsNames()
	else:
		use_obsNames=obsNames.copy()

	if not isinstance(paramV,pd.Series):
		paramV=pd.Series(paramV)

	paramF=pd.to_numeric(paramV,'coerce') # convert param values to float. Everything else is NaN # alas doesn't work at vn 0.18..
	# this is a bad hack -- we just remove RUNID, RUN_TARGET & START_TIME
	# from the paramV to make paramF
	names=[n for n in paramV.index if n not in ['RUNID', 'RUN_TARGET', 'START_TIME']]
	paramF=pd.to_numeric(paramV[names]) # generate parameters with ones that cause trouble out
	paramNames=paramF.index # index is names
	deterministic = True # hard wire -- set False if want truly random numbers
	# code below is largely cut-n-paste from optimise with name change and use of values
	if deterministic: # initialise RNG based on paramV
		# this is tricky because we need to convert them to integer arrays...
		values=paramF.values
		nonZero= (np.abs(values) > 1e-9)
		int_v=np.int_(values)
		scale = 10.0**np.floor(np.log10(np.abs(values[nonZero]))-3) # allow a thousand different values
		int_v= np.int_(np.floor(np.abs(values[nonZero]/scale)+1)) #now got it as integer
		intSim=int_v.sum()
		if (intSim < 0):
			raise ValueError("intSum stuffed")
		np.random.seed(intSim) # set the seed up
		if trace: print ": Seed set to ",int_v

	standardObs = studyCfg.standardObs(obsNames=use_obsNames,scale=False) # standard obs
	nobs=standardObs.shape[0]
	linTerm=np.array([[0.1]*10,[0.2]*10,[0.3]*10]).flatten()[0:nobs]/studyCfg.scales(obsNames=obsNames)
	sqrTerm=linTerm*2.
	cubeTerm=linTerm*4.
	pwr=pd.DataFrame([standardObs.values,
					linTerm, # linear term
					sqrTerm, # square term
					cubeTerm], # cubic term
					columns=standardObs.index)
	#cov=studyCfg.Covariances(scale=False)
	#noise = cov['CovIntVar'] # noise.

	standardParam = studyCfg.standardParam(paramNames=paramNames) # standard values

	rangep = studyCfg.paramRanges(paramNames=paramNames)
	delta_p = (paramF-standardParam)/rangep.loc['rangeParam',:] # scale parameters
	delta_p=delta_p.dropna() # remove any NaN's
	noise = np.diag(np.repeat(0.1, nobs))
	result = pd.Series(np.random.multivariate_normal(standardObs.loc[use_obsNames].values, noise),index=use_obsNames) # 0.1 here is the noise of the model (change for MOPS)
	# initialise with standard values + noise realisation.
	# hardwire for moment # need to make this a matrix -- converting from fn of parameters to fn of obs,
	# iterate over parameters

	for index,param in enumerate(delta_p.index):
		obsIndex = index%nobs
		obs=use_obsNames[obsIndex] # those are the obs this parameter affects
		##pandas nhx, shx and tropics to them
		for p in range(1,4): # hardwired limit on powers -- up to Cubic.
			result.loc[obs] += pwr.ix[p,obsIndex]*(delta_p[param]**p)

	return result
## end of fakeModel.

def arcSubmit(model_list, config,rootDir, verbose=False,  resubmit=None, runCode=None):
	"""
	Submit models to arc,  and the next iteration in the algorithm.
	:param model_list: a list of model configuration to submit
	:param config: study configuration file
	:param rootDir: root directory where files etc are to be created and found
	:param args -- the arguments that the main script was called with.
	:param resubmit -- default None. If not None next iteration cmd to be submitted. 
		Normally should be the script the user ran so next stage of model running happens.
	:param runCode -- default None, If not none specifies the project code.
	:return: status of submission

	Does the following:
		2) Submits (if provided) resubmit so once the array of post processing jobs has completed the next bit of the algorithm gets ran.
		3) Submits the model simulations

	This algorithm is not particularly robust to failure -- if anything fails the various jobs will be sitting around
	Releasing them will be quite tricky! You can always kill and run again!
	"""
	jobID = []
	for model in model_list:
		# put some dummy data in the ouput file
		modelSubmitName=model.submit()
		if verbose: print "Submitting ",modelSubmitName
		with cd(model.dirPath):
			jID = subprocess.check_output("sbatch -J %s --export=ALL %s" % (model.name(), modelSubmitName), shell=True) # submit the script (change devel after, and shouldn't have to ssh in)
		jobID.append(jID[20:-1])
		
	jobIDstr=':$'.join(jobID) # make single string appropriately formatted of job ids..
	# now re-run this entire script so that the next iteration in the algorithm.
	# can be run
	if resubmit is not None:
		# Submit the next job in the iteration. runOptimise is very quick so no need to submit to ARC again - just run on the front end.
		
		jobName='RE'+config.name()
		# TODO move to better python syntax for var printing. Think can use named vars in...
		cmd = ["sbatch -p devel --export=ALL --time=10 --dependency=afterany:%s -J %s "%(jobIDstr,jobName)]
		cmd.extend(resubmit) # add the arguments in including the programme to run..
		#cmd = resubmit
		cmd=' '.join(cmd) # convert to one string.
		cmd = cmd + " &>progressResubmit.txt"
		if verbose: print "Next iteration cmd is ", cmd
		jid = subprocess.check_output(cmd, shell=True) # submit the script. Good to remove shell=True 
		#subprocess.check_output(cmd, shell=True)
		if verbose: print "Job ID for next iteration is %s"%jid[20:-1]

	return True

def arcSubmit_oneJob(model_list, config,rootDir, verbose=False,  resubmit=None, runCode=None):
	"""
	Run model on the front end.
	:param model_list: a list of model configuration to submit
	:param config: study configuration file
	:param rootDir: root directory where files etc are to be created and found
	:param args -- the arguments that the main script was called with.
	:param resubmit -- default None. If not None next iteration cmd to be submitted. 
		Normally should be the script the user ran so next stage of model running happens.
	:param runCode -- default None, If not none specifies the project code.
	:return: status of submission

	This algorithm is not particularly robust to failure -- if anything fails the various jobs will be sitting around
	Releasing them will be quite tricky! You can always kill and run again!
	"""
	
	#jobID = []
	for model in model_list:
		# put some dummy data in the ouput file
		modelSubmitName=model.submit()
		if verbose: print "Submitting ",modelSubmitName
		with cd(model.dirPath):
			subprocess.check_output(modelSubmitName, shell=True) # submit the script

	return True


def eddieSubmit(model_list, config,rootDir, verbose=False,  resubmit=None, runCode=None):
    """
    Submit models to eddie, the post processing and the next iteration in the algorithm.
    :param model_list: a list of model configuration to submit
    :param config: study configuration file
    :param rootDir: root directory where files etc are to be created and found
    :param args -- the arguments that the main script was called with.
    :param resubmit -- default None. If not None next iteration cmd to be submitted. 
          Normally should be the script the user ran so next stage of model running happens.
        post processing will be submitted the main script will not be submitted. For dummy submission this doesn't do much!
    :param runCode -- default None, If not none specifies the project code.
    :return: status of submission

    Does the following:
        1) Submits the post processing jobs as a task array in held state.
        2) Submits (if provided) resubmit so once the array of post processing jobs has completed the next bit of the algorithm gets ran.
        3) Submits the model simulations -- which once each one has run will release the appropriate post processing task

    This algorithm is not particularly robust to failure -- if anything fails the various jobs will be sitting around
    Releasing them will be quite tricky! You can always kill and run again!
    """
    
    outputDir=os.path.join(rootDir,'jobOutput') # directory where output goes. 
    # try and create it. 
    try: 
        os.makedirs(outputDir)
    except OSError:
        if not os.path.isdir(outputDir):
            raise
    
    sshCmd='ssh login01.ecdf.ed.ac.uk " cd %s ; '%(os.getcwd()) # need to ssh to a login node to do things to Q's and cd to current dir
    #
    modelDirFile=os.path.join(rootDir,'tempDirList.txt') # name of file containing list of directories for post processing stage
    with open(modelDirFile, 'w') as f:
        for m in model_list:
            f.write(m.dirPath+','+m.ppExePath()+','+m.ppOutputFile()+'\n') # write out info for post processing job.
    # submit the following.. Need path to postProcess.sh
    jobName='PP'+config.name()
    ## work out postprocess script path
    postProcess=os.path.expandvars('$OPTCLIMTOP/eddie/postProcess.sh')
    scriptName=os.path.expandvars('$OPTCLIMTOP/eddie/qsub.sh')
    # TODO move to better python syntax for var printing. Think can use named vars in below.
    qsub_cmd='qsub -l h_vmem=2G -l h_rt=00:10:00 -V -cwd -e %s -o %s'%(outputDir,outputDir) # std stuff for submission
    # means        #  2 Gbyte Mem   10 min run, cur env, curr wd, output (error & std) in OutputDir
    # deal with runCode
    if runCode is not None: qsub_cmd += ' -P %s '%(runCode)
    cmd = qsub_cmd+' -t 1:%d -h -N %s '%(len(model_list),jobName)
    cmd += postProcess
    cmd += " %s %s "%(modelDirFile, config.fileName())
    if verbose: print "postProcess task array cmd is ",cmd
    # run the post process and get its job id
    jid = subprocess.check_output(sshCmd+cmd+'"', shell=True)
    #  '"' and shell=True seem necessary. Would be good to avoid both
    postProcessJID=jid.split()[2].split('.')[0] # extract the actual job id.
    if verbose: print "postProcess array job id is %s"%postProcessJID
    # TODO wrap this in a try/except block.
    # write the jobid + N into the model -- for later when 
    #  model gets some processing.
    for indx in range(len(model_list)):
        model_list[indx].jid=postProcessJID+'.%d'%(indx+1)

    # now submit this entire script so that the next iteration in the algorithm.
    # can be run
    if resubmit is not None:
        # submit the next job in the iteration. -hold_jid jid means the post processing job will only run after the
        # arry of post processing jobs has ran.
        jobName='RE'+config.name()
        # TODO move to better python syntax for var printing. Think can use named vars in...
        cmd = [qsub_cmd,'-hold_jid %s -N %s  %s'%(postProcessJID,jobName, scriptName)]
        cmd.extend(resubmit) # add the arguments in including the programme to run..
        cmd=' '.join(cmd) # convert to one string.
        if verbose: print "Next iteration cmd is ", cmd
        jid = subprocess.check_output(sshCmd+cmd+'"', shell=True) # submit the script. Good to remove shell=True and '"'
        jid = jid.split()[2]  # extract the actual job id.
        if verbose: print "Job ID for next iteration is %s"%jid
    # now submit the models
    for m in model_list:
        # need to put the post processing job release command in the model somehow. Depends on the model
        # but we have a mark and a file. So will modify the file. The model should define this..
        # and insert the mark into the file. Would I think be easier to keep the line no and goto that.
        for line in fileinput.input(m.postProcessFile, inplace=1, backup='.bak2'):
            # if m.postProcessFile does not exist then  get an error which is what we want!
            # fix your model method!
            print line[0:-1] # just print the line out.
            if m.postProcessMark in line: # got the mark so add some text.
                print sshCmd,'qrls ',m.jid,'"' # this releases the post processing job.
        # dealt with modifying main file.
        modelSubmitName=m.submit()
        if verbose: print "Submitting ",modelSubmitName
        subprocess.check_output(sshCmd+modelSubmitName+'"',shell=True) # submit the script

    return True

# end of eddieSubmit
