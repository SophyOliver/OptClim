"""
Class to support mops model 18-May-2018 Sophy Oliver & Simon Tett
"""

import six
import collections
import pandas as pd
import os


from OptClimSO import ModelSimulation

class MOPS(ModelSimulation.ModelSimulation):
	"""
	This is a class to run mops model as part of OptClimVn2
	"""
	# if want to model init file get from HadCM3 (def init) and edit inputs
	# chosen to inherit everything for now
	
	def __init__(self, dirPath, create=False, refDirPath=None, name=None, 
						obsNames=["mis1_PO4", "mis2_PO4", "mis3_PO4", "mis4_PO4", "mis5_PO4", "mis6_PO4", "mis7_PO4", "mis8_PO4", "mis9_PO4", "mis10_PO4", "mis11_PO4", "mis12_PO4", "mis13_PO4", "mis14_PO4",
						"mis1_NO3", "mis2_NO3", "mis3_NO3", "mis4_NO3", "mis5_NO3", "mis6_NO3", "mis7_NO3", "mis8_NO3", "mis9_NO3", "mis10_NO3", "mis11_NO3", "mis12_NO3", "mis13_NO3", "mis14_NO3",
						"mis1_O2",  "mis2_O2",  "mis3_O2",  "mis4_O2",  "mis5_O2",  "mis6_O2",  "mis7_O2",  "mis8_O2",  "mis9_O2",  "mis10_O2",  "mis11_O2",  "mis12_O2",  "mis13_O2",  "mis14_O2"],
						ppOutputFile='misfit_output.txt',# options for creating new study
						update=False, # options for updating existing study
						verbose=False,  parameters=[], ppExePath=None, runscriptName='runscript'): #if ever want to specify param names then change this to pass dictionary
		"""
		Create an instance of MOPS class. Default behaviour is to read from dirPath and prohibit updates.
		:param dirPath -- path to directory where model simulation exists or is to be created
		:param obsNames -- names for the ocean tracers we're comparing to
		:param create (optional with default False). If True create new directory and populate it.
			Afterwards the ModelSimulation will be readOnly.
			These options should be specified when creating a new study otherwise they are optional and ignored
			:param refDirPath -- reference directory. Copy all files from here into dirPath
			:param name -- name of the model simulation. If not provided will be taken from dirPath
			:param ppOutputFile -- File name of output of model executable. Default is misfit_output.txt
		:param update -- allow updates to the simulation information.
		:param verbose -- provide  verbose output. (See individual methods). Default is False.
		:param parameters -- list of parameter values to optimise in the order required by model
		:returns initialised object.
		"""
		
		# no parameters should be provided unless create or update provided
		if len(parameters) >0 and  not (create  or update):
			raise ValueError("Provided parameters but not specified create or update")
		
		# code to error check the parameter values (TO DO)
		#
		#
		
		
		#print("init params_to_pass", parameters)
		super(MOPS, self).__init__(dirPath,
				obsNames=obsNames, create=create, refDirPath=refDirPath, name=name,
				ppOutputFile=ppOutputFile, parameters=parameters,  # options for creating new study
				update=update,  # options for updating existing study
				verbose=verbose, runscriptName=runscriptName)
		if create:
			self.setReadOnly(False)
			print("init params_to_pass", parameters)
			self.setParams(parameters, verbose=verbose)
			self.setReadOnly(True)
			
	def readObs(self,verbose=False):
		"""
		Read the post processed data. This default implementation reads netcdf data and stores it in the configuration
		If the desired observation  does not exist then the method attempts to load data from _failJsonFile (fail.jsno)
		and use simulatedObservations entry to get data from.
		:param verbose (default False). If true print out helpful information
		:return: a ordered dict of observations wanted. Values not found when requested will be set to None.
		"""

		obsFile=os.path.join(self.dirPath, self.ppOutputFile())
		
		# Code to read in cmaes-specific file if it exists. Directs us to the best obs file of the cmaes population
		cmaesFile = os.path.join(self.dirPath, 'best_cmaes.txt')
		if os.path.isfile(cmaesFile):
			with open(cmaesFile) as f:
				nI = f.readline().strip()
				besti = f.readline().strip()
			ppOut = self.ppOutputFile()
			ppOut = ppOut[0:-4]
			obsFile=os.path.join(self.dirPath, ppOut + '_' + str(nI) + '_' + str(besti) + '.txt')
		
		obsJsonFile=os.path.join(self.dirPath, self._failJsonFile)
		obs=collections.OrderedDict()
		varsWant = self.get(['observations']).keys()  ## extract keys from current obs
		# see if obs file exists or json file exists. If not return dict with Nones.
		if not (os.path.exists(obsFile) or os.path.exists(obsJsonFile)):
			if verbose: print "File %s not found. Returning empty dict"%(obsFile)
			for v in varsWant:
				obs[v]=None
			return obs
		# now maybe set obsFile to obsJsonFile
		if not (os.path.exists(obsFile)):
			print("WARNING: %s not found. Trying %s instead. " % (obsFile, obsJsonFile))
			obsFile=obsJsonFile # make obs the json file if the obs file does not exist.

		fileType=os.path.splitext(obsFile)[1]  # type of file


		if fileType == '.txt':
			# text file so read in data allowing python stuff to raise exception if it doesn't exist or is corrupt.
			# current implementation is missing the constraint..
			if verbose: # provide some helpful info.
				print "Reading text data from %s "%obsFile
				print "For: ",obs.keys()
			
			obsdf = pd.read_csv(obsFile, header=None, index_col=False)
			#print(obsdf)
			for n,k in enumerate(varsWant):
				obs[k] = obsdf.iloc[n].values[0]
			
			
		elif fileType == '.json':
			# this should be a configuration readable file..
			config=StudyConfig.readConfig(obsFile,ordered=True) # read the config.
			obsData=config.simulatedObservations
			if verbose: print "json file got ", obsData.keys()
			for var in varsWant:  # loop over variables
				obs[var]=obsData.get(var,None) # get variable returning None if not defined.

		else: # don't know what to do. So raise an error
			raise NotImplementedError("Do not recognize %s"%fileType)
		self.set({'observations':obs},write=False)
	
	def setParams(self,params, addParam=True, write=True, verbose=False, fail=True):
		"""
		Set the parameter values and write them to the configuration file
		and modify the parameters in the current directory. Calls the superclass to do standard stuff first then
		uses existing code to modify parameters in HadCM3 namelists.
		:param params -- dictionary (or ordered dict or pd.Series) of the parameter values
		:param addParam (default False) -- if True add to existing parameters
		:param write (default True) -- if True update configuration file.
		:param verbose (default False) -- if True provide more verbose output
		:param fail (default True) -- if True fail if parameter not defined.
		:return:
		"""
		print("MOPS setParams def running ***********")
		if addParam == False:
			# need to find a way of resetting namelist files.
			# one option would be to copy the namelist files from the refdir. That would require working out all the files
			# that is moderately tricky and not yet needed. So raise anNotImplementedError if tried.
			raise exceptions.NotImplementedError("Not yet implemented addParam")
			
		super(MOPS,self).setParams(params,addParam=addParam, write=write,verbose=verbose) # use super class setParams
		#To do: write text file and tests for this in test script
		# same as misfit output
		
		paramFile = os.path.join(self.dirPath,'parameters_input.txt') # file only for BOBYQA/DFOLS
		paramNamesFile = os.path.join(self.dirPath,'parameters_names.txt')
		num_bgc_params_file = os.path.join(self.dirPath,'num_bgc_params.txt')
		#print(paramFile, params)
		keys=params.keys();
		
		# write parameter names and number of params files
		with open(paramNamesFile, "w") as text_file_names, open(num_bgc_params_file, "w") as num_params:
			np = 0 # Count parameters
			vals2write = 0 # Count parameter values specifically given
			for k in keys:
				if k != "RUNID":
					text_file_names.write("%s\n"%k)
					if params[k] is not None:
						vals2write = vals2write + 1
					np = np + 1
			num_params.write("%d\n"%np)
			
		# Write out parameter values to parameters_input.txt as well (if we have any! This should not happen if we use cmaes)
		if vals2write > 0:
			with open(paramFile, "w") as text_file_vals:
				for k in keys:
					if k != "RUNID":
						if params[k] != "RUN_CMAES":
							text_file_vals.write("%.15f\n"%params[k])
	
	def submit(self):
		"""
		Submit
		"""
		return os.path.join(self.dirPath,self.runscriptName())