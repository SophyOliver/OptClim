"""
Class to support mini mops model 18-May-2018 Sophy Oliver & Simon Tett
"""

import six
import collections
import pandas as pd
import os


from OptClimVn2 import ModelSimulation

class MiniMOPS(ModelSimulation.ModelSimulation):
	"""
	This is a class to run mini mops model as part of OptClimVn2
	"""
	# if want to model init file get from HadCM3 (def init) and edit inputs
	# chosen to inherit everything for now
	
	def __init__(self, dirPath, create=False, refDirPath=None, name=None, 
						obsNames=["mis1", "mis2", "mis3"], # for dfols: obsNames=["mis1", "mis2", "mis3"],
						ppOutputFile='misfit_output.txt',# options for creating new study
						update=False, # options for updating existing study
						verbose=False,  parameters={}, ppExePath=None): #if ever want to specify param names then change this to pass dictionary
		"""
		Create an instance of MiniMOPS class. Default behaviour is to read from dirPath and prohibit updates.
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
		
		
		if verbose:
			print("init params_to_pass", parameters)
		super(MiniMOPS, self).__init__(dirPath,
				obsNames=obsNames, create=create, refDirPath=refDirPath, name=name,
				ppOutputFile=ppOutputFile, parameters=parameters,  # options for creating new study
				update=update,  # options for updating existing study
				verbose=verbose)
		if create:
			self.setReadOnly(False)
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
			print(obsdf)
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
		print("MiniMOPS setParams def running ***********")
		if addParam == False:
			# need to find a way of resetting namelist files.
			# one option would be to copy the namelist files from the refdir. That would require working out all the files
			# that is moderately tricky and not yet needed. So raise anNotImplementedError if tried.
			raise exceptions.NotImplementedError("Not yet implemented addParam")

		super(MiniMOPS,self).setParams(params,addParam=addParam, write=write,verbose=verbose) # use super class setParams
		#To do: write text file and tests for this in test script
		# same as misfit output
		
		# Edit with Simon if necessary:
		# I got confused because wasn't sure where I was writing the parameters to so I hard coded it to parameters_input,txt?
		paramFile = os.path.join(self.dirPath,'parameters_input.txt')
		print(paramFile, params)
		keys=[k for k in params.keys() if (k.find('param',0) ==0)]
		keys.sort()
		with open(paramFile, "w") as text_file:
			for k in keys:
				text_file.write("%f\n"%params[k])
	
	
	def submit(self):
		"""
		Submit
		"""
		return os.path.join(self.dirPath,'runscript')