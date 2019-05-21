import collections
import os
import netCDF4
import shutil
import stat
import json
import pickle
import tempfile
import copy
import f90nml # available from http://f90nml.readthedocs.io/en/latest/
import exceptions
import numpy as np
import pandas as pd
import StudyConfig # needed for dealing with JSON files...

def errorRemoveReadonly(func, path, exc):
    """
    Function to run when error found in rmtree.
    :param func: function being called
    :param path: path to file being removed
    :param exc: failure status
    :return: None
    """
    import errno
    import stat
    excvalue = exc[1]
    if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
        # change the file to be readable,writable,executable: 0777
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
        print "changed permission on ",path
        # retry
        try:
            func(path)
        except WindowsError: # windows error is hard to do much about...
            pass

class modelEncoder(json.JSONEncoder):
    def default(self, obj):
        """
        Convert obj into something that can be converted to JSON
        :param obj -- object to be converted
        :return "primitive" objects that can be converted + type of object
        """
        stype = str(type(obj)) # string rep of data-type
        # try various approaches to generate json file
        fnsToTry=['to_dict','tolist']
        for fn in fnsToTry:
            f=getattr(obj,fn, None)
            if f is not None:
                return dict(stype=stype, data=f()) # run the method we've found
        # failed to run fns that made stuff can turn into JSON test
        if 'dtype' in dir(obj):
            return dict(data=str(obj), dtype=stype)
        else:
            return json.JSONEncoder(self, obj)  # Let the base class default method try to convert and raise the TypeError

_namedTupClass=collections.namedtuple('TupClassxx',('var','namelist','file')) # named tupple for namelist information
class ModelSimulation(object):

    """
    Class to define Model Simulations. This class is top class that
       provides basic functionality. In particular the configuration as updated
       is written to disk as a pickle file. Each set does this.
       Define new model simulations that do more by sub-classing this class.
       To get the observed values do model.getObs()
       To get parameters do model.getParams()
       To see what you have do model.get() and examine the orderedDict you have got.
       model.readObs() reads in data from data on disk. This is provided so that when a model
       has been run (and postprocessed) data can be read in and manipulated by functions outside this module.
       To support usecase of model failure readObs method will fall back to reading a json StudyConfig file and extracting
       data from it. See readObs for more details.

       TODO -- write more extensive documentation  esp on naemlist read/writes.

    """


    _simConfigPath= 'simulationConfig.cfg' # default configuration filename for the Model Simulation.


    def __init__(self, dirPath, obsNames=None,
                 create=False, refDirPath=None, name=None, ppExePath=None, ppOutputFile=None, parameters=None,  # options for creating new study
                 update=False, # options for updating existing study
                 verbose=False):
        """
        Create an instance of ModelSimulation class. Default behaviour is to read from dirPath and prohibit updates.
        :param dirPath -- path to directory where model simulation exists or is to be created. Shell variables and ~ will be expanded
        :param create (optional with default False). If True create new directory and populate it.
            If directory exists it will be deleted. 
            Afterwards the ModelSimulation will be readOnly.
            These options should be specified when creating a new study otherwise they are optional and ignored
            :param refDirPath -- reference directory. Copy all files from here into dirPath
            :param name -- name of the model simulation. If not provided will be taken from dirPath
            :param ppExePath --  path to post proessing executable
            :param ppOutputFile -- Filename for  output of post processing executable
            :param params -- dict of parameter names and values.
            :param obsNames -- list of observations to be readin. (see readObs())
        :param update -- allow updates to the simulation information.
        :param verbose -- provide  verbose output. (See individual methods). Default is False.
        :returns initialised object.
        """
        # stage 1 common initialisation
        self._convNameList = collections.OrderedDict() # information on the variable to Namelist mapping
        self._metaFn = collections.OrderedDict() # functions used for meta parameters.
        self.config=collections.OrderedDict()
        self.dirPath = os.path.expandvars(os.path.expanduser(dirPath))  # expand dirPath and store it
        # convert to an absolute path -- needed when jobs get submitted as don't know where we might be running from
        self.dirPath=os.path.abspath(self.dirPath)
        self._readOnly=True # default is read only.
        self._configFilePath = os.path.join(self.dirPath, self._simConfigPath)  # full path to configuration file
        self._failJsonFile = 'fail.json' # name of last ditch file to read for observations. Really here for case when model fails to run
        # but algorithm requires some data. Making it  json file means one can easily hack it in...
        # verify not got create and update both set
        if create and update: raise  Exception("Don't specify create and update")
        # create new model Simulation ?
        if create:
            self.createModelSimulation(parameters=parameters, ppExePath=ppExePath, obsNames=obsNames, name=name,
                                       ppOutputFile=ppOutputFile, refDirPath=refDirPath, verbose=verbose)
        else:
            self.readModelSimulation(update=update, verbose=verbose)
        # done reading in information

    def readConfig(self,verbose=False):
        """
        Read configuration file
        :param verbose (Default = False) if set then print out more information
        :return: configuration
        """

        with open(self._configFilePath, 'r') as fp:
            config = pickle.load(fp)

        return config

    def readModelSimulation(self, update=False, verbose=False):
        """
        Read in information on modelSimulation
        :param update (default = False) -- this model sim can be updated if true
        :param  verbose(default = False) be verbose
        :return:
        """

        # need to read config file
        config=self.readConfig(verbose=verbose)
        # for update I think failure to read is fine. We'd just return null...
        self.set(config, write=False,verbose=verbose)
        self.readObs(verbose=verbose)  # and read the observations
        self._readOnly = not update


    def set(self,keys_values,write=True, verbose=False): # suspect can do this with object methods...
        """
        sets value in configuration information and writes configuration out.
        :param keys_values:  dict/ordered dict of keys and values
        :param write (default = True). If true writes configuration information out to disk
        :param verbose (default = False). If true produces verbose output.

        :return:
        """
        # raise exception if readOnly.
        if self._readOnly and write: raise Exception("Config is read only")
        # set value in config
        for k, v in keys_values.iteritems():
            self.config[k] = v
        # write whole of config out to dir.
        if write:
            if os.path.isfile(self._configFilePath):
                os.chmod(self._configFilePath,stat.S_IWUSR) # make it write only so we can delete it
                os.remove(self._configFilePath) # remove file

            with open(self._configFilePath, 'w') as fp: # create it
                pickle.dump(self.config, fp)
            mode= os.stat(self._configFilePath).st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
            if verbose: print "Written current config to %s with mode %o " % (self._configFilePath, mode)
            os.chmod(self._configFilePath, mode)  # set its mode

    def get(self, keys=None, verbose=False):
        """
        get values from configuration.
        :param keys: keys to extract from configuration. If not defined returns the entire configuration.
        :param verbose (Optional -- default = False). If true print out some information.
        :return: orderedDict indexed by keys
        """

        if keys is None:
            return self.config # just return the entire config

        if type(keys) in (str, unicode):
            # keys is a single string or unicode
            return self.config[keys]

        # otherwise assume list and work through the elements of it.
        result = []
        for k in keys: result.append(self.config[k])
        # convert 1 element list to scalar
        if len(result) == 1: result=result[0]
        return result

    def createModelSimulation(self, parameters, ppExePath=None, obsNames=None, name=None,
                              ppOutputFile='observations.nc', refDirPath=None, verbose=False):
        # TODO make ppExePath & obsNames be optional. They should be or none be specified.
        """
        Create (in filesystem) a model simulation. After creation the simulation will be read only.
        :param parameters -- dict of parameter names and values OR pandaas series.
        :param ppExePath --  path to post processing executable -- Default None
        :param obsNames -- list of observations being used. -- Default None
        :param (optional)name -- name of the model simulation. If not provided will be taken from dirPath
        :param (optional) ppOutputFile -- name of output file where output from postprcessing is (default = observations.nc) 
        :param (optional) refDirPath -- reference directory. Copy all files from here into dirPath
        :param (optional) verbose -- if true be verbose. Default is False
        """
        # general setup
        self._readOnly=False # can write if wanted.

        if refDirPath is not None:
            refDirPath=os.path.expandvars(os.path.expanduser(refDirPath))
        #  fill out configuration information.
        config=collections.OrderedDict()
        if name is None:
            config['name'] = os.path.basename(self.dirPath)
        else:
            config['name'] = name

        obs=collections.OrderedDict()
        try:
            for k in obsNames: obs[k]=None
        except TypeError:
            pass

        config['ppExePath'] = ppExePath
        config['ppOutputFile'] = ppOutputFile
        if refDirPath is not None: config['refDirPath'] = refDirPath

        config['observations'] = obs
        config['parameters'] = parameters
        # check everything is defined. Commented out as don't think needed..
        #for k, v in config.iteritems():
        #    if v is None: raise Exception("Specify %s" % k)

        if verbose:   print "Config is ", config



        if os.path.exists(self.dirPath): # delete the directory (if it currently exists)
            shutil.rmtree(self.dirPath,onerror=errorRemoveReadonly)


        if refDirPath is not None:  # copy all files and directories from refDirPath to dirPath and make then user writable
            mode = os.umask(0)  # get current vale of umask
            os.umask(mode)  # and reset it
            mode = mode | stat.S_IWUSR  # mode we want is umask + write

            shutil.copytree(refDirPath, self.dirPath)  # copy everything
            # now change the permissions so user can write
            if verbose: print "Copied files from %s to %s " % (refDirPath, self.dirPath)
            for root, dirs, files in os.walk(refDirPath, topdown=True):  # iterate
                for name in dirs:  # iterate over  the directories
                    dname = os.path.join(root, name)
                    mode=os.stat(dname).st_mode |  stat.S_IWUSR
                    os.chmod(dname, mode)
            for name in files:  # iterate over files in directory
                fname = os.path.join(root, name)
                mode=os.stat(fname).st_mode |  stat.S_IWUSR # mode has user write permission.
#                print "Setting mode on %s  to %d"%(fname,mode)
                os.chmod(fname, mode)
        else:
            # Make target directory.
            os.makedirs(self.dirPath)  # create directory
        self.set(config) # set (and write) configuration
        # and no longer able to write to it.
        self._readOnly=True
    # end of createModelSimulation

    def name(self):
        """
        Return the name of the Model Configuration
        """

        return self.config['name']

    def ppOutputFile(self):
        """
        Return the name of the file to be produced by the post Processing 
        :return: 
        """

        return self.get(['ppOutputFile'])

    def ppExePath(self):
        """
        Return the name of the post processing executable file. 
        :return: 
        """

        return self.get(['ppExePath'])


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


        if fileType == '.nc':
            # netcdf file so read in data allowing python stuff to raise exception if it doesn't exist or is corrupt.
            # current implementation is missing the constraint..
            if verbose: # provide some helpful info.
                print "Reading netcdf data from %s "%obsFile
                print "For: ",obs.keys()

            with netCDF4.Dataset(obsFile, "r")  as ofile: # open it up for reading
                if verbose: print "nc files got ",ofile.variables.keys()
                for var in varsWant: # loop over variables
                    if var in ofile.variables: # got it?
                        obs[var]=ofile.variables[var][0] # set it NB this  will break if obs is not a scalar.
                    else:
                        obs[var]=None # set to None if we don't have it.
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

    def getObs(self,verbose=False):
        """
        Extract the observations
        :param verbose (optional -- default False) If True print out information
        :return: an ordered  dict of the observations
        """
        obs = self.get("observations", verbose=verbose)
        return obs

    def setParams(self,params, addParam=False, write=True, verbose=False):
        """
        Set the parameter values and write them to the configuration file
        :param params -- dictionary (or ordered dict) of the parameter values
        :param addParam (default False) -- if True add to existing parameters
        :param write (default True) -- if True update configuration file.
        :param verbose (default False) -- if True provide more verbose output
        :return:
        """
        if params is None:
            raise Exception
        if addParam:
            if verbose: print "Updating existing parameters"
            oldParams=self.getParams()
            for k,v in params.iteritems():
                oldParams[k]=v # update
            self.set({'parameters':oldParams},verbose=verbose, write=write)
        else:
            self.set({'parameters':params},verbose=verbose, write=write) # set the parameters


    def getParams(self, verbose=False, params=None):
        """
        Extract the parameter values
        :param verbose (default = False) if True print out information
        :param params (default is None) If not None then only return those parameters in the same order.
        :return: the parameter values
        """

        p=self.get('parameters',verbose=verbose)
        if params is not None: # enforce particular order.
            t=collections.OrderedDict()
            for pp in params: t[pp]=p[pp]
            p=t
        
        return p

    def registerMetaFn(self,varName, function, verbose=False):
        """
        Register a function to process meta parameters
        The function should take three  argument -- a value for forward call and a dict with names and values for inverse call
        , a named argument inverse with default False, and a named argument naemlist. value should be optional and default valaue is something sensible.
        It should return a dict of values keyed by namelist tuple.
        TODO: Add an example of how this works -- see test_modelSimulation.py
        :param varName: name of meta parameter
        :param function: Function to register
        :param verbose: be verbose. Default is False
        :return: nothing
        """
        if verbose:
            try:
                print "Registering function %s for %s"%(function.func_name, varName)
            except AttributeError:
                pass
        # verify fn works as expected.
        res=function() # run it with default value.
        keys=res.keys()

        # verify namelist keyword works
        nlKeys=function(namelist=True)
        if not ((set(nlKeys) == set(keys)) and (len(nlKeys) == len(keys))):
            print("nlKeys are",set(nlKeys))
            print("keys are ",set(keys))
            raise exceptions.Exception("namelist keys and expected keys differ")
        # verify nlKeys and keys are the same.

        a = function(res, inverse=True)  # make sure inverse works.
        res2=function(a) # having got the default arg back we run again,
        for k,v in res2.iteritems():
            if type(v) is np.ndarray:
                l = np.all(v != res[k])
            else:
                l=(v != res[k])
            if l:
                raise exceptions.Exception("Fn not invertable")

        self._metaFn[varName]=function


    def genVarToNameList(self, param, nameListVar,  nameListName, nameListFile,  verbose=False):
        """
        Generate a conversion list for use in converting model parameter names to namelist values
        :param param: the name of the parameter (as a string)
        :param nameListVar  :  variable name in the namelist
        :param nameListName: name of namelist
        :param nameListFile: Paths *relative* to the configuration dir for file containing namelist
        :param verbose (optional -- default = False) be verbose if True
        :return: a named tuple containing variable, namelist and file.
        """
        nt=_namedTupClass(var=nameListVar,namelist=nameListName,file=nameListFile)
        self._convNameList[param]=[nt]

        if verbose: print "var %s ->"%(param), self._convNameList[param]


    def applyMetaFns(self, verbose=False, fail=False, **params):
        """
        Apply transformation functions to meta parameters
          :param verbose: optional default is False. If True be verbose.
          :param fail: optional default is False. If true fail if parameter not found
          keywords are parameters and value.
        :return: returns dict of keys (which should be named tuple defining namelist var, namelist and file)  and their values,
            and list of meta parameters that were used.



        """
        result=dict()
        metaParamsUsed=[]
        for k,v in params.iteritems():
            if k in self._metaFn or fail:
                # got fn (or not fail) so run it
                fnResult=self._metaFn[k](v)
                metaParamsUsed.append(k)
                for fk, fv in fnResult.iteritems():
                    result[fk]=fv

        # sort metaParamsUsed so order is deterministic...
        return  result,sorted(metaParamsUsed)

    # TODO (if ever need it) add a NameList method that returns the namelist info possibly nicely printed.

    def writeNameList(self, verbose=False, fail=False, **params ):
        """
        Modify existing namelist files using information generated via genConversion
        Existing files will be copied to .bak
        :param verbose (optional -- default is False). If True provide more information on what is going on.
        :param fail (optional default is False). If True fail if a parameter not found.
        :keyword arguments are parameters and values.
        :return:  ordered dict of parameters and values used.
        """
        if self._readOnly:
            raise exceptions.IOError("Model is read only")

        params_used=collections.OrderedDict() #
        files=collections.OrderedDict() # list of files to be modified.
        for param,value in params.iteritems(): # extract data from conversion indexed by file --
            # could this code be moved into genVarToNameList as really a different view of the same data.
            # NO as we would need to do this only once we've finished generating namelist translate tables.
            # potential optimisation might be to cache this and trigger error in writeNameList if called after genNameList
            # search functions first
            if param in self._metaFn: # got a meta function.
                if verbose: "Running function %s"%self._metaFn[param].func_name
                metaFnValues = self._metaFn[param](value) # call the meta param function which returns a dict
                params_used[param] = metaFnValues  # and update return var
                for conv,v in metaFnValues.iteritems(): # iterate over result of fn.
                    if conv.file not in files:
                        files[conv.file] = []  # if not come across the file set it to empty list
                    files[conv.file].append((v, conv))  # append the  value  & conversion info.
            elif param in self._convNameList: # got it in convNameList ?
                for conv in self._convNameList[param]:
                    if conv.file not in files:
                        files[conv.file]=[]    # if not come across the file set it to empty list
                    files[conv.file].append((value,conv)) # append the value  & conversion
                    params_used[param]=value # and update return var
            elif fail:
                raise exceptions.KeyError("Failed to find %s in metaFn or convNameList "%param)
            else:
                pass

        # now have conversion tuples ordered by file so let's process the files
        for file in files.keys(): # iterate over files
            # need to create backup? Only do if no back up exists. This allows generateNameList to be run multiple times
            # doing updates. First time it runs we assume we have a directory ready to be modified.
            filePath=os.path.join(self.dirPath,file) # full path to namelist file
            # check file exists if not raise exception
            if not os.path.isfile(filePath):
                #raise exceptions.IOError("file %s does not exist"%(filePath))
                continue # skip this file.
            backup_file=filePath+"_nl.bak" # and full path to backup fie.
            if not os.path.isfile(backup_file):
                shutil.copyfile(filePath,backup_file)
            # now create the namelist file. Need a temp file
            with tempfile.NamedTemporaryFile(dir=self.dirPath,delete=False) as tmpNL:
                # Now construct the patch for the  namelist file for all conversion tuples.
                nlPatch=collections.OrderedDict() # path to exisiting namelist file
                for (value, conv) in files[file]:
                    if conv.namelist not in nlPatch:
                        nlPatch[conv.namelist]=collections.OrderedDict() # dom't have ordered dict so make it
                    if type(value) is np.ndarray: # convert numpy array to list for writing.
                        value=value.tolist()
                    elif isinstance(value,unicode):
                        value=str(value) # f90nml can't cope with unicode so convert it to string.
                    nlPatch[conv.namelist][conv.var]=copy.copy(value) # copy the variable to be stored rather than the name.
                    if verbose: print "Setting %s,%s to %s in %s"%(conv.namelist,conv.var,value,filePath)
                try:
                    p=f90nml.patch(filePath,nlPatch,tmpNL.name) # patch the namelist file
                    tmpNL.close() # close the temp file once done.
                except StopIteration:
                    print "Problem in f90nml for %s writing to %s"%(filePath,tmpNL.name), nlPatch
                    raise # raise exception.

                if verbose: print "Patched %s to %s"%(filePath,tmpNL.name)
                shutil.move(tmpNL.name,filePath) # and copy the patched file back in place.

        return  params_used

    def readNameList(self, params, fail=False, verbose=False):
        """
        Read parameter value from registered namelist
        :param fail: If True fail if param not found
        :param verbose (Optional -- default False). Provide verbose information.
        :param params -- a list of parameters.
        :example self.readNameList(['RHCRIT', 'VF1'])
        :return:An OrderedDict with the values indexed by the param names
        """

        result=collections.OrderedDict()
        for param in params:
            # have it as meta function?
            if param in self._metaFn: # have a meta funcion -- takes priority.
                result[param]=self.readMetaNameList(param, verbose=verbose)
            elif param in self._convNameList:  # in the conversion index
                nlValue=self.readNameListVar(self._convNameList[param],verbose=verbose)
                if len(nlValue) != 1: raise exceptions.ValueError("Should only have one key")
                for k in nlValue.keys(): # in this case shoudl only have one key and we don't want it as a list
                    result[param]= nlValue[k]
                # just want the value and as single param SHOULD return
            elif fail: # not found it and want to fail
                raise  exceptions.KeyError("Param %s not found"%param)
            else:
                pass
        return result # return the result.

    def readNameListVar(self,nameListVars, verbose=False):
        """
        Read single parameter specified via named tuple defining namelist variable

        :param verbose: default False. If True be verbose
        :param NameListVars: list of namelist variables to be retrieved
        :return: an ordered dict indexed by namelist info (if found) with values retrieved.
        """

        result=collections.OrderedDict()
        namelists={}
        # iterate over all parameters reading in namelists.
        for var in nameListVars:
            if var.file not in namelists: # not already read this namelist.
                namelists[var.file]=f90nml.read(nml_path=os.path.join(self.dirPath, var.file))
        # now having read all needed namelists extract the values we want
        for var in nameListVars:
            nlvalue = namelists[var.file][var.namelist][var.var]
            result[var]=nlvalue
        return result


    def readMetaNameList(self,param, verbose=False):
        """
        Retrieve value of meta parameter  by reading namelists and running inverse function.
        :param param:  name of meta-parameter
        :param verbose: be verbose
        :return:  value of meta-parameter
        """
        # work out what fn is and run it with default value to work out what namelist values we need to retrieve.
        fn=self._metaFn[param]  # should generate an error if not found
        nlInfo=fn(namelist=True) # get teh namelist info by askugn the function for it!
        # retrieve appropriate values from namelist
        if verbose: print "Retrieving ",vars
        var=self.readNameListVar(nlInfo, verbose=verbose) # read from namelist.
        if verbose: print "Retrieved",var
        result=fn(var, inverse=True) # run function in reverse
        return  result

    def setReadOnly(self,readOnly=None):
        """
        Modify read only flag
        :param readOnly: value to set (default is None -- which doesn't modify status. True means any attempts to write to dir will trigger error
        :return: current status
        """

        if readOnly is not None:
            self._readOnly=readOnly

        return self._readOnly

    def allParamNames(self):
        """
        
        :return: list of all parameter names -- suitable  for passing through.
        """

        names=self._convNameList.keys()
        names.extend(self._metaFn.keys())
        return names


class EddieModel(ModelSimulation):
    """
    A simple model suitable for running on Eddie -- this allows testing of the whole approach.
    """

    def __init__(self, dirPath, obsNames=None,
                 create=False, refDirPath=None, name=None, ppExePath=None, ppOutputFile=None, studyConfig=None, # options for creating new study
                  update=False, # options for updating existing study
                  verbose=False, parameters={}):
        """
        Create an instance of HadCM3 class. Default behaviour is to read from dirPath and prohibit updates.
        :param dirPath -- path to directory where model simulation exists or is to be created
        :param create (optional with default False). If True create new directory and populate it.
            Afterwards the ModelSimulation will be readOnly.
            These options should be specified when creating a new study otherwise they are optional and ignored
            :param refDirPath -- reference directory. Copy all files from here into dirPath
            :param name -- name of the model simulation. If not provided will be taken from dirPath
            :param ppExePath --  path to post proessing executable
            :param ppOutputFile -- output file for  post processing executable
            :param obsNames -- list of observations to be readin. (see readObs())
            :param studyConfig -- written into directory.
            :param parameters -- a dict on pandas series specifying hte parameter values
        :param update -- allow updates to the simulation information.
        :param verbose -- provide  verbose output. (See individual methods). Default is False.
        :returns initialised object.
        """


        # no parameters should be provided unless create or update provided
        if len(parameters) >0 and  not (create  or update):
            raise exceptions.ValueError("Provided parameters but not specified create or update")

        # call superclass init
        super(EddieModel, self).__init__(dirPath,
                obsNames=obsNames, create=create, refDirPath=refDirPath, name=name, ppExePath=ppExePath,
                ppOutputFile=ppOutputFile, parameters=parameters,  # options for creating new study
                update=update,  # options for updating existing study
                verbose=verbose)
        if studyConfig is not None:
            studyConfig.save(filename=os.path.join(self.dirPath,"config.json")) # write study configuration for fakemodel.



    def submit(self):
        """
        Provides full path to submit script.
        :return: path to submit (this get's run)
        """

        return os.path.join(self.dirPath,'submit.sh')



