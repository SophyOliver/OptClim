"""
Provides classes and methods suitable for manipulating study configurations.  Includes two useful classes:
    fileDict which designed to provide some kind of permanent store across  invocations of the framework. 

    studyConfig which inherits from fileDict but "knows" about configuration for study and provides methods to
       manipulate it. Idea being to insulate rest of code from details of configuration file and allow complex processing 
       if necessary.

Changes at time of writing (14/1/16) are partial. Much of the code uses classes in this module
but there is still code in optClimFwk which defines its own study/liveset which could be removed.

"""
import json
import collections
import numpy as np
import pandas as pd
import os
import xarray
import six
import numpy.linalg as linalg

# functions available to everything.

def readConfig(filename, ordered=False):
    """
    Read a configuration and return object of the appropriate version.
    :param filename: name of file to read.
    :param ordered:  read in as a ordered dict rather than a dict. 
    :return: Configuration of appropriate type
    """
    path=os.path.expanduser(os.path.expandvars(filename))

    if os.path.isfile(path) is False:
        raise IOError("File %s not found"%filename)
    config=dictFile(filename=path, ordered=ordered) # read configuration using rather dumb object.
    vn = config.getv('version', default=None)
    if isinstance(vn,six.string_types):
        vn = float(vn)  # convert to float as version stored as string

    # use appropriate generator fn
    if vn is None or vn < 2: # original configuration definition
        config=OptClimConfig(config)
    elif vn < 3: # version 2 config
        config=OptClimConfigVn2(config)
    else:
        raise Exception("Version must be < 3")

    return config

def getDefault(dct,key,default):
    """
    :param dct: dictionary to read value from
    :param key: key to use
    :param default: default value to use if not set or None
    :return: value from dct if provided, default if not provided or None
    """
    value=dct.get(key, default)
    if value is None: # value is none then use default
        value = default
    return value # and return the value

# class to deal with numpy in JSON.

class NumpyEncoder(json.JSONEncoder):
    # this code a fairly straight import from Mike's numpyJSON.
    # which provides numpy aware JSON encoder and decoder.
    def default(self, obj):
        """
        if input object is a ndarray it will be converted into a OrderedDict
        holding dtype, shape and the data
        """
        if isinstance(obj, np.ndarray):
            data_list = obj.tolist()
            return collections.OrderedDict(__ndarray__= data_list, dtype=str(obj.dtype), shape=obj.shape)
        elif 'dtype' in dir(obj):
             return collections.OrderedDict(__npdatum__=str(obj), dtype=str(obj.dtype))
        # Let the base class default method raise the TypeError
        return json.JSONEncoder(self, obj)


def decode(dct):
    """
    Decodes a previously encoded numpy ndarray
    with proper shape and dtype
    :param dct: (dict) json encoded ndarray
    :return: (ndarray) if input was an encoded ndarray
    """
    if isinstance(dct, dict):
        if  '__ndarray__' in dct:
            data = np.array(dct['__ndarray__']).reshape(dct['shape']) ## extract data
            return data  # already been converted etc so just return it
    elif '__npdatum__' in dct:
        data = dct['__npdatum__']
        return data

    return collections.OrderedDict(dct)


class dictFile(dict):
    """
    extends dict to prove save and load methods.
    """
    def __init__(self, filename = None, ordered=False):
        """
        Initialise dictFile object from file
        :param (optional) filename -- name of file to load
        :param (optional) ordered -- if set True then load data as orderedDict. This is incompatable with
           decoding numpy data. TODO: Understand how to decode numpy data AND make an orderedDict
        :param *arg: arguments
        :param kw: keyword arguments
        :return: just call dict init and returns that.
        """
        if filename is not None:
            fname=os.path.expandvars(filename)
            fname=os.path.expanduser(fname)
            try:
                with open(fname, 'r') as fp:
                    if ordered:
                        dct=json.load(fp, object_pairs_hook=collections.OrderedDict)
                    else:
                        dct=json.load(fp, object_hook=decode)##,object_pairs_hook=collections.OrderedDict)
            except IOError: # I/O problem
                dct=collections.OrderedDict()

        else:
            dct=collections.OrderedDict() # make it empty ordered dict.

        self.Config=dct
        self._filename=filename

    def __getattr__(self, name):
        """
        Default way of getting values back. 
        :param name: attribute name to be retrieved
        :return: value 
        """
        # TODO -- consider making __getattr__ work if name not found -- when it returns None.
        # Disadvantage of returning None is that if key mistyped will get None so makes life more Error prone
        # One approach might be to have a list of names and then if key in there return None otherwise fail?
        try:
            return self.Config[name] # see if we have it
        except KeyError: # no we do not so raise AttributeError. Alt might be to return None...
            raise AttributeError("Failed to find %s in self.Config"%(name)) # spec for __getattr__ says raise AttributeError

    def print_keys(self):
        """
        Print out top-level keys
        :return: None
        """

        for k,v in self.Config.iteritems():
            print "%s: %s"%(k,v)

    def save(self,filename = None, verbose=False):
        """
        saves dict to specified filename.
        :param filename to save file to. Optional and if not provided will use
           private variable filename in object. If provided filename that self uses
           subsequently will be this filename
        :param verbose (optional, default False). If true be verbose

        :return: status of save.
        """
        if filename is not None:
            self._filename=filename # overwrite filename
        if filename is None and hasattr(self,"_filename") :
            filename=self._filename
        if filename is None: raise ValueError("Provide a filename")
        # simply write object to JSON file. Note no error trapping
        if os.path.isfile(filename): os.remove(filename)
        with open(filename, 'w') as fp:
            json.dump(self.Config,fp, cls=NumpyEncoder, indent=4)
        if verbose:
            print("Wrote to %s"%filename)
        #return None

    def getv(self,key,default=None):
        """
        Return values from Config component
        :return: value or None if key not defined
        """
        return self.Config.get(key,default)

    def setv(self,key,var):
        """

        :param key: Key to be set
        :param var: value to be set
        :return: modifies self
        """
        self.Config[key]=var

    def increment(self,key,increment):
        """

        :param key: key for which value we are incrementing
        :param increment: how much we increment key by. Be careful what you increment..
        :return: incremented value
        """
        default=0 # default value
        if isinstance(increment,(list,tuple)): # incrementing using a list or tuple so default becomes a list
            default=[] # empty list
        value=self.getv(key,default=default) # get the value with default
        value+=increment # increment the value
        self.setv(key,value) # set it
        return value # return it

    def fileName(self):
        """
        :return: path to filename
        """
        return self._filename



class OptClimConfig(dictFile):
    """
    Class to provide methods for working with OptclimConfig file.
    Inherits from dictFile which provides methods to load and save configurations
    All methods, where sensible, return a pandas dataframe or series.

    When extending this care needs to be taken that data in DataFrames/Series is float rather than objects.
    Objects rather than floats will be created if mixed strings/floats are made into dataframes/series.
    Most of the code avoids this by extracting only the required parameters and observations.

    """

    def __init__(self,config):
        """
        Return OptClimConfig object -- likely called using readConfig
        :param config: -- a dictFile configuration. Information will be copied from here. 
        """

        self.__dict__.update(config.__dict__) # just copy the values across!


    def version(self):
        """

        :return: the current version of the configuration file.
        """
        return self.getv('version',None)



    def name(self,name=None):
        """
        :param name -- if not None then set name in self to name
        :return: a short name for the configuration file
        """
        if name is not  None: self.setv("Name",name)
        name = self.getv("Name")
        if name is None: name='Unknown'
        return name

    def paramNames(self):
        """
        :return: a list of parameter names from the configuration files
        """
        return self.Config['study']['ParamList'][:] # return a copy of the list.

    def obsNames(self,add_constraint=None):
        """
        :param (optional default = False) add_constraint -- if True add the constraint name to the list of obs
        :return: a list  of observation names from the configuration files
        """
        obs=self.Config['study']['ObsList'][:] # return a copy of the array.
        if add_constraint:
            obs.append(self.constraintName())

        return  obs

    def paramRanges(self, paramNames= None):
        """
        :param paramNames -- a list of the parameters to extract ranges for.
        If not supplied the paramNames method will be used.
        :return: a pandas array with rows names minParam, maxParam, rangeParam
        """
        if paramNames is None: paramNames=self.paramNames()
        param = pd.DataFrame(self.Config['minmax'],
            index=['minParam', 'maxParam'])
        param = param.loc[:,paramNames] # just keep the parameters we want getting rid of comments etc in the JSON file
        param = param.astype(float)
        param.loc['rangeParam',:] = param.loc['maxParam',:]-param.loc['minParam',:] # compute range
        return param

    def standardParam(self, paramNames=None, normalise=False):
        """
        Extract standard parameter values for study
        :param paramNames: Optional names of parameters to use.
        :return: pandas series
        """
        if paramNames is None:   paramNames = self.paramNames()
        svalues = pd.Series([self.Config['Parameters']['defaultParams'].get(k,np.nan) for k in paramNames],
                            index=paramNames)
        if normalise: # normalise
            raise NotImplementedError

        return svalues.rename(self.name())

    def standardObs(self, obsNames= None, scale=False):
        """
        Extract standard obs values for study
        :param obsNames: Optional names of observations to use.
        :param scale (default False). If true scale values
        :return: pandas series
        """
        if obsNames is None:      obsNames=self.obsNames()
        svalues = pd.Series([self.Config['standardModel']['SimulatedValues'].get(k,np.nan) for k in obsNames],
                            index=obsNames, )
        if scale: svalues=svalues*self.scales(obsNames=obsNames) # maybe scale it.

        return svalues.rename(self.name())

    def standardConstraint(self,constraintName= None):
        """
        Extract constraint value from standard model values.
        :param constraintName (default None): value to extract
        :return: constrain as a pandas Series
        """
        if constraintName is None:
            constraintName=[self.constraintName()]

        series=pd.Series([self.Config['standardModel']['SimulatedValues'].get(k,np.nan) for k in constraintName],
                            index=constraintName)
        return series.rename(self.name())

    def beginParam(self,paramNames = None,scale=False):
        """
        get the begin parameter values for the study. These are specified in the JSON file in begin block
        Any values not specified use the standard values
        :param paramNames: Optional names of parameters to use.
        :param scale (default False). If True scale parameters by their range so 0 is minimum and 1 is maximum
        :return: pandas series of begin parameter values.
        """
        if paramNames is None:  paramNames=self.paramNames()
        begin = {} # empty dict
        standard=self.standardParam(paramNames=paramNames)
        scaleRange=self.Config['begin'].get("paramScale") # want to scale ranges?
        range=self.paramRanges(paramNames=paramNames) # get param range

        for p in paramNames: # list below is probably rather slow and could be sped up!
            begin[p] = self.Config['begin']['paramValues'].get(p)
            if begin[p] is not None:
                if scaleRange: # values are specified as 0-1
                    begin[p]=begin[p]*range.loc['rangeParam',p]+range.loc['minParam',p]
                if scale: # want to return them in range 0-1
                    begin[p]=(begin[p]-range.loc['minParam',p])/range.loc['rangeParam',p]

            if begin[p] is None: begin[p]=standard[p]

        begin=pd.Series(begin, dtype=float)[paramNames] # order in the same way for everything.

        # verify values are within range
        if scale: L=begin.gt(1.0) | begin.lt(0.0)
        else: L = range.loc['maxParam',:].lt(begin) | begin.lt(range.loc['minParam',:])

        if np.any(L):
            print "L  \n",L
            print "begin: \n",begin
            print "range: \n",range
            print "Parameters out of range",begin[L].index
            raise ValueError("Parameters out of range: ")

        return begin.astype(float).rename(self.name())

    def firstOptimise(self):

        firstOptimise = self.Config['begin'].get('firstOptimiseStep','GN')
        if firstOptimise is None: firstOptimise='GN'
        return firstOptimise

    def studyVars(self,studyDir):
        """
        :param studyDir the study directory
        :return: a studyVars object for the current study
        """
        file=self.studyFileStore()
        file=os.path.join(studyDir,file)
        return dictFile(filename=file)


    def studyFileStore(self):
        """
        Get the name of the file used to store and update information for the whole study.

        :return: filename relative to studyDir
        """
        fileStore = self.Config['begin'].get('studyFileStore','study_vars.json')
        if fileStore is None: fileStore='study_vars.json'
        return fileStore

    def cacheFile(self):
        """
        Get the pathname relative to the study directory of the cache file.
        This file holds information (directories in current design) on model simulation directories.
        :return:  filename relative to studyDir
        """

        fileStore = self.Config['begin'].get('studyCacheFile', 'cache_file.json')
        if fileStore is None: fileStore = 'cache_file.json'
        return fileStore



    def targets(self,obsNames=None, scale=False ):
        """
        Get the target values for specific obs names
        :param obsNames: optional list of observations to use
        :param scale: optional if True (default is False) scale target values by scaling
        :return: target values as a pandas series
        """
        if obsNames is None:  obsNames=self.obsNames()
        tvalues = pd.Series([self.Config['targets'].get(k,np.nan) for k in obsNames], index=obsNames)
        if scale: tvalues=tvalues*self.scales(obsNames=obsNames)
        return tvalues.rename(self.name())

    def constraintName(self):
        """
        Extract the name of the constraint (if any) variable from the configuration file

        :return: the name of the constraint variable
        """
        return self.Config['study'].get("constraintName")

    def constraintTarget(self, constraintName=None, scale=False):
        """
        extract the value for the constraint variable target returning it as a pandas series.
        :param constraintName (optional) if not specified then constraintName() will be used
        :param scale (optional; default False) if True then constraint value will be appropriately scaled.
        :return: constraint target as a pandas series
        """

        if constraintName is None: constraintName=self.constraintName() # get constraint name
        return self.targets(obsNames=[constraintName],scale=scale) # wrap name as list and use targets method to get value

    def scales(self,obsNames=None):
        """
        Get the scales for specified obsnamaes
        :param obsNames: optional list of observations to use
        :return: scales as a pandas series
        """
        scalings=self.Config.get('scalings',{})
        if obsNames is None: obsNames= self.obsNames()
        scales = pd.Series([scalings.get(k,1.0) for k in obsNames], index=obsNames).rename(self.name())
        # get scalings -- if not defined set to 1.
        return scales

    def maxFails(self):
        """

        :return: the maximum number of fails allowed. If nothign set then return 0.
        """
        optimise=self.Config.get('optimise',{})

        maxFails=optimise.get('maxFails',0)
        if maxFails is None: maxFails=0
        return maxFails


    def Covariances(self, obsNames=None, trace=False, dirRewrite=None, scale= False, constraint=False):
        """
        If CovObsErr and CovIntVar are both specified then CovTotal will be computed from
        CovObsErr+2*CovIntVar overwriting the value of CovTotal that may have been specified.
        Unspecified values will be set equal to None.

        :param obsNames: Optional List of observations wanted and in order expected.
        :param trace: optional with default False. If True then additional output will be generated.
        :param dirRewrite: optional ith default None. If set then rewrite directory names used in readCovariances.
        :param scale: if set true (default is false) then covariances are scaled by scaling factors derived from self.scales()
        :param constraint: is set to True  (default is False) then add constraint weighting into Covariances
        :return: a dictionary containing CovTotal,CovIntVar, CovObsErr --  the covariance matrices. None if not present.
        Could be modified to cache covariance matrices to save re-reading but not bothering..
        """
        keys=['CovTotal','CovIntVar','CovObsErr']
        if obsNames is None: obsNames=self.obsNames()
        cov={} # empty dict to return things in
        covInfo=self.Config['study']['covariance']
        # extract the covariance matrix and optionally diagonalise it.

        for k in keys:
            fname=covInfo.get(k,None)
            if fname is not None: # specified in the configuration file
                cov[k] = self.readCovariances(fname, obsNames=obsNames, trace=trace, dirRewrite=dirRewrite)
                cov[k+"File"]=fname # store the filename
                if cov[k] is not None: # got some thing to further process
                    if covInfo.get(k+"Diagonalise",False): # want to diagonalise the covariance
                        # minor pain is that np.diag returns a numpy array so we have to remake the DataFrame
                        cov[k]=pd.DataFrame(np.diag(np.diag(cov[k])),index=obsNames,columns=obsNames,dtype=float)
                        if trace: print "Diagonalising "+k

        # make total covaraince from CovIntVar and CovObsErr if both are defined.
        if cov.get('CovIntVar') is not None and cov.get('CovObsErr') is not None:  # if key not defined will "get" None
            k='CovTotal'
            cov[k]=cov['CovObsErr']+2.0*cov['CovIntVar']
            cov[k+'_info']='CovTotal generated from CovObsErr and CovIntVar'
            if trace: print "Computing CovTotal from CovObsErr and CovIntVar"
            if covInfo.get(k+"Diagonalise",False):
                if trace: print "Diagonalising "+k
                cov[k]=pd.DataFrame(np.diag(np.diag(cov['CovTotal'])),index=obsNames,columns=obsNames)
                # diagonalise total covariance if requested.
        # check have total covarince and raiseError if not
        if cov['CovTotal'] is None:
            print "No covtotal found for totalFile=",covInfo.get(k,None)
            raise ValueError

        # apply constraint
        if constraint: # want to have constraint wrapped in to total covariance
            k='CovTotal'
            consValue = self.optimise()['mu']
            consName = self.constraintName()
            obsNames.append(consName)
            c=np.pad(cov[k],(0,1),'constant') # pad with zeros
            c[-1,-1]=1.0/consValue
            cov[k]=pd.DataFrame(c,index=obsNames,columns=obsNames,dtype=float)
        # scale data
        if scale:
            scales=self.scales(obsNames=obsNames)
            cov_scale=pd.DataFrame(np.outer(scales,scales),index=scales.index,columns=scales.index)
            for k in keys:
                if k in cov and cov[k] is not None:
                    cov[k]=cov[k]*cov_scale
                    if trace: print "Scaling "+k




        return cov



    def steps(self, paramNames=None):
        """
        Compute perturbation  for all parameters supplied. If value specified use that. If not use 10% of the range.
        Quasi-scientific in that 10% of range is science choice but needs knowledge of the structure of the JSON file
             so in this module.
        :param paramNames -- optional the parameter names for step sizes. If not defined uses self.paramNames() to work
                them out
        :return: the step sizes for the parameters as a pandas Series.
        """

        if paramNames is None: paramNames=self.paramNames()

        param=self.paramRanges(paramNames=paramNames)
        defaultStep = 0.1*(param.loc['maxParam',:] -param.loc['minParam',:]) # 10% of range and default cases
        pert={}
        for p in paramNames:
            pert[p]=self.Config['steps'].get(p,defaultStep.loc[p])

        perturbation = pd.Series(pert)[paramNames] # make sure in the same order.
        return perturbation.astype(float).rename(self.name())


    def rewriteDir(dir, dir_rewrite):
        """
         rewrite dir using keys in dir_rewrite
        :param dir: list or string of directories to be rewritten
        :param dir_rewrite:  dictionary of files to be rewritten
        :return: rewritten dict
        """
        if isinstance(dir,(str, unicode)): # this is a string
            result=dir
            for k in dir_rewrite.keys(): # iterate over keys
                if k in dir: # got some text to rewrite.
                    result = dir.replace(k,dir_rewrite[k])
                    continue # exit all processing
        else: # it is a list..
            result=[] # output list
            for s in dir: # iterate over list of strings to rewrite
                rs=s
                for k in dir_rewrite.keys(): # iterate over keys
                    if k in s: # got it
                        rs=s.replace(k,dir_rewrite[k])
                        next # no more replacement.
                result.append(rs) # append the rewritten text.


        return result # return the result escaping the loop



    def readCovariances(self,covFile,obsNames=None, trace=False, dirRewrite=None):
        """
        :param covFile: Filename for the covariance matrix. Env variables and ~ expanded
        :param olist: (optional) List of Observations wanted from covariance file
        :param trace: (optional) if set True then some handy trace info will be printed out.
        :param dirRewrite (optional) if set to something then the first key in dirRewrite that matches in covFile
              will be replaced with the element.
        :return: cov -- a covariance matrix sub-sampled to the observations

        Returns a covariance matrix from file optionally sub-sampling to named observations.
        Note if obsName is not specified ordering will be as in the file.
        """

        use_covFile=os.path.expanduser(os.path.expandvars(covFile))
        if dirRewrite is not None:
            use_covFile=self.rewriteDir(use_covFile, dirRewrite)
            if trace: print "use_covFile is ",use_covFile

        if not(os.path.isfile(use_covFile) and os.access(use_covFile, os.R_OK)):
            print "Failed to read ",use_covFile
            raise IOError
        # now read in covariance file
        cov=pd.read_csv(use_covFile) # read the data
        cov.set_index(cov.columns, drop=False, inplace=True,
                                verify_integrity=True) # provide index
        if obsNames is not None: # deal with olist
            cov=cov.loc[obsNames,obsNames] # extract the values comparing to olist
            expect_shape= (len(obsNames),len(obsNames))
            if cov.shape != expect_shape: # trigger error if missing
                print "Sub-sampled covariance shape = ",cov.shape, "And expected = ",expect_shape
                raise ValueError

        return cov


    def postProcessScript(self):
        """

        :return: the full path for the postprocessing script
        """
        ppScript=self.Config['postProcess'].get("script",
                                  "$OPTCLIMTOP/obs_in_nc/comp_obs.py") # get PostProcessScript
        ppScript=os.path.expanduser(os.path.expandvars(ppScript)) # expand shell variables and home
        return ppScript

    def postProcessOutput(self):
        """

        :return: relative  path for output from post processing script. Path is taken relative to model directory

        """

        ppOutput=self.Config['postProcess'].get("outputPath","observations.nc")
        return  ppOutput

    def referenceConfig(self,studyDir= None):
        """
        :param studyDir -- where study is. Default will be to use current working directory
        :return: full path to the reference configuration of model being used
        """
        if studyDir is None: studyDir=os.getcwd()
        modelConfigDir = getDefault(self.Config['study'],'referenceModelDirectory',os.path.join(studyDir,"start"))
        # and now expand home directories and env variables
        modelConfigDir = os.path.expanduser(os.path.expandvars(modelConfigDir))
        return modelConfigDir

    def optimise(self):
        """
        Extract and package all optimistion information into one directory
        :return: a dict
        """

        return self.Config['optimise']

    def fixedParams(self):
        """
        :return: a dict of all the fixed parameters. All names ending _comment or called comment will be excluded 
        """
        import copy

        fix=copy.copy(self.getv('Parameters').get('fixedParams',None))
        if fix is None: return fix # nothing found so return None
        # remove comment and _comment -- I suspect this would be good to do at higher level.
        for k in fix.keys():
            if k == 'comment' or k.endswith('_comment'):
                del fix[k] # delete the unwanted key

        return fix

    def runCode(self):
        """
        
        :return: the runCode (or None) 
        """
        return self.getv("runCode",None)

    def runTime(self):
        """
        
        :return: the run time (or None)
        """
        return self.getv("runTime",None)

    def GNgetset(self,name,variable=None):
        """
        Common method to set/get stuff
        :return:
        """

        GNinfo = self.getv('GNinfo', None)

        if GNinfo is None:  # no GNinfo so create it.
            GNinfo = collections.OrderedDict()
            GNinfo['comment'] = 'Gauss-Newton Algorithm information'

        if variable is None:
            variable = GNinfo.get(name, None)
            if variable is None: return None  # no variable so return None
            variable = np.array(variable)  # extract variable from GNinfo and convert to numpy array
        else:  # got variable so put it in the GNinfo
            GNinfo[name] = variable.tolist()
            self.setv('GNinfo', GNinfo)  # store it.

        return variable

    def GNjacobian(self,jacobian=None, normalise=False, constraint=False):
        """
        Return the Jacobian array as an xarray DataArray
        :param jacobian: (default None). If not None should be a numpy 3D array which will be stored in the config
        :return: xarray version of the Jacobian array
        """

        jacobian=self.GNgetset('jacobian',jacobian)
        if jacobian is None: return None
        # jacobian by default includes constraint. if constraint is False remove it.
        if not constraint:
            jacobian=jacobian[...,0:-1]


        paramNames = self.paramNames()
        obsNames = self.obsNames(add_constraint=constraint)
        iter = np.arange(0,jacobian.shape[0])
        name=self.name()
        jacobian = xarray.DataArray(jacobian,
                                    coords={'Iteration': iter, 'Parameter': paramNames, 'Observation': obsNames},
                                    dims=['Iteration', 'Parameter', 'Observation'], name=name)

        # want to normalise ?
        if normalise:
            rng = self.paramRanges(paramNames=jacobian.Parameter.values)
            rng = xarray.DataArray(rng.loc['rangeParam'], {'Parameter': rng.columns},dims=['Parameter'])
            jacobian = jacobian*rng

        return jacobian

    def GNhessian(self, hessian=None):
        """
        Return the Hessian array as an xarray DataArray
        :param hessian: (default None). If not None should be a numpy 3D array which will be stored in the config
        :return: xarray version of the Hessian array
        """

        hessian = self.GNgetset('hessian', hessian)
        if hessian is None: return None
        paramNames = self.paramNames()
        iter = np.arange(0, hessian.shape[0])
        name = self.name()
        hessian = xarray.DataArray(hessian,
                                    coords={'Iteration': iter, 'Parameter': paramNames, 'Parameter_2': paramNames},
                                    dims=['Iteration', 'Parameter', 'Parameter_2'], name=name)

        return hessian


    def GNparams(self,params=None):

        """
        Return (and optionally set) the best parameter values as an xarray from the GN optimisation.
        :param params: (default None). If not None should be a 2D numpy array which will be stored in the config.
        :return: parameter array as xarray.DataArray
        """

        params = self.GNgetset('params', params)
        if params is None: return None
        paramNames = self.paramNames()
        iterCount = np.arange(0, params.shape[0])
        name = self.name()
        params = xarray.DataArray(params, coords={'Iteration': iterCount, 'Parameter': paramNames},
                                  dims=['Iteration', 'Parameter'], name=name)

        return params


    def GNcost(self,cost=None):
        """
        Return (and optionally set) the cost values as a pandas Series from the GN optimisation.
        :param cost:  (default None)). If not None should a 1D numpy array which will be stored in the config.
        :return: cost as a pandas.Series
        """
        cost = self.GNgetset('cost', cost)
        if cost is None: return None
        iterCount = np.arange(0, cost.shape[0])
        name = self.name()
        cost = pd.Series(cost, index=iterCount,  name=name)
        cost.index.rename('Iteration', inplace=True)

        return  cost

    def GNalpha(self,alpha=None):
        """
        Return (and optionally set) the alpha values as a pandas Series from the N optimisation.
        :param alpha:  (default None)). If not None should a 1D numpy array which will be stored in the config.
        :return: alpha as a pandas.Series
        """
        alpha = self.GNgetset('alpha', alpha)
        if alpha is None: return None
        iterCount = np.arange(0, alpha.shape[0])
        name = self.name()
        alpha = pd.Series(alpha, index=iterCount,  name=name)
        alpha.index.rename('Iteration', inplace=True)

        return  alpha

    def GNparamErrCovar(self,normalise=False,constraint=True, Jac=None):
        """
        Compute the covariance for parameter error.
        Theory is that (to first order) \vect{\delta O}= \matrix{J}\vect{\delta p}
          where \delta O are perturbations in observations, \delta p are perturbations to parameters and J is the Jacobian.
          Then multiple by J^+ (the pseudo-inverse) to give \vect{\delta p} = \matrix{J}^+  \vect{\delta O}.
        Then the covariance is \matrix{J}^+ C {\matrix{J}^+}^T

         or alternatively with covariance...
         P=(J^TC^{-1}J)^{-1}J^TC^{-1}
        :param normalise (default = True) -- compute the normalised covariance error (fraction of range)
        :return: covariance of parameters
        """
        import numpy.linalg as linalg
        invFn=linalg.inv
        # need both the Jacobian and Covariance matrix.
        covar = self.Covariances(constraint=constraint)
        covar = covar['CovTotal'].values

        invCov=invFn(covar)
        if Jac is None: # compute Jacobian
            Jac = self.GNjacobian(normalise=normalise, constraint=constraint).isel(Iteration=-1).T # extract the last Jacobian
        P = invFn(Jac.values.T.dot(invCov).dot(Jac.values)).dot(Jac.values.T).dot(invCov) # transform matrix.
        #paramCovar=linalg.inv(Jac.T.dot(invCov).dot(Jac))*Jac.T*

        #JacPinv = linalg.pinv(Jac.values)
        #paramCovar= JacPinv.dot(covar).dot(JacPinv.T)
        paramCovar = P.dot(covar).dot(P.T)
        # now wrap paramCovar up as a dataframe.
        paramCovar = pd.DataFrame(paramCovar,index=Jac.Parameter,columns=Jac.Parameter)

        return paramCovar

    def optimumParams(self, paramNames=None, scale=False, **kwargs):
        """
        Set/get the optimum parameters.  (VN1 variant)
        :param values: default None -- if set then parameter values in configuration gets updated with values
        :param paramNames -- return specified parameters.
        :param scale (default False). If True scale parameters over range (0 is min, 1 is max)
        :return: values as pandas series.
        """

        # TODO merge this with vn2 code which could be done using  if self.version >= 2: etc

        if len(kwargs) > 0:  # set the values
            self.Config['study']['optimumParams'] = kwargs
        if paramNames is None:  paramNames = self.paramNames()
        values = self.Config['study'].get('optimumParams', None)
        values = {k: values[k] for k in paramNames}
        if values is not None:
            values = pd.Series(values)  # wrap it as a pandas series.
        if scale:
            range = self.paramRanges(paramNames=paramNames)  # get param range
            values = (values-range.loc['minParam',:])/range.loc['rangeParam',:]

        return values

    def GNsimulatedObs(self,set=None, obsNames=None):
        """
        Set/get the simulated observations from the best value in each iteration of the Gauss Newton algorithm.
        :param obsNames: Names of observations wanted
        :param set: a pandas array or xarray (or anything that hass a to_dict method) that will set values of the observations.
          Should only be the best values on each iteration. Perhaps this routne should have all values... 
        :return: Observations as a pandas array.

        """

        if set is not  None: # set the value
            # check we've got GNinfo
            GNinfo=self.getv('GNinfo',collections.OrderedDict()) # should modify to use GNgetset but that assumes numpy array.
            GNinfo['SimulatedObs'] = set.to_dict()
            self.setv('GNinfo',GNinfo) # and set it.
        # get the value.
        
        result = pd.DataFrame(self.getv('GNinfo')['SimulatedObs']) # should trigger an error if the obs does not exist.
        name=self.name()
        # convert if needed from unicode.
        try:
            name=name.encode('utf8')
        except TypeError:
            pass

        #result.rename(name,inplace=True) # rename it.
        
        if obsNames is None:
            obsNames=self.obsNames()
        result=result.loc[:,obsNames]
        
        return result # and return the result.
    
    
    def GNoptimumObs(self,obsNames=None):
        """
        :param obsNames -- the observations to extract
        :return: Optimum observations as a pandas series. 
        """
        
        obs=self.GNsimulatedObs(obsNames=obsNames)
        return obs.iloc[-1,:] # return the last obs. 
        
    # stuff for DFOLS

    def DFOLSinfo(self,diagnosticInfo=None):
        """

        :param diagnosticInfo: (optional -- defult None). If set add ths to the confguration
        :return: the vlues as a  dataframe
        """
        DFOLSinfo = self.getv('DFOLS',
                              collections.OrderedDict())  # should modify to use GNgetset but that assumes numpy array.
        if diagnosticInfo is not None:
            DFOLSinfo['diagnostic'] = diagnosticInfo.to_json(orient='split')

        self.setv("DFOLS",DFOLSinfo) # store the diagnostic info.

        return pd.read_json(DFOLSinfo.get('diagnostic'),orient='split')

    def simObs(self,simObs=None, best = False):
        """

        :param simObs: optional -- default None and if passed should a pandas dataframe
        :param best: optional-- defaultFalse. If true return the obs for the best iteration
        :return:  dataframe of simulated observations
        """
        if simObs is not None:
            # add fake index to preserve index.

            self.setv('simObs',simObs.to_json(orient='index'))

        #print self.get('simObs')
        sObs = pd.read_json(self.getv('simObs'),orient='index')
        if best:
            sObs = sObs.loc[self.bestEval,:]

        return  sObs


    def parameters(self,parameters=None,best=False,normalise=False):
        """

        :param parameters: optional -- default None and if passed should a pandas dataframe
        : param (optional) normalise -- default False. If True the parameters 
           returned are normalised to 0 (min allowed) to 1 (max allowed)
        :return:  dataframe of simulated observations
        """
        if parameters is not None:
            self.setv('parameters',parameters.to_json(orient='index'))


        params = pd.read_json(self.getv('parameters'),orient='index')
        if best:
            params = params.loc[self.bestEval,:]

        if normalise:
            prange = self.paramRanges(paramNames=params.columns)  # get param range
            params = (params-prange.loc['minParam',:])/prange.loc['rangeParam',:]

        return  params


    def cost(self,cost=None,best =False):
        """

        :param cost: (optional) the cost for each model evaluation as a pandas series.
        :param best (optional) -- if True return only the best value
        :return: cost as a pandas series.
        """

        if cost is not None:
            self.setv('costd',cost.to_json())
        cost = pd.read_json(self.costd,typ='series')
        cost.name = 'Cost'
        if best:
            cost=cost.loc[self.bestEval]

        return cost
    def directories(self,directories=None, best=False):
        """

        :param directories: (optional -- default = None). The directory (currently full path) where each simulation was run as a pandas series
        :param best (optional) -- if True return only the best value
        :return: directories as a pandas series
        """
        if directories is not None:
            self.setv('dirValues',directories.to_json())

        directories = self.getv('dirValues',[])
        if len(directories) is 0: return pd.Series()
        directories = pd.read_json(directories,typ='series')

        if best:
            directories = directories[self.bestEval]

        return directories

class OptClimConfigVn2(OptClimConfig):
    """
    Version 2 of OptClimConfig -- modify OptClimConfig methods and see OptClimConfig.__init__() for  
    """

    # NB __init__ method is just the superclasss OptClimConfig) __init__ method.

    def paramNames(self):
        """
        :return: a list of parameter names from the configuration files
        """
        keys=self.Config['Parameters']['initParams'].keys()  # return a copy of the list.
        # remove comment keys
        keys = [k for k in keys if 'comment' not in k]
        return  keys

    def standardParams(self, paramNames=None):
        """
        Extract standard parameter values for study
        :param paramNames: Optional names of parameters to use.
        :return: pandas series
        """
        if paramNames is None:   paramNames = self.paramNames()
        svalues = pd.Series([self.Config['Parameters']['defaultParams'].get(k,np.nan) for k in paramNames],
                            index=paramNames)

        return svalues.rename(self.name())

    def paramRanges(self, paramNames= None):
        """
        :param paramNames -- a list of the parameters to extract ranges for.
        If not supplied the paramNames method will be used.
        :return: a pandas array with rows names minParam, maxParam, rangeParam
        """
        if paramNames is None: paramNames=self.paramNames()
        param = pd.DataFrame(self.Config['Parameters']['minmax'],
            index=['minParam', 'maxParam'])
        param = param.loc[:,paramNames] # just keep the parameters we want getting rid of comments etc in the JSON file
        param = param.astype(float)
        param.loc['rangeParam',:] = param.loc['maxParam',:]-param.loc['minParam',:] # compute range
        return param

    def beginParam(self,paramNames = None,scale=False):
        #TODO -- make this a set/get (if values passed in then set them).
        """
        get the begin parameter values for the study. These are specified in the JSON file in begin block
        Any values not specified use the standard values
        :param paramNames: Optional names of parameters to use.
        :param scale (default False). If True scale parameters by their range so 0 is minimum and 1 is maximum
        :return: pandas series of begin parameter values.
        """
        if paramNames is None:  paramNames=self.paramNames()
        begin = {} # empty dict
        standard=self.standardParam(paramNames=paramNames)
        scaleRange=self.Config['Parameters'].get("initScale") # want to scale ranges?
        range=self.paramRanges(paramNames=paramNames) # get param range

        for p in paramNames: # list below is probably rather slow and could be sped up!
            begin[p] = self.Config['Parameters']['initParams'].get(p,standard.get(p))
            if begin[p] is None: begin[p]=standard[p]
            if scaleRange: # values are specified as 0-1
                begin[p]=begin[p]*range.loc['rangeParam',p]+range.loc['minParam',p]
            if scale: # want to return them in range 0-1
                    begin[p]=(begin[p]-range.loc['minParam',p])/range.loc['rangeParam',p]



        begin=pd.Series(begin, dtype=float)[paramNames] # order in the same way for everything.

        # verify values are within range
        if scale: L=begin.gt(1.0) | begin.lt(0.0)
        else: L = range.loc['maxParam',:].lt(begin) | begin.lt(range.loc['minParam',:])

        if np.any(L):
            print "L  \n",L
            print "begin: \n",begin
            print "range: \n",range
            print "Parameters out of range",begin[L].index
            raise ValueError("Parameters out of range: ")

        return begin.astype(float).rename(self.name())

    def optimumParams(self,paramNames = None, scale=False, **kwargs):
        """
        Set/get the optimum parameters.
        :param normalise (default False). If True then normalise parameters.
        :param values: default None -- if set then parameter values in configuration gets updated with values
        :param paramNames -- name of parameters/
        :return: values as pandas series.
        """

        if paramNames is None:  paramNames = self.paramNames()

        if len(kwargs) >0 : # set the values
            self.Config['Parameters']['optimumParams']=kwargs
            # add defaults for ones we have not got.
        default=self.standardParam(paramNames=paramNames)

        values=self.Config['Parameters'].get('optimumParams',None)
        if values is None:
            return values # no optimum values so return none

        for k in paramNames:
            values[k]=self.Config['Parameters']['optimumParams'].get(k,default[k])
                    # set it to std value if we don't have it...Note will trigger error if default not got value


        ##values = {k: values.get(k,None) for k in paramNames} # extract the requested params setting values to None if we don't have them.
        values = pd.Series(values)[paramNames] # wrap it as a pandas series and order it appropriately.

        if scale:
            range = self.paramRanges(paramNames=paramNames)  # get param range
            values = (values-range.loc['minParam',:])/range.loc['rangeParam',:]


        return values.rename(self.name())


    def steps(self, paramNames=None):
        """
        Compute perturbation  for all parameters supplied. If value specified use that. If not use 10% of the range.
        Quasi-scientific in that 10% of range is science choice but needs knowledge of the structure of the JSON file
             so in this module.
        :param paramNames -- optional the parameter names for step sizes. If not defined uses self.paramNames() to work
                them out
        :return: the step sizes for the parameters as a pandas Series.
        """

        if paramNames is None: paramNames=self.paramNames()

        param=self.paramRanges(paramNames=paramNames)
        defaultStep = 0.1*param.loc['rangeParam',:] # 10% of range and default cases
        pert={}
        for p in paramNames:
            pert[p]=self.Config['Parameters']['steps'].get(p,defaultStep.loc[p])
            if pert[p] is None: pert[p]=defaultStep.loc[p]
        perturbation = pd.Series(pert)[paramNames] # make sure in the same order.
        return perturbation.astype(float).rename(self.name())

    def cacheFile(self):
        """
        Get the pathname relative to the study directory of the cache file.
        This file holds information (directories in current design) on model simulation directories.
        :return:  filename relative to studyDir
        """

        fileStore = self.Config.get('studyCacheFile', 'cache_file.json')
        if fileStore is None: fileStore = 'cache_file.json'
        return fileStore

class OptClimConfigVn3(OptClimConfigVn2):
    """
    Vn3 of OptClimConfig. Currently does nothing but is a place to grab things for next upgrade. 
    1) have runTime() and runCode() methods only work with runInfo block -- vn1/2 methods work with names in top level dict. 
    """


