"""
test_modelSimulation: test cases for modelSimulation methods.
This test routine needs OPTCLIMTOP specified sensibly. 
"""

import collections
import math
import os
import shutil
import unittest
import f90nml
import pandas as pd

from OptClimVn2 import ModelSimulation
from HadCM3 import  iceAlbedo


__author__ = 'stett2'


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
        # retry
        func(path)

class testModelSimulation(unittest.TestCase):
    """
    Test cases for modelSimulation. There should be one for every method in ModelSimulation.

    """
    def setUp(self):
        """
        Standard setup for all test cases
        :return:
        """
        import time
        parameters=collections.OrderedDict()
        parameters['one']=1.0
        parameters['two']=2.0
        testDir='$OPTCLIMTOP/test_out'
        refDir='$OPTCLIMTOP/test_in'
        testDir='test_out'
        refDir='test_in'
        self.dirPath=testDir
        self.refPath=refDir

        testDir=os.path.expanduser(os.path.expandvars(testDir))
        refDir=os.path.expandvars(os.path.expanduser(refDir))
        shutil.rmtree(testDir,onerror=errorRemoveReadonly)
        if os.path.exists(testDir):
            os.rmdir(testDir)


        self.model= ModelSimulation.ModelSimulation(testDir,
                                                    name='test', create=True, refDirPath=os.path.join(refDir,'start'),
                                                    ppExePath='postProcess.sh',
                                                    ppOutputFile='obs.nc',
                                                    parameters=parameters,
                                                    obsNames=['temp@500_nhx','temp@500_tropics','temp@500_shx'],
                                                    verbose=True)
        shutil.copy(os.path.join(refDir,'01_GN','h0101','observables.nc'),
                    os.path.join(testDir,'obs.nc')) # copy over a netcdf file of observations.


    def test_Init(self):
        """
        Test init methods works.
        Probably over-kill and could be a pain if internal details change.
        But idea is that public methods all get a work out when modelSimulation initialised.
        :return:
        """
        # using implicit run of setup.
        expectObs=collections.OrderedDict()
        for k in ['temp@500_nhx','temp@500_tropics','temp@500_shx']:expectObs[k]=None

        expectParm=collections.OrderedDict()
        expectParm['one']=1.0
        expectParm['two']=2.0

        self.assertEqual(self.model.get(['name']),'test')
        self.assertEqual(self.model.get(['ppExePath']), 'postProcess.sh')
        self.assertEqual(self.model.get(['observations']),expectObs)
        self.assertEqual(self.model.get(['parameters']),expectParm)
        self.assertEqual(self.model.get(['ppOutputFile']),'obs.nc')
        self.assertListEqual(self.model.get(['observations']).keys(),expectObs.keys())
        # test that read works. Works means no failures and have observations..

        m= ModelSimulation.ModelSimulation(self.dirPath, verbose=True)
        # observations should have changed but nothing else.
        self.assertEqual(m.get(['name']),'test')
        self.assertEqual(m.get(['ppExePath']), 'postProcess.sh')
        #self.assertEqual(m.get(['observations']),expectObs)
        self.assertEqual(m.get(['parameters']),expectParm)
        #self.assertEqual(self.model.config['refDir'], None)
        self.assertEqual(m.get(['ppOutputFile']),'obs.nc')
        self.assertListEqual(m.getObs().keys(),expectObs.keys())
        self.assertNotEqual(m.getObs(),expectObs)
        # # updating parameters should trigger an error as we are read only
        # with self.assertRaises(Exception):
        #    m.setParams({'aaa':99.9},addParam=True)

        ## now  to test that update works.
        m= ModelSimulation.ModelSimulation(self.dirPath, verbose=True, update=True)
        # observations should have changed but nothing else.
        self.assertEqual(m.get(['name']),'test')
        self.assertEqual(m.get(['ppExePath']), 'postProcess.sh')
        #self.assertEqual(m.get(['observations']),expectObs)
        self.assertEqual(m.get(['parameters']),expectParm)
        #self.assertEqual(self.model.config['refDir'], None)
        self.assertEqual(m.get(['ppOutputFile']),'obs.nc')
        self.assertListEqual(m.getObs().keys(),expectObs.keys())
        self.assertNotEqual(m.getObs(),expectObs)
        m.setParams({'aaa':99.9},addParam=True)




    def test_get(self):
        """
        Test that get works
        :return:
        """
        self.model.set({'fbf':1.0,'fbg':2.0},write=False) # set value
        self.assertEqual(self.model.get('fbf'),1.0)
        self.assertEqual(self.model.get(['fbf','fbg']),[1.0,2.0])

    def test_getObs(self):
        """
        Test that getObs works. Not really needed as init test checks this but we'll do it anyhow.
        :return:
        """
        self.assertEqual(self.model.getObs(),self.model.get('observations'))


    def test_getParams(self):
        """
        Test that getParam works. 
        :return:
        """
        self.assertEqual(self.model.getParams(),self.model.get('parameters'))
        # test if we specify parameters it works
        self.assertEqual(self.model.getParams(params=['one','two']), self.model.get('parameters'))



    def test_readModelSimulation(self):
        """
        Test that readModelSimulation works.  Not really needed as init test checks this works
          do anyway.
        :return:
        """
        # compare exisiting config with fresh one we read in.
        self.model.readObs() # should prob be init.
        m= ModelSimulation.ModelSimulation(self.model.dirPath, verbose=True)
        self.assertEqual(m.get(),self.model.get())
        # verify we fail to modify
        with self.assertRaises(Exception):
            m.set({'someInfo',{'silly':2,'silly_more':3}})
        # make an updatable version.
        m = ModelSimulation.ModelSimulation(self.model.dirPath, update=True, verbose=True)
        self.assertEqual(m.get(), self.model.get())
        # verify can  modify
        m.set({'someInfo':{'silly': 2, 'silly_more': 3}})  # this will get written through but next read should overwrite.
        # now config should be different from ref case. Read it in and verify
        m = ModelSimulation.ModelSimulation(self.model.dirPath, verbose=True)
        self.assertNotEqual(m.get(), self.model.get())

    def test_readObs(self):
        """
        Test that readObs works
        :return:
        """
        # do by changing the observations and then rereading..
        obs=self.model.getObs()
        self.assertEqual(type(obs),collections.OrderedDict)
        obs['rh@500_nhx']=False
        self.model.set({'observations':obs},write=False)
        self.model.readObs(verbose=True)
        mobs=self.model.getObs()
        self.assertEqual(mobs.keys(),obs.keys())
        self.assertNotEqual(mobs,obs)
        # test that json failure works...
        dir=self.model.dirPath
        # remove the obs.nc file and copy OPTCLIMTOP/Configurations/example.json to fail.json
        os.remove(os.path.join(dir,'obs.nc'))
        #refJson=os.path.join(os.path.expandvars(os.path.join('$OPTCLIMTOP','Configurations','example.json')))
        refJson=os.path.join(os.path.expandvars(os.path.join('Configurations','example.json')))
        reserveFile=os.path.join(dir,'fail.json')
        shutil.copy(refJson,reserveFile)
        import StudyConfig
        config=StudyConfig.readConfig(reserveFile,ordered=True)
        expectObs = config.simulatedObservations # extract obs
        varsWant=self.model.get(['observations']).keys()
        eObs=collections.OrderedDict()
        for k in varsWant:
            eObs[k]=expectObs.get(k)
        self.model.readObs(verbose=True)
        obs = self.model.getObs(verbose=True)
        self.assertEqual(eObs,obs)




    def test_set(self):
        """
        Test that set works by setting a value and then getting it back.
        set modified config file so will check that config file has been modified.
        :return:
        """
        # test without writing..
        self.model.set({"fbg":1.0}, write=False)
        self.assertEqual(self.model.get('fbg'),1.0)
        # then with writing
        m= ModelSimulation.ModelSimulation(self.model.dirPath, update=True)
        # get config file.
        config=m.readConfig()
        m.set({"fbg": 1.0})
        self.assertEqual(m.get('fbg'), 1.0)
        self.assertNotEqual(m.get(),config) # on disk config should have changed...


    def test_setParams(self):
        """
        Test that setParams works
        :return:
        """
        param=collections.OrderedDict()
        param['p1']=1
        param['p2']=2
        self.model.setParams(param,write=False)
        self.assertEqual(self.model.get('parameters'),param)
        # verify addParam works
        self.model.setParams({'p3':3}, write=False,addParam=True, verbose=True)
        param['p3']=3
        self.assertEqual(self.model.get('parameters'), param)

        # check a pandas series work
        params=pd.Series([1,2,3],index=['p1','p2','p3'])
        self.model.setParams(params,verbose=True,write=False)
        self.assertEqual(self.model.get('parameters').to_dict(), param)


    def test_genVarToNameList(self):
        """
        Test that genVarToNameList works
        :return:
        """



        self.model.genVarToNameList('RHCRIT',"RHCRIT",'CNTL','ATMOS.NL')
        # check we've got what we think we should have
        expect=[ModelSimulation._namedTupClass(var='RHCRIT',namelist='CNTL',file='ATMOS.NL')]
        self.assertEqual(self.model._convNameList['RHCRIT'],expect)
        # try two var case with list of functions

    def test_registerMetaFn(self):
        """
        Test that registering meta function works.
        :return:
        """
        def fn(x=1.0, inverse=False, namelist=False):
            if inverse:
                return x['fn1'][0]
            else:
                return {'fn1':[x]*19}

        def fn2(x=1.0, inverse=False, namelist=False):
            if inverse:
                return math.sqrt(x['fn2'][0])
            else:
                return {'fn2': [x**2] * 19}


        self.model.registerMetaFn('rhcrit',fn)
        self.model.registerMetaFn('rhcrit2', fn2)
        self.assertEqual(self.model._metaFn['rhcrit'],fn)
        self.assertEqual(self.model._metaFn['rhcrit2'], fn2)

    def test_applyMetaFn(self):
        """
        Test that applying metaFn works
        :return:
        """

        def fn(x=1.0, inverse=False, namelist=False):
            if inverse:
                return x['rhcrit'][0]
            else:
                return {'rhcrit': [x] * 19}

        def fn2(x=1.0, inverse=False, namelist=False):
            if inverse:
                return math.sqrt(x['rhcrit2'][0])
            else:
                return {'rhcrit2': [x ** 2] * 19}


        def fn3(x=2.0, inverse=False, namelist=False):
            if inverse:
                return x['rhcrit3'][2]
            else:
                return {'rhcrit':x+2, 'rhcrit2':x/2, 'rhcrit3':[x]*19}


        self.model.registerMetaFn('rhcrit',fn,verbose=True)
        self.assertEqual(self.model.applyMetaFns(rhcrit=1.0),({'rhcrit':[1.0]*19},['rhcrit']))

        # verify failure happens when we force it!
        with self.assertRaises(Exception):
           print self.model.applyMetaFns(fail=True, rhcrit2=1.0)
        # verify failure *does not* happen when we don't ask for it.
        self.assertEqual(self.model.applyMetaFns(rhcrit2=1.0),({},[]))
        # add  a 2nd fn
        self.model.registerMetaFn('rhcrit2', fn2, verbose=True)
        self.assertEqual(self.model.applyMetaFns(rhcrit=1.0,rhcrit2=2.0), ({'rhcrit': [1.0] * 19,'rhcrit2':[4.0]*19},['rhcrit','rhcrit2']))

        # verify failure happens when we force it!
        with self.assertRaises(Exception):
            print self.model.applyMetaFns(fail=True, rhcrit=1.0,rhcrit3=1.0)

            # verify multiple parameter apply works..
        self.model.registerMetaFn('rhcritAll', fn3, verbose=True)
        self.assertEqual(self.model.applyMetaFns(rhcritAll=1.0),
                         ({'rhcrit': 3.0, 'rhcrit2': 0.5, 'rhcrit3': [1.0] * 19},
                          ['rhcritAll']))








    def test_writeNameList(self):
        """
        Test that namelist creation works
        :return:
        """
        ## copy a trial namelist from configurations.
        infile=os.path.join("test_in","start","CNTLATM")
        outFile = 'CNTLATM'
        outPath=os.path.join("test_out",outFile)
        shutil.copyfile(infile,outPath) # copy test file over
        backFile=outPath+"_nl.bak"
        # now to make configuration and patch outfile
        self.model._readOnly=False # we can write to it.
        self.model.genVarToNameList('VF1',nameListVar='vf1',nameListName='slbc21',nameListFile=outFile)

        self.model.writeNameList(verbose=True, VF1=1.5)
        # should have a backup file
        self.assertEqual(os.path.isfile(backFile),True)
        # read the namelist
        with open(outPath) as nml_file:
            nml = f90nml.read(nml_file)
        self.assertEqual(nml['SLBC21']['vf1'],1.5) # check namelist modified
        # make vf1 change two variables...
        os.remove(outPath)
        os.rename(backFile, outPath)  # move backup file back again
        self.model.genVarToNameList('O2MMR', nameListVar='o2mmr', nameListName='runcnst', nameListFile=outFile)

        self.model.writeNameList( verbose=False, fail=True, VF1=1.5, O2MMR=0.21)
        #
        # should have a backup file
        self.assertEqual(os.path.isfile(backFile), True)
        # read the namelist
        with open(outPath) as nml_file:
            nml = f90nml.read(nml_file)
        self.assertEqual(nml['SLBC21']['vf1'], 1.5)  # check namelist modified
        self.assertEqual(nml['runcnst']['o2mmr'],0.21)
        # add another nl variable.
        self.model.genVarToNameList('RHCRIT', nameListVar='rhcrit', nameListName='runcnst',
                               nameListFile=outFile)
        # generate the namelist.
        rhcv=[0.65,0.7]
        rhcv.extend([0.8]*17) # there is a bug with f90nml version 0.19 which does not overwrite existing parts of the array
        self.model.writeNameList(verbose=True, VF1=1.5, RHCRIT=rhcv )
        self.assertEqual(os.path.isfile(backFile), True)
        # read the namelist
        with open(outPath) as nml_file:
            nml = f90nml.read(nml_file)
        self.assertEqual(nml['SLBC21']['vf1'], 1.5)  # check namelist modified
        self.assertEqual(nml['runcnst']['o2mmr'],0.21)
        self.assertEqual(nml['runcnst']['rhcrit'],rhcv)

    def test_readNameList(self):
        """
        Test cases for readNameListVar
        :return:
        """
        outFile = 'CNTLATM'
        self.model.setReadOnly(False)
        self.model.genVarToNameList('RHCRIT', nameListVar='rhcrit', nameListName='runcnst',
                                    nameListFile=outFile)
        self.model.genVarToNameList('VF1', nameListVar='vf1', nameListName='slbc21', nameListFile=outFile)

        rhcv = [0.65, 0.7]
        rhcv.extend([0.8] * 17)  # there is a bug with f90nml version 0.19 which does not overwrite existing parts of the array
        expect=collections.OrderedDict()
        expect['RHCRIT']=rhcv
        expect['VF1']=1.5
        self.model.writeNameList(verbose=True, VF1=1.5, RHCRIT=rhcv)
        self.model.setReadOnly(True)
        vars=self.model.readNameList( ['RHCRIT', 'VF1'], fail=True)
        self.assertDictEqual(vars,expect)

    def test_readMetaNameList(self):
        """
        test cases for readMetaNameList
        :return:
        """
        # set up meta fn and test it works... HadCM3 as will read the standard file.



        self.model.registerMetaFn('ALPHAM', iceAlbedo)  # register ice fn
        a=self.model.readMetaNameList('ALPHAM', verbose=True)
        self.assertEqual(a,0.5)


if __name__ == "__main__":
    print "Running Test Cases"
    unittest.main() ## actually run the test cases