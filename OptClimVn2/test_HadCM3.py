"""
test_HadCM3: test cases for HadCM3 methods.
This test routine needs OPTCLIMTOP specified sensibly. (ideally so that $OPTCLIMTOP/um45 corresponds to
  where this routine lives..
"""
import collections
import os
import shutil
import unittest
import re
import exceptions
import numpy as np
import pandas as pd

import HadCM3

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

def cmp_lines(path_1, path_2, ignore=None, verbose=False):
    """
    From Stack exchange -- http://stackoverflow.com/questions/23036576/python-compare-two-files-with-different-line-endings
    :param path_1: path to file 1 
    :param path_2: path to file 2
    :param ignore -- list of regep patterns to ignore
    :param verbose: default False -- if True print out lines that don't match.
    :return: True if files the same, False if not.
    """

    if ignore is None: ignore=[] # make ignore an empty list if None is provided.
    l1 = l2 = ' '
    with open(path_1, 'r') as f1, open(path_2, 'r') as f2:
        while l1 != '' and l2 != '':
            l1 = f1.readline()
            l2 = f2.readline()
            skip=False # skip when set true and set so if ignore pattern found
            for p in ignore:
                if re.search(p,l1) or re.search(p,l2):
                    skip=True
                    continue
            if skip: continue # go to next while bit of of looping over files.
            if l1 != l2:
                if verbose:
                    print ">>",l1
                    print "<<",l2
                return False
    return True



class testModelSimulation(unittest.TestCase):
    """
    Test cases for HadCM3. There should be one for every method in HadCM3.

    """
    def setUp(self):
        """
        Setup case
        :return:
        """
        parameters={"CT": 1e-4, "EACF": 0.5,"ENTCOEF": 3.0, "ICE_SIZE": 30e-6,
            "RHCRIT": 0.7, "VF1": 1.0,"CW_LAND": 2e-4,"DYNDIFF": 12.0,"KAY_GWAVE": 2e4,
            "ASYM_LAMBDA": 0.15, "CHARNOCK": 0.012, "G0": 10.0,"Z0FSEA": 1.3e-3,"ALPHAM": 0.5,
                    'START_TIME':[1997, 12, 1], 'RUNID': 'a0101','ASTART':'$MYDUMPS/fred.dmp'}

        testDir = os.path.join('$OPTCLIMTOP','test_out')
        refDir = os.path.join('$OPTCLIMTOP','Configurations','HadAM3_ed3_SL7_15m')
        simObsDir= os.path.join('$OPTCLIMTOP','test_in')
        testDir = 'test_out'
        refDir = os.path.join('Configurations','HadAM3_ed3_SL7_15m')
        simObsDir= 'test_in'
        self.dirPath = testDir
        self.refPath = refDir

        testDir = os.path.expanduser(os.path.expandvars(testDir))
        refDir = os.path.expandvars(os.path.expanduser(refDir))
        simObsDir = os.path.expandvars(os.path.expanduser(simObsDir))
        shutil.rmtree(testDir, onerror=errorRemoveReadonly)
        self.model = HadCM3.HadCM3(testDir, name='a0101', create=True, refDirPath=refDir,
                                   ppExePath='postProcess.sh', SPHERICAL_ICE=False,
                                   ppOutputFile='obs.nc', runTime=1200, runCode='test',
                                   obsNames=['temp@500_nhx', 'temp@500_tropics', 'temp@500_shx'],
                                   verbose=True, **parameters)

        shutil.copy(os.path.join(simObsDir, '01_GN', 'h0101', 'observables.nc'),
            os.path.join(testDir, 'obs.nc'))  # copy over a netcdf file of observations.


    def test_init(self):
        """
        Test init
        :return:
        """
        # possible tests
        # 1) Test that changed name and directory has been changed.
        #  Note name and directory don't need to be the same.
        """
        Test init methods works.
        Probably over-kill and could be a pain if internal details change.
        But idea is that public methods all get a work out when modelSimulation initialised.
        :return:
        """
        # using implicit run of setup.
        expectObs=collections.OrderedDict()
        for k in ['temp@500_nhx','temp@500_tropics','temp@500_shx']:expectObs[k]=None
        expectParam={"CT": 1e-4, "EACF": 0.5,"ENTCOEF": 3.0, "ICE_SIZE": 30e-6,
            "RHCRIT": 0.7, "VF1": 1.0,"CW_LAND": 2e-4,"DYNDIFF": 12.0,"KAY_GWAVE": 2e4,"SPHERICAL_ICE":False,
            "ASYM_LAMBDA": 0.15, "CHARNOCK": 0.012, "G0": 10.0,"Z0FSEA": 1.3e-3,"ALPHAM": 0.5,
            'START_TIME':[1997, 12, 1], 'RUNID': 'a0101','ASTART':'$MYDUMPS/fred.dmp'}
        self.assertEqual(self.model.get(['name']),'a0101')
        self.assertEqual(self.model.get(['ppExePath']), 'postProcess.sh')
        self.assertEqual(self.model.get(['observations']),expectObs)
        self.assertDictEqual(self.model.get(['parameters']),expectParam)
        self.assertEqual(self.model.get(['ppOutputFile']),'obs.nc')
        self.assertListEqual(self.model.get(['observations']).keys(),expectObs.keys())
        # test that read works. Works means no failures and have observations..

        m= HadCM3.HadCM3(self.dirPath, verbose=True)
        # Nothing should have changed except observations have been read in
        self.assertEqual(m.get(['name']),'a0101')
        self.assertEqual(m.get(['ppExePath']), 'postProcess.sh')
        self.assertDictEqual(m.get(['parameters']),expectParam)

        #self.assertEqual(self.model.config['refDir'], None)
        self.assertEqual(m.get(['ppOutputFile']),'obs.nc')
        self.assertListEqual(m.getObs().keys(),expectObs.keys())
        self.assertNotEqual(m.getObs(),expectObs)
        # test that consistency checks work

        with self.assertRaises(exceptions.NameError):
            m = HadCM3.HadCM3(self.dirPath, verbose=True,create=True,AINITIAL='test.dmp',ASTART='test2.dmp')

        with self.assertRaises(exceptions.NameError):
            m = HadCM3.HadCM3(self.dirPath, verbose=True,create=True,OINITIAL='test.dmp',OSTART='test2.dmp')
        # test that read state fails if param provided.

        with self.assertRaises(Exception):
            m = HadCM3.HadCM3(self.dirPath, verbose=True, VF1=2.5)

        ## now  to test that update works.
        m= HadCM3.HadCM3(self.dirPath, verbose=True, update=True, RHCRIT=0.6, RUNID='abcde')
        expectParam['RHCRIT']=0.6
        expectParam['RUNID']='abcde'
        # observations should have changed but nothing else.
        self.assertEqual(m.get(['name']),'a0101')
        self.assertEqual(m.get(['ppExePath']), 'postProcess.sh')
        #self.assertEqual(m.get(['observations']),expectObs)
        self.assertDictEqual(m.get(['parameters']), expectParam)
        #self.assertEqual(self.model.config['refDir'], None)
        self.assertEqual(m.get(['ppOutputFile']),'obs.nc')
        self.assertListEqual(m.getObs().keys(),expectObs.keys())
        self.assertNotEqual(m.getObs(),expectObs)

    def test_readMetaParams(self):
        """
        Test that HadCM3 specific meta functions all work..
        :return:
        """
        expect_values = {"CT": 1e-4, "EACF": 0.5,"ENTCOEF": 3.0, "ICE_SIZE": 30e-6,
            "RHCRIT": 0.7, "VF1": 1.0,"CW_LAND": 2e-4,"DYNDIFF": 12.0,"KAY_GWAVE": 2e4,'SPHERICAL_ICE':False,
            "ASYM_LAMBDA": 0.15, "CHARNOCK": 0.012, "G0": 10.0,"Z0FSEA": 1.3e-3,"ALPHAM": 0.5,'RUN_TARGET':[1,3,21,0,0,0],
            'START_TIME':[1997, 12, 1], 'RUNID': 'a0101',
            'OSTART':'$DATAW/$RUNID.ostart','ASTART':'$MYDUMPS/fred.dmp','AINITIAL':'$MY_DUMPS/xefqda.daj9c10', 'OINITIAL':'$DUMPS/acpdko.daj9c10'}
        got_values={}
        for k in self.model._metaFn.keys(): # test all meta functions
            got_values[k]=self.model.readMetaNameList(k, verbose=True)

        for k,v in got_values.iteritems(): # iterate over meta functions.
            if isinstance(v,basestring):
                self.assertEqual(v,expect_values[k])
            elif isinstance(v,list):
                self.assertEqual(v,expect_values[k])
            else:
                self.assertAlmostEqual(v,expect_values[k],places=2,
                                   msg='Failed to almost compare for %s got %.4g expect %.4g'%(k,v,expect_values[k]))



    def test_setParams(self):
        """
        Test setParams
        :return:
        """
        # will test that can set namelist variables, that setting something that doesn't exist fails.
        self.model.setReadOnly(False) # want to modify model.
        param={'VF1':1.5, 'ICE_SIZE':32e-6, 'ENTCOEF':0.6, 'CT':2e-4, 'ASYM_LAMBDA':0.2,
               'CHARNOCK':0.015, 'G0':20.0, 'Z0FSEA':1.5e-3,'AINITIAL':'fred.dmp','OINITIAL':'james.dmp'}
        metaParam={'KAY_GWAVE':2e4,'ALPHAM':0.65,'CW_LAND':1e-3,'RHCRIT':0.666,'EACF':0.777,
                   'DYNDIFF':11.98,'RUN_TARGET':[2,1,1,0,0,0]}
        un=collections.OrderedDict()
        for k,v in param.iteritems():
            un[k]=v
        expect=un.copy()
        # got problem here.
        for k,v in metaParam.iteritems():
            un[k]=v
            if type(v) == np.ndarray: v = v.round(3)
            expect[k] = v

        self.model.setParams(un,fail=True,verbose=True)
        # verify namelists are as expected.
        vars=self.model.readNameList(expect.keys(), verbose=True,fail=True)

        for k in expect.keys():
            msg = 'Key is %s' % (k)
            print "vars[%s]=%s got %s"%(k,vars[k],expect[k])
            if type(expect[k]) == list:
                self.assertEqual(vars[k],expect[k],msg=msg)
            else:
                self.assertAlmostEqual(expect[k],vars[k],msg=msg)

        # check pd.Series works
        series=pd.Series(un)
        series['VF1']=1.75
        expect['VF1']=1.75
        self.model.setReadOnly(False)  # want to modify model
        self.model.setParams(series,fail=True,verbose=True)
        # verify namelists are as expected.
        vars = self.model.readNameList(expect.keys(), verbose=True, fail=True)

        for k in expect.keys():
            print "vars[%s]=%s got %s" % (k, vars[k],expect[k])
            if type(expect[k]) == list:
                self.assertEqual(vars[k],expect[k],msg=msg)
            else:
                self.assertAlmostEqual(expect[k],vars[k],msg=msg)


    def test_modifyScript(self):
        """
        Test modifyScript produces expected result
        and fails if it tries to modify something already modifies
        :return: 
        """

        # testing modifyScript is tricky. I did use comparision with a file but I think this is difficult to maintain.
        # better to count the number of lines with ## modified at the end.
        modifyStr='## modified *$'
        file=os.path.join(self.model.dirPath,'SCRIPT')
        count=0
        with open(file,'r') as f:
            for line in f:
                if re.search(modifyStr,line): count += 1

        # expected changes are  exprID, jobID, 4 DATA[M,W,U,T] +2 more DATA [M,W]+ MY_DATADIR + mark (5 changes). -- this depends on config.
        # If emailing then will expect more changes.
        expect=14
        # Note config this being tested on has no MY_DATADIR/A
        self.assertEqual(count, expect, 'Expected %d %s got %d'%(expect,modifyStr,count))

        ## two lines below old implementation and commented out
        ##stat=cmp_lines(os.path.join('expectHadCM3','SCRIPT'),file,verbose=True)
        ##self.assertEqual(stat,True,'script files differ')
        self.assertRaises(Exception, self.model.modifyScript)
        #

    def test_modifySubmit(self):
        """
        Test modifySubmit
        :return: 
        """


        # no need to run modifySubmit as already ran when init happens.

        modifyStr = '## modified$'
        file = os.path.join(self.model.dirPath, 'SUBMIT')
        count = 0
        with open(file, 'r') as f:
            for line in f:
                if re.search(modifyStr, line): count += 1
        # expected changes are CJOBN, RUNID, JOBDIR, MY_DATADIR, 2 [NR]RUN+_TIME_LIMIT, ACCOUNT and 5 DATADIR/$RUNID
        expect = 12
        self.assertEqual(count, expect, 'Expected %d %s got %d' % (expect, modifyStr, count))


        self.assertRaises(Exception, self.model.modifySubmit) # run modify a 2nd time and get an error.


    def test_createWorkDir(self):
        """
        Test that createWorkDir worked as expected
        :return: 
        """
        # no need to run createWorkDIr as already ran by init
        # just check it exists and is a dir
        self.assertTrue(os.path.isdir(os.path.join(self.model.dirPath,'W')))

    def test_fixClimFCG(self):
        """
        Test that fixClimFCG works as expects 
         
         converting all CLIM_FCG_YEARS(1,..) to CLIM_FCG_YEARS(:,...)
        :return: 
        """

        # test is that have NO lines in CNTLALL that look like (1,\w[0-9])

        file = os.path.join(self.model.dirPath, 'CNTLATM')
        with open(file, 'r') as f:
            for line in f:
                if re.search('\(1\w*,\w*[0-9]*\w\)=\w*[0-9]*,', line): self.fail("Found bad namelist %s"%line)

    def test_submit(self):
        """
        Test the submit method works -- returns sensible path.
        Rather trivial test.. 
        :return: 
        """
        p=self.model.submit()
        self.assertEqual(p,os.path.join(self.model.dirPath,'SUBMIT'))

if __name__ == "__main__":
    print "Running Test Cases"
    unittest.main() ## actually run the test cases