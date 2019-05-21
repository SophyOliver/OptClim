"""
Test code for StudyConfig
"""

import unittest
import StudyConfig
import os
import pandas as pd
import collections
import numpy as np
import xarray

class testStudyConfig(unittest.TestCase):
    """
    Test cases for StudyConfig. Currently partial so as timee passes more tests should be written. 
    Ideally >1 for each method
    """
    def setUp(self):
        """
        Standard setup for all test cases
        :return:
        """


        configFile=os.path.join('Configurations','example.json')
        self.config=StudyConfig.readConfig(configFile, ordered=True)

    def test_version(self):
        """
        Test version as expected.
        :return: 
        """
        version=self.config.version()
        self.assertEqual(version,2,msg='Expected version = 2 got %s'%version)

    def test_begin(self):
        """
        Test that begin returns expected values
        :return: 
        """

        expect=pd.Series({"CT": 1.1e-4,  "EACF": 0.51,  "ENTCOEF": 3.1,  "ICE_SIZE": 31e-6, "RHCRIT": 0.71,
             "VF1":1.0 ,"CW_LAND": 2.1e-4}) # values from example.json.
        # first test the parameters are as we expect:
        params=self.config.paramNames()
        params.sort()
        expectKeys=list(expect.keys())
        expectKeys.sort()
        self.assertListEqual(params,expectKeys,msg='Paramters differ')
        got=self.config.beginParam()
        got.sort_values(inplace=True)
        expect.sort_values(inplace=True)
        self.assertEqual(expect.equals(got),True,msg='Differences in beginParam')

    def test_version(self):
        """
        Test the version is 2 as expeved
        :return: 
        """

        vn=self.config.version()
        self.assertEqual(vn,2,msg='Got vn = %d expected %d'%(vn,2))

    def test_referenceConfig(self):
        """
        Test that the ref config path works
        :return: 
        """

        ref=self.config.referenceConfig()
        #expect=os.path.expandvars("$OPTCLIMTOP/Configurations/HadAM3_ed3_SL7_15m")
        expect=os.path.join("Configurations","HadAM3_ed3_SL7_15m")
        # example json set up for unix (so convert any / to \
        ref=ref.replace('/','\\')
        self.assertEqual(ref,expect,msg='expected %s got %s'%(expect,ref))

    def test_cacheFile(self):
        """
        Test that cacheFile method works
        :return: 
        """

        expect='cache_file.json'
        got=self.config.cacheFile()
        self.assertEqual(got,expect,msg='cache_file different got %s expected %s'%(got,expect))

    def test_ranges(self):
        """
        Test ranges as expected
        :return: 
        """
        #pNames=self.config.paramNames()
        minP=pd.Series({'CT':5e-5,'EACF':0.5,'ENTCOEF':0.6,'ICE_SIZE':2.5e-5,'RHCRIT':0.6,'VF1':0.5,'CW_LAND':1e-4},
                       name='minParam')
        maxP = pd.Series({'CT': 4e-4, 'EACF': 0.7, 'ENTCOEF': 9.0, 'ICE_SIZE': 4e-5, 'RHCRIT': 0.9, 'VF1': 2.0, 'CW_LAND': 2e-3},
            name='maxParam')
        rng=(maxP-minP).rename('rangeParam',inplace=True)
        got=self.config.paramRanges(minP.index.values)
        self.assertEqual(got.loc['minParam',:].equals(minP),True,msg='Min Values differ')
        self.assertEqual(got.loc['maxParam', :].equals(maxP), True, msg='Max Values differ')
        self.assertEqual(got.loc['rangeParam', :].equals(rng), True, msg='Ranges differ')


    def test_steps(self):
        """
        Test steps are as expected
        :return: 
        """

        step=pd.Series({'CT':1e-5,'EACF':0.02,'ENTCOEF':0.15,'ICE_SIZE':1.5e-6,'RHCRIT':0.01,'VF1':0.1,'CW_LAND':2e-4},
                       name='steps')
        got=self.config.steps(step.index.values)
        self.assertEqual(got.equals(step), True, msg='Step values differ')
        # test that 10% part works -- note that test is specific to vn 2.
        self.config.Config['Parameters']['steps']['VF1']=None
        rng = self.config.paramRanges(step.index.values)
        step['VF1']=0.1*rng.loc['rangeParam','VF1']
        got = self.config.steps(step.index.values)
        self.assertEqual(got.equals(step), True, msg='Step values differ')

    def test_targets(self):
        """
        Test targets are as expected
        :return: 
        """
        # values ebl;ow are cut-n-paste from input json file!
        tgt=pd.Series({"rsr_nhx": 102.276779013,
                      "rsr_tropics": 94.172585284,
                      "rsr_shx": 108.113226173,
                      "rsr_gm": 99.5,
                      "olr_nhx": 222.981135462,
                      "olr_tropics": 259.889979681,
                      "olr_shx": 216.123659078,
                      "olr_gm": 239.6,
                      "lat_nhx": 275.837176497,
                      "lat_tropics": 297.551167045,
                      "lat_shx": 287.433252179,
                      "lprecip_nhx": 1.67274541294e-05,
                      "lprecip_tropics": 3.61223235446e-05,
                      "lprecip_shx": 2.23188435704e-05,
                      "mslp_nhx_dgm": 3.30683773931e+02,
                      "mslp_tropics_dgm": 1.78755937185e+02,
                      "temp@500_nhx": 251.407284938,
                      "temp@500_tropics": 266.733035653,
                      "temp@500_shx": 248.927897989,
                      "rh@500_nhx": 53.4221821948,
                      "rh@500_tropics": 33.9426666031,
                      "rh@500_shx": 52.6728738156})
        got=self.config.targets(obsNames=tgt.index.values)
        self.assertEqual(got.equals(tgt), True, msg='Target values differ')
        # test the constraint works too.
        got=self.config.constraintTarget()
        self.assertEqual(got.values, 0.5, msg='Constraint values differ')

    def test_Fixed(self):
        """
        Test can read fixed parameters
        :return: 
        """

        fix=self.config.fixedParams()
        expect=collections.OrderedDict([ (u'START_TIME', [1998,12,1]), (u'RUN_TARGET', [6, 3, 0])])
        self.assertEqual(expect,fix,msg='fix not as expected')

    def test_runTime(self):
        """
        Test that runtime is None
        :return: 
        """
        self.assertIsNone(self.config.getv("runTime"))

    def test_optimumParam(self):
        """
        Test that optimum param works.
        :return:
        """
        # nothing there to start with so should get None back.

        optParam=self.config.optimumParams()
        self.assertIsNone(optParam,msg='optParam should  be None')
        # now set it
        # need to set parameter list.
        res=self.config.optimumParams(RHCRIT=2.1,VF1=2.005,ENTCOEF=3.2)
        expect=pd.Series({'VF1':2.005,'ENTCOEF':3.2,'RHCRIT':2.1,'CT':1e-4,'EACF':0.5,'ICE_SIZE':30e-6,'CW_LAND':2e-4})
        expect=expect[self.config.paramNames()]
        self.assertTrue(res.equals(expect))

    def test_obsNames(self):
        """
        Test obsnames works
        :return:
        """

        obs=self.config.obsNames()
        expect=self.config.getv('study')['ObsList']
        self.assertEqual(obs,expect)

        # now with constrain wrapped in.
        obs = self.config.obsNames(add_constraint=True)
        expect.append(self.config.constraintName())
        self.assertEqual(obs, expect)


    def test_GNsetget(self):
        """
        test that GNsetget works as expects
        :return:
        """
        import numpy as np
        # first case -- nothing set then should get None back and self.getv('GNinfo') should return None
        self.config.GNgetset('slartybartfast')
        self.assertIsNone(self.config.getv('GNinfo'))
        value=np.arange(1.0,20.0)
        got=self.config.GNgetset('slartybartfast',value)
        self.assertTrue(np.array_equal(value,got))
        # different value works
        value=value*2
        got=self.config.GNgetset('slartybartfast',value)
        self.assertTrue(np.array_equal(value, got))
        # no having set it should just get it back.
        got=self.config.GNgetset('slartybartfast')
        self.assertTrue(np.array_equal(value, got))

        #self.fail()
    def test_GNjacobian(self):
        """
        test that GNjcobian works.
        :return:
        """

        # nothing to start with -- should get None
        self.assertIsNone(self.config.GNjacobian(), msg='Jacobian should  be None')
        # now check we get an xarray back with expected values.
        paramNames=self.config.paramNames()
        obsNames=self.config.obsNames(add_constraint=True)
        jac=np.random.uniform(0,1,(2,len(paramNames),len(obsNames))) # random numpy array
        expect=xarray.DataArray(jac,coords={'Iteration': np.arange(0,2), 'Parameter': paramNames, 'Observation': obsNames},
                                    dims=['Iteration', 'Parameter', 'Observation'], name=self.config.name())
        jacx=self.config.GNjacobian(jac,constraint=True)
        self.assertTrue(jacx.equals(expect),msg='failed in GNjacobian')
        # and now it is set should get back same results.
        jacx=self.config.GNjacobian(constraint=True)
        self.assertTrue(jacx.equals(expect),msg='failed in GNjacobian#2')


    def test_GNparams(self):
        """
        Test that GNparams works as expected.
        :return:
        """

        # nothing to start with -- should get None
        self.assertIsNone(self.config.GNparams(), msg='Params should  be None')
        # now check we get an xarray back with expected values.
        paramNames = self.config.paramNames()
        par = np.random.uniform(0, 1, (2, len(paramNames)))  # random numpy array
        expect = xarray.DataArray(par, coords={'Iteration': np.arange(0, 2), 'Parameter': paramNames},
                                  dims=['Iteration', 'Parameter'], name=self.config.name())
        parx = self.config.GNparams(par)
        self.assertTrue(parx.equals(expect), msg='failed in GNparams')
        # and now it is set should get back same results.
        parx = self.config.GNparams()
        self.assertTrue(parx.equals(expect), msg='failed in GNparams#2')

    def test_GNhessian(self):
        """
        Test that GNhessian works as expected
        :return:
        """

        # nothing to start with -- should get None
        self.assertIsNone(self.config.GNhessian(), msg='Hessian should  be None')
        # now check we get an xarray back with expected values.
        paramNames = self.config.paramNames()
        hes = np.random.uniform(0, 1, (2, len(paramNames), len(paramNames)))  # random numpy array
        expect = xarray.DataArray(hes, coords={'Iteration': np.arange(0, 2), 'Parameter': paramNames,'Parameter_2':paramNames},
                                  dims=['Iteration', 'Parameter','Parameter_2'], name=self.config.name())
        hesx = self.config.GNhessian(hes)
        self.assertTrue(hesx.equals(expect), msg='failed in GNhessian')
        # and now it is set should get back same results.
        hesx = self.config.GNhessian()
        self.assertTrue(hesx.equals(expect), msg='failed in GNhessian#2')

    def test_GNcost(self):
        """
        Test that GNcost works as expected
        :return:
        """
        # nothing to start with -- should get None
        self.assertIsNone(self.config.GNcost(), msg='Cost should  be None')
        # now check we get an xarray back with expected values.
        cost = np.random.uniform(0, 1, (2))  # random numpy array
        expect = pd.Series(cost,index=np.arange(0,2),  name=self.config.name())
        expect.index.rename('Iteration',inplace=True)
        costx = self.config.GNcost(cost)
        self.assertTrue(costx.equals(expect), msg='failed in GNcost')
        # and now it is set should get back same results.
        costx = self.config.GNcost()
        self.assertTrue(costx.equals(expect), msg='failed in GNcost#2')


    def test_GNalpha(self):
        """
        Test that GNalpha works as expected
        :return:
        """
        # nothing to start with -- should get None
        self.assertIsNone(self.config.GNalpha(), msg='alpha should  be None')
        # now check we get an pd.Series back with expected values.
        alpha = np.random.uniform(0, 1, (2))  # random numpy array
        expect = pd.Series(alpha,index=np.arange(0,2),  name=self.config.name())
        expect.index.rename('Iteration',inplace=True)
        alphax = self.config.GNalpha(alpha)
        self.assertTrue(alphax.equals(expect), msg='failed in GNalpha')
        # and now it is set should get back same results.
        alphax = self.config.GNalpha()
        self.assertTrue(alphax.equals(expect), msg='failed in GNalpha#2')

if __name__ == "__main__":
    print "Running Test Cases"
    unittest.main() ## actually run the test cases