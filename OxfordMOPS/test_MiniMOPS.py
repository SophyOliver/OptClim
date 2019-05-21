"""
Test MiniMOPS model
"""

import six # allows switch between python 2 and 3
import unittest
import MiniMOPS
import os
import shutil
import collections
import pandas as pd

ROOTDIR = os.path.join(os.environ["OPTCLIMTOP"], "Configurations/OxfordMiniMOPS_Configs/NoNoise/")
TESTDIR = os.path.join(os.environ["OPTCLIMTOP"], "OxfordMOPS")

class testMINIMOPS(unittest.TestCase):
	"""
	test cases for MiniMOPS
	"""
	
	def setUp(self):
		"""
		Runs every time
		"""
		print(TESTDIR)
		print(ROOTDIR)
		parameters = {"param1":175.0, "param2":2.5} # default values for now to optimise
		testDir = os.path.join(TESTDIR, "tempDirMini")
		refDir = os.path.join(ROOTDIR, "RunCode")
		try:
			shutil.rmtree(testDir)
		except OSError:
			pass #if dir doesn't exist keep going
		
		self.dirPath = testDir
		self.refDir = refDir
		
		self.model = MiniMOPS.MiniMOPS(testDir, name="test1", create=True,
				refDirPath=refDir,
				verbose=True, parameters = parameters)
		 
	def test_init(self):
		"""
		Check testMINIMOPS initialisation
		""" 
		print('test ini =================')
		#self.assertEqual(self.model.get(["name"]), "test1")
		expectParam=collections.OrderedDict()
		for k,v in zip(["param1", "param2"], [175.0, 2.5]):
			expectParam[k] = v
		expectObs=collections.OrderedDict()
		for k in ["mis1", 'mis2', "mis3"]: 
			expectObs[k] = None
			
		self.assertEqual(self.model.get(['name']),'test1')
		self.assertEqual(self.model.get(['observations']),expectObs)
		self.assertDictEqual(self.model.get(['parameters']),expectParam)
		self.assertEqual(self.model.get(['ppOutputFile']),'misfit_output.txt')
		self.assertListEqual(self.model.get(['observations']).keys(),expectObs.keys())
		
		m= MiniMOPS.MiniMOPS(self.dirPath, verbose=True)
		# Nothing should have changed except observations have been read in
		self.assertEqual(m.get(['name']),'test1')
		self.assertDictEqual(m.get(['parameters']),expectParam)
		
		self.assertEqual(m.get(['ppOutputFile']),'misfit_output.txt')
		self.assertListEqual(m.getObs().keys(),expectObs.keys())
		self.assertEqual(m.getObs(),expectObs)

	def test_readObs(self):
		"""
		Test readObs function
		"""
		print("read obs test ===============")
		expectObs=collections.OrderedDict()
		for k,v in zip(["mis1", 'mis2', "mis3"], [2.0, 4.0, 3.0]): 
			expectObs[k] = v
		
		obsFile = os.path.join(self.model.dirPath,self.model.get(['ppOutputFile']))
		
		with open(obsFile, "w") as text_file:
			for k,v in expectObs.items():
				text_file.write("%f\n"%v)
		m= MiniMOPS.MiniMOPS(self.dirPath, verbose=True)
		self.assertListEqual(m.getObs().keys(),expectObs.keys())
		#import pdb; pdb.set_trace()
		self.assertDictEqual(m.getObs(),expectObs)
		
	def test_setParams(self):
		"""
		Test setParams function
		We want to read in and check the parameters_input.txt which setParams function should have written
		"""
		print("read params test ===============")
		expectParams=[175.0, 2.5]
		
		paramFile = os.path.join(self.model.dirPath,'parameters_input.txt')
		actualParams = [0,0]
		
		paramsdf = pd.read_csv(paramFile, header=None, index_col=False)
		print(paramsdf)
		for n,k in enumerate(expectParams):
			actualParams[n] = paramsdf.iloc[n].values[0]

		self.assertListEqual(expectParams, actualParams)
		
if __name__=="__main__":
	print("Running Test Case")
	unittest.main()
	
	
		
				
