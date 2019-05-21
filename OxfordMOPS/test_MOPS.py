"""
Test MOPS model
"""

import six # allows switch between python 2 and 3
import unittest
import MOPS
import os
import shutil
import collections
import pandas as pd

ROOTDIR = os.path.join(os.environ["OPTCLIMTOP"], "Configurations/OxfordMOPS_Configs/")
TESTDIR = os.path.join(os.environ["OPTCLIMTOP"], "OxfordMOPS")

class testMOPS(unittest.TestCase):
	"""
	test cases for MOPS
	"""
	
	def setUp(self):
		"""
		Runs every time
		"""
		print(TESTDIR)
		print(ROOTDIR)
		parameters = {"ro2ut":170.0, "ACik":24.0}; # default values for now to optimise
		testDir = os.path.join(TESTDIR, "tempDir")
		refDir = os.path.join(ROOTDIR, "RunCode")
		try:
			shutil.rmtree(testDir)
		except OSError:
			pass #if dir doesn't exist keep going
		
		self.dirPath = testDir
		self.refDir = refDir
		
		self.model = MOPS.MOPS(testDir, name="test1", create=True,
				refDirPath=refDir,
				verbose=True, parameters = parameters)
		 
	def test_init(self):
		"""
		Check testMOPS initialisation
		""" 
		print('test ini =================')
		#self.assertEqual(self.model.get(["name"]), "test1")
		expectParam=collections.OrderedDict()
		for k,v in zip(["ro2ut", "ACik"], [170.0, 24.0]):
			expectParam[k] = v
		expectObs=collections.OrderedDict()
		for k in ["mis1_PO4", "mis2_PO4", "mis3_PO4", "mis4_PO4", "mis5_PO4", "mis6_PO4", "mis7_PO4", "mis8_PO4", "mis9_PO4", "mis10_PO4", "mis11_PO4", "mis12_PO4", "mis13_PO4", "mis14_PO4",
						"mis1_NO3", "mis2_NO3", "mis3_NO3", "mis4_NO3", "mis5_NO3", "mis6_NO3", "mis7_NO3", "mis8_NO3", "mis9_NO3", "mis10_NO3", "mis11_NO3", "mis12_NO3", "mis13_NO3", "mis14_NO3",
						"mis1_O2",  "mis2_O2",  "mis3_O2",  "mis4_O2",  "mis5_O2",  "mis6_O2",  "mis7_O2",  "mis8_O2",  "mis9_O2",  "mis10_O2",  "mis11_O2",  "mis12_O2",  "mis13_O2",  "mis14_O2"]: 
			expectObs[k] = None
			
		self.assertEqual(self.model.get(['name']),'test1')
		self.assertEqual(self.model.get(['observations']),expectObs)
		self.assertDictEqual(self.model.get(['parameters']),expectParam)
		self.assertEqual(self.model.get(['ppOutputFile']),'misfit_output.txt')
		self.assertListEqual(self.model.get(['observations']).keys(),expectObs.keys())
		
		m= MOPS.MOPS(self.dirPath, verbose=True)
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
		for k,v in zip(["mis1_PO4", "mis2_PO4", "mis3_PO4", "mis4_PO4", "mis5_PO4", "mis6_PO4", "mis7_PO4", "mis8_PO4", "mis9_PO4", "mis10_PO4", "mis11_PO4", "mis12_PO4", "mis13_PO4", "mis14_PO4",
						"mis1_NO3", "mis2_NO3", "mis3_NO3", "mis4_NO3", "mis5_NO3", "mis6_NO3", "mis7_NO3", "mis8_NO3", "mis9_NO3", "mis10_NO3", "mis11_NO3", "mis12_NO3", "mis13_NO3", "mis14_NO3",
						"mis1_O2",  "mis2_O2",  "mis3_O2",  "mis4_O2",  "mis5_O2",  "mis6_O2",  "mis7_O2",  "mis8_O2",  "mis9_O2",  "mis10_O2",  "mis11_O2",  "mis12_O2",  "mis13_O2",  "mis14_O2"],
						[2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0,
						 2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0,
						 2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0, 3.0, 2.0, 4.0]): 
			expectObs[k] = v
		
		obsFile = os.path.join(self.model.dirPath,self.model.get(['ppOutputFile']))
		
		with open(obsFile, "w") as text_file:
			for k,v in expectObs.items():
				text_file.write("%f\n"%v)
		m= MOPS.MOPS(self.dirPath, verbose=True)
		self.assertListEqual(m.getObs().keys(),expectObs.keys())
		#import pdb; pdb.set_trace()
		self.assertDictEqual(m.getObs(),expectObs)
		
	def test_setParams(self):
		"""
		Test setParams function
		We want to read in and check the parameters_input.txt which setParams function should have written
		"""
		print("read params test ===============")
		expectParams=[170.0, 24.0]
		expectNames = ["ro2ut", "ACik"]
		expectNum = 2
		
		paramFile = os.path.join(self.model.dirPath,'parameters_input.txt')
		nameFile = os.path.join(self.model.dirPath,'parameters_names.txt')
		numFile = os.path.join(self.model.dirPath,'num_bgc_params.txt')
		actualParams = [0,0]
		actualNames = ["",""]
		actualNum = 0;
		
		paramsdf = pd.read_csv(paramFile, header=None, index_col=False)
		namesdf = pd.read_csv(nameFile, header=None, index_col=False)
		numdf = pd.read_csv(numFile, header=None, index_col=False)
		print(paramsdf)
		print(namesdf)
		print(numdf)
		for n,k in enumerate(expectParams):
			actualParams[n] = paramsdf.iloc[n].values[0]
			actualNames[n] = namesdf.iloc[n].values[0]
		actualNum = numdf.iloc[0].values[0]

		self.assertListEqual(expectParams, actualParams)
		self.assertListEqual(expectNames, actualNames)
		self.assertEqual(expectNum, actualNum)
		
if __name__=="__main__":
	print("Running Test Case")
	unittest.main()
	
	
		
				
