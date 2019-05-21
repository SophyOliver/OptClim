"""
Place to put tests for Submit. Currently empty.. 
"""
import unittest
from OptClimVn2 import Submit
import collections
import os
import shutil

# this routine should be available in a library! 
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



class testSubmit(unittest.TestCase):
    """
    Test cases for Submit There should be one for every method in Submit.

    """
    def setUp(self):
        """
        Standard setup for all test cases
        :return:
        """
        testDir='$OPTCLIMTOP/test_out'
        refDir='$OPTCLIMTOP/test_in'

        self.dirPath=testDir
        self.refPath=refDir

        testDir=os.path.expanduser(os.path.expandvars(testDir))
        refDir=os.path.expandvars(os.path.expanduser(refDir))
        shutil.rmtree(testDir,onerror=errorRemoveReadonly)
        
        # now do stuff. 


