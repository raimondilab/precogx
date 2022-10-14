'''
A set of test suites to check
the code using pytest
'''
#from http import client
#from urllib import response
import json, sys, os
import pytest
from flask import Flask
#sys.path.insert(0, os.getcwd())
from static.predictor import precogx

class TestClass():
    '''A class to test precogx'''
    def test_precogx(self):
        '''error code 1 = not GPCR; error code 2 = longer than 1024 AA
        In case of "P046", it must run except as that uniprot acc does not exist'''
        data = [('TP53', 1), ('P32249', 0), ('MC1R/M294L', 0), ('P046', 1), ('AGRV1_HUMAN', 2)]
        for input, eC in data:
            try:
                uniq_id, errorCode, flaggedGPCR = precogx.main(15, input, None, 'all', os.getcwd())
                if eC != errorCode:
                    pytest.fail("Error code mismatch. Expected "+str(eC)+", got "+str(erroCode)+" in "+input)
            except:
                if input != 'P046':
                    pytest.fail("PRECOGx failed to run for "+input)
