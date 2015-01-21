# is run: python setup_py2exe.py py2exe
from distutils.core import setup
import py2exe

setup(
    console=['WeightedSimilarityMatrix.py', 'AbundanceCorrection.py'],
    options={
        'py2exe':{
            'includes' : ['scipy.sparse.csgraph._validation']}}) # necessary to prevent error: http://stackoverflow.com/questions/14215303/scipy-with-py2exe