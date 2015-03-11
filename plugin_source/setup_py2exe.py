# is run: python setup_py2exe.py
from distutils.core import setup
import py2exe, sys

sys.argv.append('py2exe') # so we do not have to add py2exe to the arguments

setup(
    console=[
        'WeightedSimilarityMatrix.py',
        'AbundanceCorrection.py',
        'ExcludePeptides.py',
        'CountPeptides.py',
        'ConcatenateFiles.py',
        'IdxmlToFasta.py',
        'SixFrameTranslation.py',
        'SubtractIdxml.py'],
    options={
        'py2exe':{
            'includes' : ['scipy.sparse.csgraph._validation']}}) # necessary to prevent error: http://stackoverflow.com/questions/14215303/scipy-with-py2exe