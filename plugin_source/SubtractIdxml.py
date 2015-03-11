import xml.etree.ElementTree as ET
from optparse import OptionParser

op = OptionParser()

op.add_option('-m','--minuend', dest='minuendFile', default=None, help='the input idxml file from which hits get removed')
op.add_option('-s','--subtrahend', dest='subtrahendFile', default=None, help='the input idxml file with hits that get gemove from the minuend files')
op.add_option('-o','--output', dest='outputFile', default=None, help='the genetrated idxml file of remaining hits from the minend file')

opts, args = op.parse_args()

subtrahendTree = ET.parse(opts.subtrahendFile)
subtrahendRoot = subtrahendTree.getroot()

subtrahendIons = []
for run in subtrahendRoot.findall('IdentificationRun'):
    for peptideIdentification in run.findall('PeptideIdentification'):
        mz = peptideIdentification.get('MZ')
        rt = peptideIdentification.get('RT')
        charge = peptideIdentification.find('PeptideHit').get('charge')
        subtrahendIons.append((mz, rt, charge))
        
print len(subtrahendIons)
        
minuendTree = ET.parse(opts.minuendFile)
minuendRoot = minuendTree.getroot()

for run in minuendRoot.findall('IdentificationRun'):
    for peptideIdentification in run.findall('PeptideIdentification'):
        mz = peptideIdentification.get('MZ')
        rt = peptideIdentification.get('RT')
        charge = peptideIdentification.find('PeptideHit').get('charge')
        
        if (mz, rt, charge) in subtrahendIons:
            run.remove(peptideIdentification)
            
minuendTree.write(opts.outputFile, xml_declaration = True, method='xml', encoding = 'UTF-8')