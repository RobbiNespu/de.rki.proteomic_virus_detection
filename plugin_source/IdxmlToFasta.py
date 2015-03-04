from re import sub
import xml.etree.ElementTree as ET
from optparse import OptionParser

op = OptionParser()

op.add_option('-i','--idxml', dest='idxmlFile', default=None)
op.add_option('-f','--fasta', dest='fastaFile', default=None)

opts, args = op.parse_args()

idxmlTree = ET.parse(opts.idxmlFile)
idxmlRoot = idxmlTree.getroot()

peptides = set()
for run in idxmlRoot.findall('IdentificationRun'):
    for peptide_id in run.findall('PeptideIdentification'):
        for hit in peptide_id.findall('PeptideHit'):
            peptide = hit.get('sequence')
            peptide = sub('\(.+?\)', '', peptide)
            peptides.add(peptide)
            
with open(opts.fastaFile, 'w+') as fastaFile:
    for peptide in peptides:
        fastaFile.write('>%s\n' % peptide)
        fastaFile.write('%s\n' % peptide)
        