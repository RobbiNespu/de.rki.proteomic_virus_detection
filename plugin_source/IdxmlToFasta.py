import xml.etree.ElementTree as ET
from optparse import OptionParser
from Bio import SeqIO

op = OptionParser()

op.add_option('-i','--idxml', dest='idxmlFile', default=None, help='the input idxml file with peptide hits and associated proteins')
op.add_option('-p','--proteome', dest='proteomeFile', default=None, help='the input fasta proteome file with the proteins whose sequeces are retrieved')
op.add_option('-f','--fasta', dest='fastaFile', default=None, help='the output fasta file with the found proteins')

opts, args = op.parse_args()

idxmlTree = ET.parse(opts.idxmlFile)
idxmlRoot = idxmlTree.getroot()

# load the found proteins from the idxml file
idxmlProteins = []
for run in idxmlRoot.findall('IdentificationRun'):
    for proteinIdentification in run.findall('ProteinIdentification'):
        for proteinHit in proteinIdentification.findall('ProteinHit'):
            idxmlProteins.append(proteinHit.get('accession'))

foundProteins = []
for protein in SeqIO.parse(open(opts.proteomeFile), 'fasta'):
    if protein.id in idxmlProteins:
        foundProteins.append(protein)
        
SeqIO.write(foundProteins, open(opts.fastaFile, 'w+'), 'fasta')
        
if not len(idxmlProteins) == len(foundProteins):
    print 'Only %i of %i proteins from the idxml could be found in the genomic database' % (len(foundProteins), len(idxmlProteins))
    print 'Make sure to use the right genomic database'
    pids = [p.id for p in foundProteins]
    for ip in idxmlProteins:
        if not ip in pids:
            print ip
else:
    print 'The script run with out problems, %i proteins were identified and copied to the new fasta file' % len(idxmlProteins)

        