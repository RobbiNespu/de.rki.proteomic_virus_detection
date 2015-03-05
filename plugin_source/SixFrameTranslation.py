from Bio import SeqIO
from Bio.Alphabet import IUPAC, ProteinAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser

def translate_to_six_frame(dnaRecord, translationTable):
    '''
    translate a Bio.SeqRecord of a DNA sequence via the given translation table into the six possible translations
    
    dnaRecord = Bio.SeqRecord of DNA sequence (or other)
    translationTable = the codon table for translating base triplets into amino acids (number between 1 and 25 based on http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes)
    '''
    
    translations = []
    
    for frame in range(3):
        for direction in ['forward', 'reverse']:
            if direction == 'forward':
                sequence = dnaRecord.seq[frame:]
            else:
                sequence = dnaRecord.seq.reverse_complement()[frame:]
            aaSeq = Seq(str(sequence.translate(translationTable)), alphabet = ProteinAlphabet)

            aaRecord = SeqRecord(aaSeq, dnaRecord.name)
            aaRecord.id = '%s_%s%i' % (aaRecord.id, direction[0], frame)
            aaRecord.description = '%s|translation %s frame %i' % (aaRecord.description, direction, frame)
            
            translations.append(aaRecord)
            
    return translations
    
if __name__ == '__main__':

    op = OptionParser()

    op.add_option('-g','--genomes', dest='genomeFilenames', action='append', default=None, help='the input genome in fasta format')
    op.add_option('-o','--outputs', dest='outputFilenames', action='append', default=None, help='the output fasta file with the six frame translation')
    op.add_option('-t','--translTable', dest='translationTableNumber', default=1, type='int', help='a translation table number according to http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes')

    opts, args = op.parse_args()
    
    for genomeFilename, outputFilename in zip(opts.genomeFilenames, opts.outputFilenames):
    
        translations = []
        for dnaRecord in SeqIO.parse(open(genomeFilename), 'fasta'):
            translations.extend(translate_to_six_frame(dnaRecord, opts.translationTableNumber))
            
        SeqIO.write(translations, open(outputFilename, 'w+'), 'fasta')

    
