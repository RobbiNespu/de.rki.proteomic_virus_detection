from Bio import SeqIO
from Bio.Alphabet import IUPAC, ProteinAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser

def translate_to_six_frame(dnaRecord, translationTable):
    
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

    op.add_option('-g','--genome', dest='genomeFilename', default=None, help='the input genome in fasta format')
    op.add_option('-o','--output', dest='outputFilename', default=None, help='the output fasta file with the six frame translation')
    op.add_option('-t','--translTable', dest='translationTableNumber', default=1, type='int', help='a translation table number according to http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes')

    opts, args = op.parse_args()
    
    translations = []
    for dnaRecord in SeqIO.parse(open(opts.genomeFilename), 'fasta'):
        translations.extend(translate_to_six_frame(dnaRecord, opts.translationTableNumber))
        
    SeqIO.write(translations, open(opts.outputFilename, 'w+'), 'fasta')

    
