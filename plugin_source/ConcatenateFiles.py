'''
This tool is intended to concatenate two (or more) files into one file.

Purpose: connect several proteoms to one
'''
from optparse import OptionParser
import fileinput

if __name__ == '__main__':

    op = OptionParser()

    op.add_option('-i','--input', dest='input_files', action='append', default=None, help='input fasta files')
    op.add_option('-o','--output', dest='output_file', default=None, help='the output text file with the peptide counts')

    opts, args = op.parse_args()
    
    with open(opts.output_file, 'w') as fout:
        for line in fileinput.input(opts.input_files):
            fout.write(line)