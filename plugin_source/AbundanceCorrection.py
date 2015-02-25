from optparse import OptionParser
import numpy as np
from PASiC import PASiC

def determineNumberOfSpectra(spectra_file):
    with open(spectra_file, 'r') as f:
        for line in f.readlines(): #generator, does not load the whole file as we only need to access a certain line and the file might be huge
            if 'spectrumList' in line: # does probably not work very well with files that contain multiple spectrumLists
                spectra_count = int(line.split('spectrumList count="')[1].split('"')[0]) # '... <spectrumList count="10001" d ...'
                break
    return spectra_count

def countPeptides(peptide_file):
    with open(peptide_file, 'r') as f:
        peptide_count = len(f.read().split('<PeptideIdentification ')) -1
    return peptide_count

def calculateAbundance(spectra_file, peptide_files, matrix_file, relative_abundance_file = None, corrected_counts_file = None):
    # Count the number of peptide hits for the various species
    peptide_counts = []
    for peptide_file in peptide_files:
        peptide_counts.append(countPeptides(peptide_file))
    peptide_counts = np.array(peptide_counts)

    print peptide_counts

    # count/determine the number of spectra = sample size
    spectra_count = determineNumberOfSpectra(spectra_file)

    print spectra_count

    # load the similarity matrix
    similarity_matrix = np.loadtxt(matrix_file, delimiter=',\t')

    # run the (Pi-)Pasic correction
    wtd_Abunds, relative_abundancies = PASiC(simiMat=similarity_matrix, counts=peptide_counts, N=spectra_count)
    corrected_counts = wtd_Abunds * spectra_count

    # output the specified files
    if relative_abundance_file:
        np.savetxt(relative_abundance_file, relative_abundancies, delimiter=',\t')
    if corrected_counts_file:
        np.savetxt(corrected_counts_file, corrected_counts, delimiter=',\t')

if __name__ == '__main__':

    op = OptionParser()

    op.add_option('-s','--spectra', dest='spectra_files', default=None, help='the original spectra file(.mzml) to count the number of spectras')
    op.add_option('-p','--peptides', dest='peptide_files', action='append', default=None, help='(.idXML) files of the found peptides for the different species to calculate the numbers')
    op.add_option('-m','--matrix', dest='matrix_file', default=None, help='the calculated matrix (pickle object)')
    op.add_option('-r','--relativeAbundance', dest='relative_abundance_file', default = None, help='returns file with relative abundancies')
    op.add_option('-c','--correctedCounts', dest='corrected_counts_file', default = None, help='returns file with relative abundancies')

    opts, args = op.parse_args()

    calculateAbundance(opts.spectra_files, opts.peptide_files, opts.matrix_file, opts.relative_abundance_file, opts.corrected_counts_file)
