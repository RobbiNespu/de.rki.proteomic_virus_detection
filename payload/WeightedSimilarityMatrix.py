from optparse import OptionParser
import numpy as np
from re import match
from Bio import SeqIO

usage = '''
%prog [options]

Creates a weighted similarity matrix for the pipasic tool
'''

def constructWeightedSimilarityMatrix(digest_files, peptide_files, database_files, output_file):
    '''
    Objective:
    
    construct a weigthed similarity matrix
    (return: proteome similarity matrix weighted by dataset in output_file)
    '''
    
    # make sure we have the same amount of datasets from digest and identified peptides
    assert len(digest_files) == len(peptide_files)
    
    similarity_matrix = np.zeros([len(digest_files),len(digest_files)]) # proteome similarity matrix
    
    for i, files in enumerate(zip(peptide_files, digest_files)):
        peptide_file, digest_file = files
        
        # peptides are already found => identifiedPeptides
        # and the databases are already digested
        weighted_peptide_file = weightPeptides(peptide_file, digest_file)
        
        for j, database_file in enumerate(database_files):
            pMatch, cMatch, cPept = weightedMatching(weighted_peptide_file, database_file) # may need more arguments
            similarity_matrix[i, j] = pMatch
    
    # write the similarity matrix into a file (human unreadable format)
    np.savetxt(output_file, similarity_matrix, delimiter=',\t')
            
def weightPeptides(identifiedPeptides, trypticPeptides, outfile='weighthedPeptides.txt', init=0, verbose=True):
    """
    search identified peptides in tryptic peptide list of a proteome in order 
    to weight each peptide
    (return: path to output file with line structure: "weight    sequence")
    
    identifiedPeptides: inspectparser.py output (experimental spectra)
    trypticPeptides:    typticCut result of proteome in question
    init:               initial weight for all peptides
    """
    
    ############### Read the digested Peptides ###############
    
    trypList = [] # list of tryptic peptides
    trypDic = dict() # dictionary of all peptide sequences and corresponding indices in trypList
    protDic = dict() # dictionary of protein names and corresponding indices in trypList
    
    if verbose: print "reading sequences of tryptic peptides ..."
    
    for idx, pept in enumerate(SeqIO.parse(open(trypticPeptides, "rU"), "fasta")) :
        
        # dictionary of protein names (descriptions)
        p_name = str(pept.description)
        if p_name in protDic: # add peptide (by index) to its protein
            protDic[p_name].append(idx)
        else: # add new protein name (+ first entry) to dictionary
            protDic[p_name] = [idx]
            
        # peptide sequence
        pept = str(pept.seq)
        trypList.append(pept)
        
        # dictionary of peptide sequences indicating unique and shared peptides
        if pept in trypDic: # add index of current peptide to already known sequence
            trypDic[pept].append(idx)
        else: # add new entry to sequence dictionary
            trypDic[pept] = list([idx])

    if verbose: print len(trypDic), "different sequences out of", len(trypList), "peptides from", len(protDic), "proteins"
    
    # array for peptide match counts
    counts = [init]*len(trypList)
    
    # counter for peptides not in tryptic peptide list
    notfound = 0
    
    ############### Comparison to identified peptides ###############
    
    if verbose: print "reading and searching sequences of identified peptides ..."
    
    with open(identifiedPeptides,"r") as inF:
        
        #file processing to extract the sequences
        text = inF.read()
        peptide_hits = text.split('<PeptideIdentification ')[1:] # => '... sequence="KJAHGSNNK" ....'
        peptides = [peptide_hit.split('sequence="')[1].split('"')[0] for peptide_hit in peptide_hits] # => 'KJAHGSNNK'
        
        for pept in peptides:
            ### split identified peptide into tryptic peptides (if neccessary)
            
            # cut after lysine:
            pepsLys = [pep+"K" for pep in pept.split("K")]
            pepsLys[-1] = pepsLys[-1][:-1]
            # remove cut before proline:
            pepsLysP = [pepsLys[0]]
            for pe in range(1,len(pepsLys)):
                if match("P",pepsLys[pe]): # is this correct, don't we only look for starting proline
                    pepsLysP[-1] += pepsLys[pe]
                else: pepsLysP.append(pepsLys[pe])
            # cut after arginine: 
            for pept in pepsLysP:
                pepsArg = [pep+"R" for pep in pept.split("R")]
                pepsArg[-1] = pepsArg[-1][:-1]
                # remove cut before proline:
                pepsArgP = [pepsArg[0]]
                for pe in range(1,len(pepsArg)):
                    if match("P",pepsArg[pe]):
                        pepsArgP[-1] += pepsArg[pe]
                    else: pepsArgP.append(pepsArg[pe])
                for pep in pepsArgP:
                    # only accept peptides of length 7 - 25
                    if len(pep) > 6 and len(pep) < 26: 
                        # search each peptide in list of proteome tryptic peptides
                        if pep in trypList:
                            counts[trypList.index(pep)] += 1
                        else:
                            notfound += 1
                                    
    # calculate weights from counts
    if verbose: print "not found: %i of %i"%(notfound,sum(counts)+notfound)
    if sum(counts) > 0:
        if verbose: print "recalculating counts for shared peptides ..."
        for t in trypList: # equally distribute counts for peptides of same sequence
            t_ind = trypDic.get(t)
            x = sum([counts[c] for c in t_ind])/float(len(t_ind))
            for i in t_ind:
                counts[i] = x
    
        sum_c = sum(counts)
        for p in protDic: # equally distribute counts for peptides of same protein
            # normalized by the sum of all counts (sum_c)
            p_ind = protDic.get(p)
            x = sum([counts[c] for c in p_ind])/float(len(p_ind))/sum_c
            for i in p_ind:
                counts[i] = x
    else:
        print "Warning: no peptides identified!"
    if verbose: print "writing resulting weights ..."
    cPos = 0 # number of positively weighted peptides
    outF = open(outfile,"w+")
    # output file line structure: "weight    sequence"
    for idx, seq in enumerate(trypList):
        outF.write(str(counts[idx])+"\t"+seq+"\n")
        if not counts[idx] == 0:
            cPos += 1
    outF.close()
    if verbose: print cPos,"of",len(trypList),"peptides got positive weight."
    return outfile
    
def weightedMatching(weightedPeptFile, dbFile, sum_all=True, verbose=True):
    """
    (return: weighted sum of matches, number of matches, number of peptides)
    search all weighted peptide sequences from a given file in a given database
    
    weightedPeptFile: textfile of weighted peptide sequences 
             (weightPeptides output file line structure: "count    sequence")
    dbFile:  protein database (fasta-format)
    sum_all: add up weights of matching peptides; if false, set count of 
             matching peptides with weight > 0 to 1 (only suitable for init=0)
    """
    
    DB = [] #list of protein sequences from dbFile
    cPept = 0 # number of peptides
    cMatch = 0 # number of matches (weighted peptide <-> DB peptide)
    pMatch = float(0) # weighted sum of matches
    
    if verbose: print "reading DB sequences ..."
    for record in SeqIO.parse(open(dbFile, "rU"), "fasta") :
        DB.append(str(record.seq).upper())
    if not sum_all: prot_count = [0]*len(DB)
    if verbose: print "reading and mapping peptide sequences ..."
    
    for line in open(weightedPeptFile, "rU"):
        if float(line.split("\t")[0]) > 0: #this line is NEW (irritates cmatches and cPept)
            cPept += 1
            peptide_seq = line.split("\t")[1].rstrip().upper() # another minor (speed) change
            for i,DBi in enumerate(DB): # search process (weighted peptide in DB)
                if peptide_seq in DBi:
                    cMatch += 1
                    if sum_all: pMatch += float(line.split("\t")[0])
                    elif (float(line.split("\t")[0]) > 0): prot_count[i] = 1
                    break
                
    if not sum_all: #normalization (neccessary for count data)
        pMatch = float(sum(prot_count))/len(DB)
        
    if verbose: print "\r%i of %i peptides found" %(cMatch,cPept)
    return pMatch, cMatch, cPept # weighted sum of matches, #matches, #peptides
    
if __name__ == '__main__':

    op = OptionParser(usage=usage)
    
    op.add_option('-d','--digest', dest='digest_files', action='append', default=None, help='files with the (tryptic) digest of the .fasta protein libaries')
    op.add_option('-p','--peptides', dest='peptide_files', action='append', default=None, help='(.idXML) files of the found peptides for the different species')
    op.add_option('-b','--db', dest='database_files', action='append', default=None, help='(.fasta) database files for the different species')
    op.add_option('-o','--output', dest='output_file', default = None, help='the file the similarity matrix is output to (currently a python pickle file)')
    
    op.add_option('-i','--init', default=0, help = 'initial weight for all peptides')
    
    opts, args = op.parse_args()
    
    constructWeightedSimilarityMatrix(opts.digest_files, opts.peptide_files, opts.database_files, opts.output_file)