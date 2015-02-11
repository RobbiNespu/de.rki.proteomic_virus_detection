import xml.etree.ElementTree as ET
from ntpath import basename # determine name of a file
from optparse import OptionParser

def filename(file):
    '''
    Return just the name of a file without the path and file format
    '''
    
    file = basename(file)
    return file.split('.')[0]
    
def countPeptides(idxml_file):
    '''
    count the number of identified peptides in a idxml file and returns it (together with database name)
    '''

    idxml_tree = ET.parse(idxml_file)
    idxml_root = idxml_tree.getroot()

    # count peptides
    num_peptides = 0
    for run in idxml_root.findall('IdentificationRun'):
        for peptide_id in run.findall('PeptideIdentification'):
            num_peptides += 1
        
    # retrieve the database name
    db = idxml_root.find('SearchParameters').get('db')
    if db:
        db_name = filename(db)
    else:
        db_name = 'Unknown'
    
    return db_name, num_peptides
   
    
def main(idxml_files, output_file=None):
    
    counts = []
    
    print 'spectra name\tdb name\thit counts'
    for idxml_file in idxml_files:
        db_name, count = countPeptides(idxml_file)
        spectra_name = filename(idxml_file)
        counts.append([spectra_name, db_name, count])
        print '%s\t%s\t%i hits'%(spectra_name, db_name, count)
    
    if output_file:
        with open(output_file, 'w+') as f:
            f.write('The following numbers of peptides were detected:\n\n')
            f.write('spectra name\tdb name\thit counts\n')
            for count in counts:
                f.write('%s\t%s\t%i\n'%(count[0], count[1], count[2]))
            

if __name__ == '__main__':

    op = OptionParser()

    op.add_option('-i','--idxml', dest='idxml_files', action='append', default=None, help='input idxml files')
    op.add_option('-o','--output', dest='output_file', default=None, help='the output text file with the peptide counts')

    opts, args = op.parse_args()

    main(opts.idxml_files, opts.output_file)       
    
