import xml.etree.ElementTree as ET
from re import match
from optparse import OptionParser
    
def idxmlSpectraNums(idxml_file):

    idxml_tree = ET.parse(idxml_file)
    idxml_root = idxml_tree.getroot()

    peptides = []

    for run in idxml_root.findall('IdentificationRun'):
        for peptide_id in run.findall('PeptideIdentification'):
            mz = float(peptide_id.get('MZ'))
            rt = float(peptide_id.get('RT'))

            # assume that there is only one PeptideHit (nothing else detected but theoretically/by definition possible)
            charge = int(peptide_id.find('PeptideHit').get('charge'))

            peptides.append([mz, rt, charge])
    return peptides

def main(idxml_file, mzml_file, output_file, mz_precision = 0, rt_precision = 0):

    peptides = idxmlSpectraNums(idxml_file)

    mzml_tree = ET.parse(mzml_file) # load the input mzml file and parse the xml structure
    mzml_root = mzml_tree.getroot() # access point to the xml structure

    # load the namespace via regex (group 0 for first hit)
    namespace = match('\{.*\}', mzml_root.tag).group(0) # should be "http://psi.hupo.org/ms/mzml"

    # remove those spectra, that also appear (with their mz, rt and charge) in the idxml file
    num = 0
    for run in mzml_root.findall('%srun' % namespace): # there should be only one but you never know
        for spectrum_list in run.findall('%sspectrumList' % namespace): # same here
            for spectrum in spectrum_list.findall('%sspectrum' % namespace):
                # scan the spectra and retrieve the parameters rt, mz and charge
                rt = float(spectrum.find('%sscanList' % namespace).find('%sscan' % namespace).find('%scvParam' % namespace).get('value'))
                ion = spectrum.find('%sprecursorList' % namespace).find('%sprecursor' % namespace).find('%sselectedIonList' % namespace).find('%sselectedIon' % namespace)
                for cvParam in ion.findall('%scvParam' % namespace):
                    if cvParam.get('name') == 'selected ion m/z':
                        mz = float(cvParam.get('value'))
                    elif cvParam.get('name') == 'charge state':
                        charge = int(cvParam.get('value'))

                # remove if the spectrum appears in the peptide list
                for peptide in peptides:
                    if abs(mz - peptide[0]) <= mz_precision and abs(rt - peptide[1]) <= rt_precision and charge == peptide[2]:
                        spectrum_list.remove(spectrum)
                        num += 1
                        break

    print '%i spectra removed' % num
    
    # refresh the spectra ids and counts:
    for run in mzml_root.findall('%srun' % namespace): # there should be only one but you never know
        for spectrum_list in run.findall('%sspectrumList' % namespace): # same here
            index = 0
            for spectrum in spectrum_list.findall('%sspectrum' % namespace):
                spectrum.set('index', str(index))
                index += 1
            spectrum_list.set('count', str(index))
            
    #refresh software list (to include this modification step)
    software_list = mzml_root.find('%ssoftwareList' % namespace)
    # update the count
    software_count = int(software_list.get('count')) +1
    software_list.set('count', str(software_count))
    # create software subelement in the list
    software_element = ET.SubElement(software_list, 'software', attrib={'id' : 'unknown id', 'version' : '2.0.0'})
    cvParam = ET.SubElement(software_element, 'cvParam', attrib={'cvRef':"MS", 'accession':"MS:1000799", 'name':"custom unreleased software tool"})
    
    
    if num != len(peptides):
        print '%i spectra where found in the idxml file, so there seems to be something strange.\nIt maybe be of interest to adjust the precision parameters for rt and mz values' % len(peptides)
    else:
        print 'Run as expected'
    # central namespace, so it does not get added to every (sub-) element when writing the new file ( 1:-1 to exclude the brackets)
    ET.register_namespace('', namespace[1:-1])
    
    # write the new mzml file (obey the encoding as everything else leads to errors)
    mzml_tree.write(output_file, xml_declaration = True, method='xml', encoding = 'ISO-8859-1')

if __name__ == '__main__':

    op = OptionParser()

    op.add_option('-m','--mzml', dest='mzml_file', default=None, help='the input mzml file from which the peptides are removed')
    op.add_option('-p','--peptides', dest='idxml_file', default=None, help='the idXML file that contains the detected peptides')
    op.add_option('-o','--output', dest='output_file', default=None, help='the output mzml file with the remaining spectra')
    op.add_option('-r','--rt', dest='rt_precision', type="float", default=0, help='some peptide search engine seem to have floating point errors. To still backtrack peptides you can specify how much rt values can deviate to overwrite such errors')
    op.add_option('-z','--mz', dest='mz_precision', type="float", default=0, help='as with rt but for mz')

    opts, args = op.parse_args()

    main(opts.idxml_file, opts.mzml_file, opts.output_file, opts.mz_precision, opts.rt_precision)
