""" msgfdb2pepxml - converter from MS-GFDB output to pepXML

Parsing of command-line argument requires Python 2.7 or additional
libraries. The converter does not work with default installations
of Python 2.6.x or 3.2.1.

COPYRIGHT NOTICE

Copyright 2011 Boris Nagaev, Ksenia Yashina and Magnus Palmblad

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

## \mainpage MS-GFDB to pepXML converter
#
# \ref msgfdb2pepxml - documentation on converter
#

import argparse
import re
import os.path
from cStringIO import StringIO
from decimal import Decimal as D
from xml.etree.cElementTree import ElementTree, Element, SubElement, parse
from datetime import datetime
from warnings import warn

aa2mass = {
'G': D('57.02147'),
'A': D('71.03712'),
'S': D('87.03203'),
'P': D('97.05277'),
'V': D('99.06842'),
'T': D('101.04768'),
'C': D('103.00919'),
'I': D('113.08407'),
'L': D('113.08407'),
'N': D('114.04293'),
'D': D('115.02695'),
'Q': D('128.05858'),
'K': D('128.09497'),
'E': D('129.04260'),
'M': D('131.04049'),
'H': D('137.05891'),
'F': D('147.06842'),
'R': D('156.10112'),
'Y': D('163.06333'),
'W': D('186.07932'),
}

atom2mass = {
'H': D('1.0078250'),
'C': D('12'),
'N': D('14.0030740'),
'O': D('15.9949146'),
'P': D('30.9737615'),
'S': D('31.9720707'),
}

enzyme2name = {
'ArgC': 'Arg-C|R|P|C',
'LysC': 'Lys-C|K|P|C',
'Trypsin': 'Trypsin|KR|P|C',
'AspNambic': 'Asp-N_ambic|DE||N',
'AspN': 'Asp-N|ND||N',
'CNBr': 'CNBr|M||C',
'Formic_acid': 'Formic_acid|D||NC',
'Chymotrypsin': 'Chymotrypsin|FYWL|P|C',
'LysCP': 'Lys-C/P|K||C',
'TrypsinP': 'Trypsin/P|KR||C',
'PepsinA': 'PepsinA|FL||C',
'V8E': 'V8-E|EQ|P|C',
'V8DE': 'V8-DE|EQND|P|C',
'TrypChymo': 'TrypChymo|FYWLKR|P|C',
}

H_plus = D('1.007276')
H2O = D('18.01056')

score_fields = ['DeNovoScore', 'MSGFScore', 'SpecProb', 'P-value', 'EFDR']

ArgC = ('(?P<ArgC>^((R\.[^P])|(_\.)))','(?P<ArgC>((R([+-]?\d+\.\d+)?\.[^P])|(\._))$)')
LysC = ('(?P<LysC>^((K\.[^P])|(_\.)))','(?P<LysC>((K([+-]?\d+\.\d+)?\.[^P])|(\._))$)')
Trypsin = ('(?P<Trypsin>^(([KR]\.[^P])|(_\.)))','(?P<Trypsin>(([KR]([+-]?\d+\.\d+)?\.[^P])|(\._))$)')
AspNambic = ('(?P<AspNambic>^((.\.[DE])|(_\.)))','(?P<AspNambic>\.[DE_]$)')
AspN = ('(?P<AspN>^((.\.[ND])|(_\.)))','(?P<AspN>\.[ND_]$)')
CNBr = ('(?P<CNBr>^[M_]\.)','(?P<CNBr>((M([+-]?\d+\.\d+)?\..)|(\._))$)')
Formic_acid = ('(?P<Formic_acid>^(([D_]\.)|(.\.D)))','(?P<Formic_acid>((D([+-]?\d+\.\d+)?\..)|(\.[D_]))$)')
Chymotrypsin = ('(?P<Chymotrypsin>^(([FYWL]\.[^P])|(_\..)))','(?P<Chymotrypsin>(([FYWL]([+-]?\d+\.\d+)?\.[^P])|(\._))$)')
LysCP = ('(?P<LysCP>^[K_]\.)','(?P<LysCP>((K([+-]?\d+\.\d+)?\..)|(\._))$)')
TrypsinP = ('(?P<TrypsinP>^[KR_]\.)','(?P<TrypsinP>(([KR]([+-]?\d+\.\d+)?\..)|(\._))$)')
PepsinA = ('(?P<PepsinA>^[FL_]\.)','(?P<PepsinA>(([FL]([+-]?\d+\.\d+)?\..)|(\._))$)')
V8E = ('(?P<V8E>^(([EQ]\.[^P])|(_\.)))','(?P<V8E>(([EQ]([+-]?\d+\.\d+)?\.[^P])|(\._))$)')
V8DE = ('(?P<V8DE>^(([NDEQ]\.[^P])|(_\.)))','(?P<V8DE>(([NDEQ]([+-]?\d+\.\d+)?\.[^P])|(\._))$)')
TrypChymo = ('(?P<TrypChymo>^(([FYWLKR]\.[^P])|(_\.)))','(?P<TrypChymo>(([FYWLKR]([+-]?\d+\.\d+)?\.[^P])|(\._))$)')

default_modifications = ['C2H3N1O1,C,fix,any,Carbamidomethylation']

def remove_file_extention(name):
    return '.'.join(name.split('.')[:-1])

def mass_of_peptide(peptide):
    """ Return neutral mass of peptide (without taking into account mass diffs)

    * peptide -- string with one-letter codes of aminoacids with optional mass diffs
    """
    return sum(aa2mass[aa] for aa in peptide if aa in aa2mass) + H2O

single_diff_pattern = re.compile("([" + "".join(atom2mass.keys()) + "])([-+]?\d+)?")
modified_aa_pattern = re.compile("(\w)([-+]?\d+\.\d+)")

def mass_of_molecule_diff(molecule_diff):
    """ Return mass of molecule difference (string of element names with numbers) """
    mass_diff = D('0')
    for single_diff in re.findall(single_diff_pattern, molecule_diff):
        atom = single_diff[0]
        number = single_diff[1] or 1
        number = int(number)
        mass_diff += atom2mass[atom] * number
    return mass_diff
	
def what_enzyme(enz_list,s_list,pep):
    """ Find enzymes which could have cleaved the protein
    
    * enz_list -- list of possible enzymes for this peptide
                  enzymes of previous peptides are taken into account
    * pep -- peptide sequence

    """
    enz_list_new = []
    semi_list_new = []
    for enz in s_list:
        if re.search(enz[0],pep) or re.search(enz[1],pep):
            semi_list_new.append(enz)
    for enz in enz_list:
        if re.search(enz[0],pep) and re.search(enz[1],pep):
            enz_list_new.append(enz)
        elif re.search(enz[0],pep) or re.search(enz[1],pep):
            semi_list_new.append(enz)
    return enz_list_new, semi_list_new

def find_modifications(modification, peptide):
    """ Return list of modification instance

    modification instance = (mass, mass_diff, aa_number, is_opt)
    mass of modification and numbers of affected residues: all and opt

    * modification -- string like "C2H3N1O1,C,fix,any,Carbamidomethylation".
        See "Mods.txt" for explanation
    * peptide -- see mass_of_peptide for explanation
    """
    molecule_diff, residues, mod_type, positions, name = modification.split(',')
    modification_mass = mass_of_molecule_diff(molecule_diff) \
        if molecule_diff[0].isalpha() else D(molecule_diff)
    result = []
    if residues == '*' and positions == 'any':
        raise Exception('Error in modifications file: * should not be "anywhere" modification')
    if mod_type.lower() == 'opt':
        for modified_aa in re.findall(modified_aa_pattern, peptide):
            aa, mass_diff_str = modified_aa
            mass_diff = D(mass_diff_str)
            if residues == '*' or aa in residues:
                if abs(mass_diff - modification_mass) <= D('0.001'):
                    peptide_before = [c for c in peptide.partition(aa+mass_diff_str)[0] if c.isalpha()]
                    aa_number = len(peptide_before) + 1
                    result.append((aa2mass[aa]+modification_mass, modification_mass, aa_number, True))
    elif mod_type.lower() == 'fix':
        for aa_number, aa in enumerate((c for c in peptide if c.isalpha()), 1):
            if residues == '*' or aa in residues:
                result.append((aa2mass[aa]+modification_mass, modification_mass, aa_number, False))
    else:
        raise Exception('Mods.txt: ModType can be "fix" or "opt"')
    return result

def read_modifications(modifications_file):
    """ Return (num_mods, modifications) from modifications file """
    num_mods = None
    modifications = []
    for line in modifications_file:
        line = line.partition('#')[0].strip()
        if 'NumMods' in line:
            if num_mods:
                raise Exception('Mods.txt: multiple NumMods')
            num_mods = line.partition('=')[2]
        elif len(line.split(',')) == 5:
            modifications.append(line)
    return (num_mods, modifications)

def modified_peptide_mass(modification_instances, peptide, num_mods=None):
    """ Return the final mass of peptide with all modifications
    
    * modification_instances -- list returned by find_modifications
    * peptide -- string representation of peptide, maybe with mass modifiers
        Example of mass modifier: "PEPTIDE+12.345SEQUENCE"
    * num_mods -- maximum allowable number of modifications in this peptide
        AssertationError is thrown if the real number of modifications is larger
    """
    mass = mass_of_peptide(peptide)
    real_num_mods = 0
    real_opt_residues = 0
    real_num_mods = len(modification_instances)
    for mass_aa, added_mass, aa_number, is_opt in modification_instances:
        real_opt_residues += is_opt
        mass += added_mass
    if num_mods is not None and real_num_mods > num_mods:
        raise Exception('Too many modifications in peptide %s' % peptide)
    if real_opt_residues != len(re.findall(modified_aa_pattern, peptide)):
        warn('Not all non-default optional modifications of peptide %s were found in the modifications file or no Mods.txt provided' % peptide, RuntimeWarning, 2)
    return mass

def xmlnsless_parse(xml_file):
    """ Parse xml_file, ignoring xml namespace, dummy for ElementTree.parse """
    stream = StringIO()
    stream.write(re.sub(r'xmlns=".*"', '', xml_file.read(), 1))
    stream.seek(0)
    return parse(stream)

def get_msrun(mzxml):
    """ Parse mzXML file abd return msRun element """
    mzxml_tree = xmlnsless_parse(mzxml)
    msrun = mzxml_tree.getroot().find('msRun')
    return msrun

xsd_duration_pattern = \
    re.compile(r"(?P<minus>-?)P((?P<years>[\d\.]+)Y)?((?P<months>[\d\.]+)M)?((?P<days>[\d\.]+)D)?"+
    r"(T((?P<hours>[\d\.]+)H)?((?P<minutes>[\d\.]+)M)?((?P<seconds>[\d\.]+)S)?)?")

def read_xsd_duration(string):
    """ Read time duration from xsd format and return number of seconds """
    d = re.match(xsd_duration_pattern, string).groupdict()
    for field in ['years', 'months', 'days', 'hours', 'minutes', 'seconds']:
        d[field] = D(d[field]) if d[field] else D()
    p = d['years']*365*24*3600 + d['months']*30*24*3600 + d['days']*24*3600
    p += d['hours'] * 3600 + d['minutes'] * 60 + d['seconds']
    if d['minus']:
        p = -p
    return p

def set_retention_times(msrun, msms_run_summary):
    """ Get retention times from msRun and set to children of msms_run_summary """
    scan2query = {}
    for spectrum_query in msms_run_summary.findall('spectrum_query'):
        scan2query[spectrum_query.get('start_scan')] = spectrum_query
    def get_scans(scan_container):
        for scan in scan_container.iterfind('scan'):
            yield scan
            for s in get_scans(scan):
                yield s
    for scan in get_scans(msrun):
        num = scan.get('num')
        if num in scan2query:
            spectrum_query = scan2query[num]
            retention_time_sec = str(read_xsd_duration(scan.get('retentionTime')))
            spectrum_query.set('retention_time_sec', retention_time_sec)

def create_msms_pipeline_analysis():
    msms_pipeline_analysis = Element('msms_pipeline_analysis')
    msms_pipeline_analysis.set('date', datetime.today().strftime("%Y-%m-%dT%H:%M:%S"))
    msms_pipeline_analysis.set('xmlns', "http://regis-web.systemsbiology.net/pepXML")
    msms_pipeline_analysis.set('xmlns:xsi', "http://www.w3.org/2001/XMLSchema-instance")
    msms_pipeline_analysis.set('xsi:schemaLocation',
        "http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v115.xsd")
    return msms_pipeline_analysis

def set_hit_ranks(msms_run_summary):
    """ Set hit_rank attribute to every search_hit """
    for spectrum_query in msms_run_summary.findall('spectrum_query'):
        search_result = spectrum_query.find('search_result')
        search_hits = search_result.findall('search_hit')
        search_hits.sort(key=lambda hit: D([score for score in hit.findall('search_score') \
            if score.get('name') == 'SpecProb'][0].get('value')))
        for search_hit in search_hits:
            search_result.remove(search_hit)
        for search_hit in search_hits:
            search_result.append(search_hit)
        for hit_rank, search_hit in enumerate(search_hits):
            search_hit.set('hit_rank', str(hit_rank + 1))

def read_msgfdb(in_file):
    """ Read output file of MS-GFDB, return list of dicts of fields """
    header = in_file.readline().strip().split()
    for line in in_file:
        line = line.strip()
        if line:
            parts = line.split('\t')
            f = dict(zip(header, parts))
            f['Scan#'] = int(f['Scan#'])
            f['Precursor'] = D(f['Precursor'])
            f['Charge'] = int(f['Charge'])
            yield f

def apply_msgfdb(in_file, msms_run_summary, modifications, num_mods):
    """ Read output file of MS-GFDB and add child elements to msms_run_summary """
    spectrum2element = {}
    enzyme_list = [ArgC,LysC,Trypsin,LysCP,Chymotrypsin,TrypChymo,TrypsinP,PepsinA,
                   CNBr,V8E,AspN,Formic_acid,AspNambic,V8DE]
    semi_list = []
    sample_enzyme = msms_run_summary.find('sample_enzyme')
    for f in read_msgfdb(in_file):
        spectrum = '%(name)s.%(scan)05i.%(scan)05i.%(charge)i' % \
            {'name': remove_file_extention(f['#SpecFile']),
             'scan': f['Scan#'], 'charge': f['Charge']}
        enzyme_list, semi_list = what_enzyme(enzyme_list, semi_list, f['Peptide'])
        peptide_prev_aa = f['Peptide'][0]
        if peptide_prev_aa == '_':
            peptide_prev_aa = '-'
        peptide_middle = f['Peptide'][2:-2]
        peptide_next_aa = f['Peptide'][-1]
        if peptide_next_aa == '_':
            peptide_next_aa = '-'
        if ' ' in f['Protein']:
            protein_name, protein_descr = f['Protein'].split(' ', 1)
        else:
            protein_name = f['Protein']
            protein_descr = ''
        precursor_neutral_mass = f['Precursor'] * f['Charge'] - f['Charge'] * H_plus

        if spectrum not in spectrum2element:
            spectrum_query = SubElement(msms_run_summary, 'spectrum_query')
            spectrum2element[spectrum] = spectrum_query
            spectrum_query.append(Element('search_result'))
            spectrum_query.set('spectrum', spectrum)
            spectrum_query.set('start_scan', str(f['Scan#']))
            spectrum_query.set('end_scan', str(f['Scan#']))
            spectrum_query.set('assumed_charge', str(f['Charge']))
            spectrum_query.set('precursor_neutral_mass', str(precursor_neutral_mass))

        spectrum_query = spectrum2element[spectrum]
        search_result = spectrum_query.find('search_result')
        search_hit = SubElement(search_result, 'search_hit')
        search_hit.set('peptide', "".join(aa for aa in peptide_middle if aa.isalpha()))
        search_hit.set('peptide_prev_aa', peptide_prev_aa)
        search_hit.set('peptide_next_aa', peptide_next_aa)
        search_hit.set('protein', protein_name)
        search_hit.set('protein_descr', protein_descr)

        modification_instances = sum((find_modifications(mod, peptide_middle) for mod in modifications), [])
        calc_neutral_pep_mass = modified_peptide_mass(modification_instances, peptide_middle, num_mods)
        if modification_instances:
            modification_info = SubElement(search_hit, 'modification_info')
            for mass, mass_diff, aa_number, is_opt in modification_instances:
                maam = SubElement(modification_info, 'mod_aminoacid_mass')
                maam.set('position', str(aa_number))
                maam.set('mass', str(mass))
        search_hit.set('calc_neutral_pep_mass', str(calc_neutral_pep_mass))
        search_hit.set('massdiff', str(precursor_neutral_mass - calc_neutral_pep_mass))
        for field in score_fields:
            if field in f:
                SubElement(search_hit, 'search_score', name=field, value=f[field])
#    sample_enzyme.set('fidelity',flag)
    if enzyme_list == []:
        if semi_list == []:
            sample_enzyme.set('name','NoEnzyme')
            sample_enzyme.set('fidelity','nonspecific')
        else:
            sample_enzyme.set('fidelity','semispecific')
            enzyme = re.split("\|",enzyme2name[re.search(r'<(\w+)>',semi_list[0][0]).group(1)])
    else:
        sample_enzyme.set('fidelity','specific')
        enzyme = re.split("\|",enzyme2name[re.search(r'<(\w+)>',enzyme_list[0][0]).group(1)])
    if not(enzyme_list == [] and semi_list == []):
        sample_enzyme.set('name',enzyme[0])
        specificity = SubElement(sample_enzyme, 'specificity')
        specificity.set('cut',enzyme[1])
        if enzyme[2]:
            specificity.set('no_cut',enzyme[2])
        specificity.set('sense',enzyme[3])

def msgfdb2pepxml(in_file, out_file, modifications_file=None, mzxml=None, fasta=None):
    """ Convert tab-separated ms-gfdb output to pepXML

    * in_file -- output file of MS-GFDB
    * out_file -- pepXML file to write to
    * modifications_file -- modifications file of MS-GFDB (Mods.txt)
    * mzxml -- input mzXML file of MS-GFDB
    * fasta -- input fasta file of MS-GFDB with database of peptides
    """
    if not out_file:
        out_filename = re.sub(r"\.msgfdb$", ".pep.xml", in_file.name)
        if out_filename == in_file.name:
            raise Exception("Provide output file or input file with extension .msgfdb")
        out_file = open(out_filename, 'w')
    modifications = default_modifications
    num_mods = None
    if modifications_file:
        num_mods, modifications = read_modifications(modifications_file)
    msms_pipeline_analysis = create_msms_pipeline_analysis()
    tree = ElementTree(msms_pipeline_analysis)
    msms_run_summary = SubElement(msms_pipeline_analysis, 'msms_run_summary')
    base_name = remove_file_extention(os.path.abspath(in_file.name))
    msms_run_summary.set('base_name', base_name)
    SubElement(msms_run_summary,'sample_enzyme')
    search_summary = SubElement(msms_run_summary, 'search_summary')
    search_summary.set('search_engine', "X! Tandem (k-score)")
    if fasta:
        search_database = SubElement(search_summary, 'search_database')
        search_database.set('local_path', os.path.abspath(fasta.name))
    apply_msgfdb(in_file, msms_run_summary, modifications, num_mods)
    set_hit_ranks(msms_run_summary)
    if mzxml:
        msrun = get_msrun(mzxml)
        ms_instrument = msrun.find('msInstrument')
        if ms_instrument is not None:
            for field in ['msManufacturer', 'msModel', 'msIonisation', 'msDetector']:
                param = ms_instrument.find(field)
                if param is not None:
                    msms_run_summary.set(field, param.get('value'))
        set_retention_times(msrun, msms_run_summary)
    out = StringIO()
    tree.write(out, encoding='UTF-8', xml_declaration=True)
    out_file.write(out.getvalue().replace('>', '>\n'))
    in_file.close()
    out_file.close()
    if modifications_file:
        modifications_file.close()
    if mzxml:
        mzxml.close()
    if fasta:
        fasta.close()

def main():
    p = argparse.ArgumentParser(description='MS-GFDB output to pepXML converter')
    p.add_argument('-i', help='MS-GFDB output', metavar='FILE', type=argparse.FileType('r'), required=True)
    p.add_argument('-c', help='MS-GFDB configuration (modifications)', metavar='FILE', type=argparse.FileType('r'))
    p.add_argument('-s', help='source mzXML', metavar='FILE', type=argparse.FileType('r'))
    p.add_argument('-o', help='pepXML file to write', metavar='FILE', type=argparse.FileType('w'))
    p.add_argument('-f', help='FASTA database file', metavar='FILE', type=argparse.FileType('r'), required=False)
    args = p.parse_args()
    print "Convert %s from MS-GFDB output format to pepXML file %s" % \
        (args.i.name, args.o.name if args.o else "")
    msgfdb2pepxml(args.i, args.o, args.c, args.s, args.f)

if __name__ == '__main__':
    try:
        main()
    except Exception, e:
        print e

