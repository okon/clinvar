import re
import sys
import gzip
import argparse
from collections import defaultdict
import xml.etree.ElementTree as ET
import functools
import operator
import pandas as pd

# functions
def replace_semicolons(s, replace_with=":"):
    return s.replace(";", replace_with)

def remove_newlines_and_tabs(s):
    return re.sub("[\t\n\r]", " ", s)

def reduce_list(l):
    return functools.reduce(operator.concat, l) 

def join_entries(d):
    for key,value in d.items():
        if type(value) == set:
            value_sorted=sorted(value)
            d[key]=remove_newlines_and_tabs(';'.join(map(replace_semicolons, value_sorted)))
        elif type(value) == list:
            d[key]=remove_newlines_and_tabs(';'.join(map(replace_semicolons, value)))
        if value is None:
            d[key]=''
    return d

def update_record_by_header(d, header):
    d_updated = {}
    for column in header:
        if column in d:
            d_updated[column]= d[column]
        else:
            d_updated[column]= ''
    return d_updated

def get_handle(file):
    if file[-3:] == '.gz':
        handle = gzip.open(file)
    else:
        handle = open(file)
    return handle

def write_to_out(dest, d=None, multi=None, incomplete=None, tally = None, measure_len=None, step='header'):
    # define file header
    HEADER = ['chrom', 'pos', 'ref', 'alt', 'start', 'stop', 'strand', 'variation_type', 'variation_id', 'rcv', 'scv',
        'allele_id', 'symbol',
        'hgvs_c', 'hgvs_n', 'hgvs_p', 'molecular_consequence',
        'clinical_significance', 'clinical_significance_ordered', 'pathogenic', 'likely_pathogenic',
        'uncertain_significance',
        'likely_benign', 'benign', 'review_status', 'review_status_ordered',
        'last_evaluated', 'all_submitters', 'submitters_ordered', 'all_traits',
        'all_pmids', 'inheritance_modes', 'age_of_onset', 'prevalence',
        'disease_mechanism', 'origin', 'xrefs', 'dates_ordered']
    if step == 'header':
        # Start writing
        dest.write(('\t'.join(HEADER) + '\n').encode('utf-8'))
        if multi is not None:
            multi.write(('\t'.join(HEADER) + '\n').encode('utf-8'))
        if incomplete is not None:
            incomplete.write(('\t'.join(HEADER) + '\n').encode('utf-8'))  
        # start tally
        tally = defaultdict(int)
        tally = {'single': 0,
                'multi': 0,
                'total': 0,
                'skipped':{
                'MeasureSet - multiple or none':0,
                'SequenceLocation - missing':0,
                'Other':0}}
        return tally
    elif step == 'records':
        # Convert all lists and sets to strings
        d = join_entries(d)
        # Convert record to list and match with header
        d = update_record_by_header(d, header=HEADER)
        
        if d['chrom'] != '' and measure_len == 1:
            dest.write(('\t'.join([str(d[column])  for column in HEADER]) + '\n').encode('utf-8'))
            tally['single'] += 1
        elif d['chrom'] != '' and measure_len > 1:
            if multi is not None:
                multi.write(('\t'.join([str(d[column]) for column in HEADER]) + '\n').encode('utf-8'))
                for key in ['allele_id','chrom','pos','ref','alt','hgvs_c','hgvs_n','hgvs_p']:
                    d[key]='' #reset values to get from next measure
            tally['multi'] += 1
        elif d['chrom'] == '' and measure_len == 0:
            if incomplete is not None:
                incomplete.write(('\t'.join([str(d[column])  for column in HEADER]) + '\n').encode('utf-8'))
            tally['skipped']['MeasureSet - multiple or none'] += 1
        elif d['chrom'] == '' and measure_len > 0:
            if incomplete is not None:            
                incomplete.write(('\t'.join([str(d[column])  for column in HEADER]) + '\n').encode('utf-8'))
            tally['skipped']['SequenceLocation - missing'] += 1
        else:
            if incomplete is not None:            
                incomplete.write(('\t'.join([str(d[column])  for column in HEADER]) + '\n').encode('utf-8'))
            tally['skipped']['Other'] += 1   
        
        tally['total'] += 1 
        
        assert tally['total'] == (tally['single'] + tally['multi'] + sum(tally['skipped'].values()))
        
        if tally['single'] % 100 == 0:
            dest.flush()
        if multi is not None and tally['multi'] % 100 == 0:
            multi.flush()
        if incomplete is not None and sum(tally['skipped'].values()) % 100 == 0:
            incomplete.flush()
            
        return (tally, d)

def get_accession(elem, field=None, record_id=None):
    acc_list = [x.find('./ClinVarAccession') for x in elem]
    try: 
        acc_set = set()
        for acc in acc_list:
            if acc.attrib.get('Type') == field:
                acc_set.add(acc.attrib.get('Acc'))
        #return ';'.join(acc_set)  
        return acc_set
    except:
        print("Record", record_id, "has no", field, "record")

def get_measureset(elem, field=None, record_id=None):
    measureset= elem.findall('.//MeasureSet')
    if len(measureset) == 1:
        if field == 'id':
            return measureset[0].attrib.get('ID')
        elif field == 'type':
            return measureset[0].get('Type')
        elif field == 'measure':
            return measureset[0].findall('.//Measure')
        if field == 'symbol':
            try: 
                var_name = measureset[0].find(".//Name/ElementValue").text
                match = re.search(r"\(([A-Za-z0-9]+)\)", var_name) # not good with regex, so leaving it...
                genesymbol = match.group(1)
                return(genesymbol)
            except:
                return ''
    elif len(measureset) > 1:
        print("Record", record_id, "has more than one measure set.")
        return ''
    elif len(measureset) == 0:
        print("Record", record_id, "has no measure set type")
        return ''
       
def get_pmids(elem):
    # group(1) will be all the text after the word PubMed or PMID
    mentions_pubmed_regex = '(?:PubMed|PMID)(.*)' 
    # group(1) will be the first PubMed ID, group(2) will be all remaining text 
    extract_pubmed_id_regex = '[^0-9]+([0-9]+)[^0-9](.*)'  
    pmids =set()
    for citation in elem.findall('.//Citation'):
        for id_node in citation.findall('.//ID'): 
            if id_node.attrib.get('Source') == 'PubMed':
                pmids.add(id_node.text)
    for comment in elem.findall('.//Comment'):
        mentions_pubmed = re.search(mentions_pubmed_regex, comment.text)
        if mentions_pubmed is not None and mentions_pubmed.group(1) is not None:
            remaining_text = mentions_pubmed.group(1)
            while True:
                pubmed_id_extraction = re.search(extract_pubmed_id_regex, remaining_text)
                if pubmed_id_extraction is None:
                    break
                elif pubmed_id_extraction.group(1) is not None:
                    pmids.add(pubmed_id_extraction.group(1))
                    if pubmed_id_extraction.group(2) is not None:
                        remaining_text = pubmed_id_extraction.group(2)
    return(pmids)

def get_submitters(elem, ordered=True, record_id=None):
    # find any/all submitters
    submitters= []
    for submitter_node in elem.findall('.//ClinVarSubmissionID'):
        try:
            submitters.append(submitter_node.attrib['submitter'].replace(';', ','))
        except:
            print("Record", record_id, "has no submitters")
    if ordered:
        return submitters 
    else:
        return set(submitters)

def get_clinsig(elem, field=None, count=False, record_id=None):
    clinsig_list = reduce_list([x.findall('./ClinicalSignificance') for x in elem])
    if field == "status": 
        review_status = ''
        if len(clinsig_list) == 0:
            print("Record", record_id, "has ClinicalSignificance in", elem[0].tag)
        else:
            try:
                review_status = clinsig_list[0].find('./ReviewStatus').text
            except:
                print("Record", record_id, "has no ReviewStatus field in", elem[0].tag)                
        return review_status
    elif field == "status_ordered":
        status_ordered = []
        if len(clinsig_list) == 0:
            print("Record", record_id, "has ClinicalSignificance in", elem[0].tag)
        else:
            try:
                status_ordered = [x.find('./ReviewStatus').text for x in clinsig_list]
            except:
                print("Record", record_id, "has no ReviewStatus field in", elem[0].tag)                  
        return status_ordered
    elif field == "description":
        desc_ordered = []
        if len(clinsig_list) == 0:
            print("Record", record_id, "has no ClinicalSignificance in", elem[0].tag)
        else:
            try:   
                desc_ordered = [x.find('./Description').text for x in clinsig_list]
            except:
                print("Record", record_id, "has Description field in", elem[0].tag)  
        if count:
            desc_ordered_low = [x.lower() for x in desc_ordered] #needs to be lower case
            desc_count = {}
            for k in ['pathogenic','likely_pathogenic','uncertain_significance','benign','likely_benign']:
                k_count = k.replace('_',' ') #key for count
                desc_count[k] = str(desc_ordered_low.count(k_count))
            return desc_count
        elif not count:
            return desc_ordered
    elif field == "last_eval":
        if len(clinsig_list) == 0:
            print("Record", record_id, "has ClinicalSignificance in", elem[0].tag)
        else:            
            date_ordered = [x.attrib.get('DateLastEvaluated', '0000-00-00') for x in clinsig_list]
            if len(date_ordered) == 0:
                return '0000-00-00'
            else:
                return date_ordered

def get_traits(elem, field=None):
    # now find the disease(s) this variant is associated with
    for traitset in elem.findall('.//TraitSet'):
        if field == "traits":
            trait_values = []
            disease_name_nodes = traitset.findall('.//Name/ElementValue')
            trait_values = [x.text for x in disease_name_nodes if x.attrib.get('Type') == 'Preferred']
            return trait_values
        elif field == "mechanism":
            attribute_nodes = traitset.findall('.//AttributeSet/Attribute')
            disease_mechanism = {x.text.strip() for x in attribute_nodes if x.attrib.get('Type') == 'disease mechanism'}
            return(disease_mechanism)
        elif field == "xrefs":
            # put all the cross references one column, it may contains NCBI gene ID, conditions ID in disease databases.
            xref_nodes = traitset.findall('.//XRef')
            xrefs = {f"{x.attrib.get('DB')}:{x.attrib.get('ID')}".replace(';',',') for x in xref_nodes}
            return(xrefs)
        
def get_disease_attrib(elem, field=None):
    att_set=set()
    for node in elem:
        for node_elem in node.findall('./AttributeSet/Attribute'):
            if ((field == "inheritance" and node_elem.attrib.get('Type') == 'ModeOfInheritance') or 
                (field == "age" and node_elem.attrib.get('Type') == 'AgeOfOnset')): 
                att_set.add(node_elem.text.strip())
    return(att_set)

def get_origin(elem):
    for node in elem:
        origin_set = {x for x in node.findall('./ObservedIn/Sample/Origin')}
    return (origin_set)

def get_measure_info(elem, field=None, current_symbol=None):
    if field == "symbol":
        symbol = ''
        for node in elem.findall('.//Symbol'):
             if node.find('ElementValue').attrib.get('Type') == 'Preferred':
                 symbol = node.find('ElementValue').text
        if symbol == '':
            return current_symbol
        else:
            return symbol
    elif field == "allele_id":  # find the allele ID (//Measure/@ID)
        att_id = elem.attrib.get('ID')
        return att_id
        
def get_measure_loc(elem, symbol=None, build='GRCh37', record_id=None):
    loc_list = []
    measure_symbol = ''
    lookup_dict={'chrom':['Chr'],
    'pos':['positionVCF','start'],
    'ref':['referenceAlleleVCF','referenceAllele'],
    'alt':['alternateAlleleVCF','alternateAllele'],
    'start':['start'],
    'stop':['stop'],
    'strand':['Strand'],
    'Assembly':['Assembly'],
    'MeasureRelationship':['MeasureRelationship']}
    for node in elem:
        if node.tag == 'MeasureRelationship':
            measure_symbol = node.find(".//Symbol/ElementValue").text
            try:
                temp_dict=node.find(".//SequenceLocation").attrib
                temp_dict['MeasureRelationship']='yes'
                loc_list.append(temp_dict)
            except:
                print(f"Record {record_id} has no sequence locations")
        if node.tag == 'SequenceLocation':
            temp_dict=node.attrib
            temp_dict['MeasureRelationship']='no'
            loc_list.append(node.attrib)
    loc_df = pd.DataFrame(loc_list)
    for key,values in lookup_dict.items():
        for value in values:
            if value in loc_df.columns and key not in loc_df.columns:
                loc_df[key]=loc_df[value]
        if key not in loc_df.columns:
            loc_df[key]=''
    loc_df = loc_df.fillna('')
    full_loc_df = loc_df[(loc_df.MeasureRelationship == 'no') & 
                        (loc_df.Assembly == build) & 
                        (loc_df.chrom != '') &
                        (loc_df.pos != '') &
                        (loc_df.ref != '') &
                        (loc_df.alt != '')]
    if full_loc_df.empty and loc_df.empty:
        print(f"Record {record_id} has no sequence locations")
        return full_loc_df.to_dict('r')
    elif full_loc_df.empty and not loc_df.empty:
        print(f"Record {record_id} has no full sequence locations with chrom, pos, ref and alt")
        print(loc_df.to_string())
        return full_loc_df.to_dict('r')
    elif full_loc_df.shape[0] > 1:
        try:
            full_loc_df.duplicated(subset=['chrom','pos','ref','alt'], keep=False).all() 
        except:
            print(f"Record {record_id} chrom positions don't match, taking first one")    
    full_loc_df = full_loc_df.iloc[0]
    if symbol and symbol == measure_symbol:
        try:
            full_loc_df.strand = loc_df[(loc_df.strand != '') & 
                                        (loc_df.Accession == full_loc_df.Accession) &
                                        (loc_df.MeasureRelationship == 'yes') ].strand.values[0]
        except:
            pass
    return full_loc_df.to_dict()

def get_measure_cons(elem):
    d=defaultdict(set)
    lookup_dict = {'HGVS, coding, RefSeq':['c.','hgvs_c'],
                    'HGVS, non-coding':['n.','hgvs_n'],
                    'HGVS, protein, RefSeq':['p.','hgvs_p']}
    for node in elem.findall('./AttributeSet'):
        att_type=node.find('./Attribute').attrib.get('Type')
        att_val=node.find('./Attribute').text
        for key,value in lookup_dict.items():
            if key in att_type and value[0] in att_val:
                d[value[1]].add(att_val)
        if att_type == 'MolecularConsequence':
            for xref in node.findall('.//XRef'):
                if xref.attrib.get('DB') == "RefSeq":
                    d['molecular_consequence'].add(":".join([xref.attrib.get('ID'), att_val]))
    return d

def parse_clinvar_tree(handle, dest_file , multi_file=None, incomplete_file=None, log_file=None, genome_build='GRCh37'):
    """Parse clinvar XML
    Args:
        handle: Open input file handle for reading the XML data
        dest_file: Output file for simple variants
        multi_file: Output file for complex non-single-variant clinvar records
            (eg. compound het, haplotypes, etc.)
        incomplete_file: Output file for incomplete clinvar records (that could not be completely parsed)
        log_file: Log file
        genome_build: Either 'GRCh37' or 'GRCh38'
    """
    # open files for writing
    dest = open(dest_file,'wb')

    if multi_file is not None:
        multi= open(multi_file,'wb')
    if incomplete_file is not None:
        incomplete = open(incomplete_file,'wb')
    if log_file is not None:    
        log = open(log_file, "w")
        sys.stdout = log

    # establish tally and write out header
    tally = write_to_out(dest=dest,multi=multi, incomplete=incomplete,step='header')

    # get an iterable
    context = ET.iterparse(handle, events=("start", "end"))

    # temp_count=0 #troubleshooting
    for event, elem in context:
        if elem.tag != 'ClinVarSet' or event != 'end':
            continue
        
        # temp_count +=1
        # if temp_count % 100 == 0:
        #     print("Number of records processed", temp_count)
        # if temp_count < 1000: # troubleshooting
        #     continue

        root = elem
        current_row=defaultdict()
        current_id=root.attrib.get('ID')
        # Extract all ReferenceClinVarAssertion and ClinVarAssertion fields
        assert root.findall("./ReferenceClinVarAssertion"), "No ReferenceClinVarAssertion field"
        assert root.findall("./ClinVarAssertion"), "No ClinVarAssertion field"
        rcva=root.findall("./ReferenceClinVarAssertion")
        cva=root.findall("./ClinVarAssertion")
        
        # Pull out relevant fields from this ClinVarSet record
        current_row['rcv'] = get_accession([rcva[0]], field='RCV', record_id=current_id) 
        current_row['variation_id'] = get_measureset(rcva[0], field='id', record_id=current_id)
        current_row['variation_type'] = get_measureset(rcva[0], field='type', record_id=current_id)
        current_row['scv'] = get_accession(cva, field='SCV')
        current_row['all_pmids'] = get_pmids(root)
        current_row['submitters_ordered'] = get_submitters(root, record_id=current_id)
        current_row['all_submitters'] = get_submitters(root, ordered=False, record_id=current_id)
        current_row['review_status'] = get_clinsig([rcva[0]], field='status', record_id=current_id)
        current_row['clinical_significance'] = get_clinsig([rcva[0]], field='description', record_id=current_id)
        current_row['last_evaluated'] = get_clinsig([rcva[0]], field='last_eval', record_id=current_id)
        current_row['review_status_ordered'] = get_clinsig(cva, field='status_ordered', record_id=current_id)
        current_row['clinical_significance_ordered'] = get_clinsig(cva, field='description', record_id=current_id)
        current_row.update(get_clinsig(cva, 
                                        field='description',
                                        count=True,
                                        record_id=current_id)) 
        current_row['dates_ordered'] = get_clinsig(cva, field='last_eval', record_id=current_id) 
        current_row['all_traits'] = get_traits(root, field='traits') 
        current_row['disease_mechanism'] = get_traits(root, field='mechanism') # can't find multiple yet
        current_row['xrefs'] = get_traits(root, field='xrefs')
        current_row['inheritance_modes'] = get_disease_attrib(cva, field='inheritance')  # can't find multiple yet
        current_row['age_of_onset'] = get_disease_attrib(cva, field='age') # can't find any, need to check
        current_row['symbol'] = get_measureset(rcva[0], field='symbol', record_id=current_id)
        measures = get_measureset(rcva[0], field='measure', record_id=current_id)
        for measure in measures:
            if measure != '':        
                current_row['symbol'] = get_measure_info(measure, field='symbol', current_symbol=current_row['symbol'])
                current_row['allele_id'] = get_measure_info(measure, field='allele_id')
                current_row.update(get_measure_loc(measure,
                                        symbol=current_row['symbol'],
                                        build=genome_build,
                                        record_id=current_id))
                current_row.update(get_measure_cons(measure))
                    
            # Convert all lists and sets to strings
            tally,current_row = write_to_out(d=current_row,dest=dest, multi=multi,incomplete=incomplete,
                                                tally=tally, measure_len=len(measures),step='records')
            
            if tally['total'] % 100 == 0:
                print("Number of records processed", tally['total'])
                
            
        root.clear()
        
    print('total records processed:', tally['total'])
    print('single measure records:',tally['single']) 
    print('multi measure records (each measure counted as record):',tally['multi']) 
    for key,value in tally['skipped'].items():
        print("skipped due to", key,":",value)
    
    dest.close()
    multi.close()
    incomplete.close()
    log.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract ClinVar records from the ClinVar XML dump')
    parser.add_argument('-g', '--genome-build', choices=['GRCh37', 'GRCh38'],
                        help='Genome version (either GRCh37 or GRCh38)', required=True)
    parser.add_argument('-x', '--xml', dest='xml_file', type=str, help='ClinVar XML dump file (can be .gz)', required=True)
    parser.add_argument('-o', '--out', type=str, help="Output file name for single alleles",required=True)
    parser.add_argument('-m', '--multi',type=str, help="Output file name for complex alleles")
    parser.add_argument('-i', '--incomplete',type=str, help="Output file name for incompletely parsed records")
    parser.add_argument('-l', '--log',type=str, help="Log file name")
    
    args = parser.parse_args()
    parse_clinvar_tree(get_handle(args.xml_file), dest_file=args.out, 
                       multi_file=args.multi, incomplete_file=args.incomplete,
                       log_file=args.log, genome_build=args.genome_build)


