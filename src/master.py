"""
Alternate implementation of master.bash with improved logging, skipping of commands if output files are up-to-date, parallelization, etc.

Run with -h to see all options.
"""

import configargparse
from datetime import datetime
import ftplib
import os
import sys
from distutils import spawn

try:
    import pypez
    import pysam
    import pandas   # make sure all dependencies are installed
except ImportError as e:
    sys.exit("ERROR: Python module not installed. %s. Please run 'pip install -r requirements.txt' " % e)
for executable in ['wget', 'Rscript', 'tabix', 'vt']:
    assert spawn.find_executable(executable), "Command %s not found, see README" % executable

p = configargparse.getArgParser()
g = p.add_argument_group('main args')
g.add("--b37-genome", help="b37 .fa genome reference file", required=True)
g.add("--b38-genome", help="b38 .fa genome reference file. NOTE: chromosome names must be like '1', '2'.. 'X', 'Y', 'MT'.", required=True)
g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, grab the latest from NCBI.")
g.add("-S", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, grab the latest from NCBI.")
g.add("-E", "--exac-sites-vcf",  help="ExAC sites vcf file. If specified, a clinvar table with extra ExAC fields will also be created.")

pypez.init_command_line_args()
args = p.parse_args()

reference_genomes = {'b37': args.b37_genome, 'b38': args.b38_genome}
clinvar_xml = args.clinvar_xml
exac_sites_vcf = args.exac_sites_vcf
clinvar_variant_summary_table = args.clinvar_variant_summary_table

for key, path in reference_genomes.items():
    if not os.path.isfile(path):
        p.error("%s genome reference: file not found: %s" % (key, path))

if exac_sites_vcf:
    if not os.path.isfile(exac_sites_vcf):
        p.error("ExAC sites vcf: file not found: %s" % exac_sites_vcf)
    if not os.path.isfile(exac_sites_vcf + ".tbi"):
        p.error("ExAC sites vcf: tabix index not found: %s.tbi" % exac_sites_vcf)


def get_remote_file_changed_time(ftp_host, ftp_path):
    """Returns time modified in seconds since the epoch"""

    print("Retrieving last-changed time for %s" % os.path.join(ftp_host, ftp_path))
    try:
        ftp = ftplib.FTP(ftp_host, timeout=3)  # timeout = 3 seconds
        ftp.login()
        response = ftp.sendcmd("MDTM " + ftp_path)
        last_changed_time = datetime.strptime(response[4:], "%Y%m%d%H%M%S")
        return int(last_changed_time.strftime("%s"))  #.strftime("%d %B %Y %H:%M:%S")
    except Exception as e:
        print("ERROR: retrieving last-changed time for %s: %s" % (os.path.join(ftp_host, ftp_path), e))
        return 0


def download_if_changed(job_runner, local_path, ftp_host, ftp_path):
    remote_changed_time = get_remote_file_changed_time(ftp_host, ftp_path)
    local_changed_time = os.path.getmtime(local_path) if os.path.isfile(local_path) else 0
    ftp_address = "ftp://%s/%s" % (ftp_host, ftp_path)
    if remote_changed_time > local_changed_time:
        print("Local copy of %s is out of date. The remote version changed on %s" % (ftp_address, datetime.fromtimestamp(remote_changed_time)))
        job_runner.add_parallel(pypez.Job("wget %s -O OUT:%s" % (ftp_address, local_path)))
    else:
        print("Local copy of %s is up to date. The remote version hasn't changed since %s" % (ftp_address, datetime.fromtimestamp(remote_changed_time)))
        #ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

jr = pypez.JobRunner()

if clinvar_xml:
    if not os.path.isfile(clinvar_xml):
        p.error("ClinVar XML specified but not found: %s" % clinvar_xml)
    if not clinvar_xml.endswith('.gz'):
        p.error("ClinVar XML expected to be gzipped: %s" % clinvar_xml)
    clinvar_xml = clinvar_xml
else:
    print("Checking for new clinvar release")
    clinvar_xml = "ClinVarFullRelease_00-latest.xml.gz"
    download_if_changed(jr, clinvar_xml,  "ftp.ncbi.nlm.nih.gov", "/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz")

if clinvar_variant_summary_table:
    if not os.path.isfile(clinvar_variant_summary_table):
        p.error("ClinVar variant summary table specified but not found: %s" % clinvar_variant_summary_table)
    if not clinvar_variant_summary_table.endswith('.gz'):
        p.error("ClinVar variant summary table expected to be gzipped: %s" % clinvar_variant_summary_table)
    variant_summary_table = clinvar_variant_summary_table
else:
    print("Checking for new clinvar release")
    variant_summary_table = "variant_summary.txt.gz"
    download_if_changed(jr, variant_summary_table,  "ftp.ncbi.nlm.nih.gov", "/pub/clinvar/tab_delimited/variant_summary.txt.gz")

jr.run()

job = pypez.Job()

# normalize (convert to minimal representation and left-align)
# the normalization code is in a different repo (useful for more than just clinvar) so here I just wget it:
job.add("wget -N https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py")

for genome_build in ('b37', 'b38'):
    # extract the GRCh37 coordinates, mutant allele, MeasureSet ID and PubMed IDs from it. This currently takes about 20 minutes.
    genome_build_id = genome_build.replace('b', 'GRCh')
    reference_genome = reference_genomes[genome_build]

    job.add(("python -u IN:parse_clinvar_xml.py "
            "-x IN:%(clinvar_xml)s "
            "-g %(genome_build_id)s "
            "-o OUT:clinvar_table_raw.single.%(genome_build)s.tsv "
            "-m OUT:clinvar_table_raw.multi.%(genome_build)s.tsv") % locals())

    for is_multi in (True, False):  # multi = clinvar submission that describes multiple alleles (eg. compound het, haplotypes, etc.)
        single_or_multi = 'multi' if is_multi else 'single'
        fsuffix = "%(single_or_multi)s.%(genome_build)s" % locals() # file suffix
        output_dir = '../output/%(genome_build)s/%(single_or_multi)s' % locals()
        os.system('mkdir -p ' + output_dir)

        job.add("python -u normalize.py -R IN:%(reference_genome)s < IN:clinvar_table_raw.%(fsuffix)s.tsv > OUT:clinvar_table_normalized.%(fsuffix)s.tsv" % locals())

        #remove empty rows
        #job.add("ex -s +'bufdo!v/\S/d' -cxa clinvar_table_normalized.tsv") #remove the empty lines

        #sort
        job.add(("(cat IN:clinvar_table_normalized.%(fsuffix)s.tsv | head -1 > OUT:clinvar_allele_trait_pairs.%(fsuffix)s.tsv ) && " + # header row
            "(cat IN:clinvar_table_normalized.%(fsuffix)s.tsv | tail -n +2 | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> OUT:clinvar_allele_trait_pairs.%(fsuffix)s.tsv ) && " + # numerically sort chroms 1-22
            "(cat IN:clinvar_table_normalized.%(fsuffix)s.tsv | tail -n +2 | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> OUT:clinvar_allele_trait_pairs.%(fsuffix)s.tsv )") % locals())   # lexicogaraphically sort non-numerical chroms at end

        job.add("cat IN:clinvar_allele_trait_pairs.%(fsuffix)s.tsv | grep -v ^$ | bgzip -c > OUT:clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz" % locals())
        job.add("tabix -S 1 -s 1 -b 2 -e 2 IN:clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz" % locals(), output_filenames=["clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz.tbi" % locals()])
        job.add("cp IN:clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz IN:clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz.tbi %(output_dir)s/" % locals(), output_filenames=[
                "%(output_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz" % locals(),
                "%(output_dir)s/clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz.tbi" % locals()
                ])

        # join information from the tab-delimited summary to the normalized genomic coordinates
        job.add(("Rscript IN:join_data.R "
                    "IN:%(variant_summary_table)s "
                    "IN:clinvar_allele_trait_pairs.%(fsuffix)s.tsv "
                    "OUT:clinvar_combined.%(fsuffix)s.tsv") % locals())

        # now sort again by genomic coordinates (because R's merge function ruins this)
        # and group by allele
        job.add(("(cat IN:clinvar_combined.%(fsuffix)s.tsv | head -1 > OUT:clinvar_sorted.%(fsuffix)s.tsv ) && " + # header row
            "(cat IN:clinvar_combined.%(fsuffix)s.tsv | tail -n +2 | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> OUT:clinvar_sorted.%(fsuffix)s.tsv ) && " + # numerically sort chroms 1-22
            "(cat IN:clinvar_combined.%(fsuffix)s.tsv | tail -n +2 | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> OUT:clinvar_sorted.%(fsuffix)s.tsv )") % locals())   # lexicogaraphically sort non-numerical chroms at end

        job.add("python -u IN:group_by_allele.py -i IN:clinvar_sorted.%(fsuffix)s.tsv | tee clinvar_alleles.%(fsuffix)s.tsv | bgzip -c > OUT:clinvar_alleles.%(fsuffix)s.tsv.gz" % locals())


        job.add("tabix -S 1 -s 1 -b 2 -e 2 IN:clinvar_alleles.%(fsuffix)s.tsv.gz" % locals(),
                output_filenames=[
                    "clinvar_alleles.%(fsuffix)s.tsv.gz.tbi" % locals()
                ])
        job.add("cp IN:clinvar_alleles.%(fsuffix)s.tsv.gz IN:clinvar_alleles.%(fsuffix)s.tsv.gz.tbi %(output_dir)s/"  % locals(),
                output_filenames=[
                    "%(output_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz" % locals(),
                    "%(output_dir)s/clinvar_alleles.%(fsuffix)s.tsv.gz.tbi" % locals()
                ])


        # create vcf
        job.add(("python -u IN:clinvar_table_to_vcf.py IN:clinvar_alleles.%(fsuffix)s.tsv | "
                    "bgzip -c > OUT:clinvar_alleles.%(fsuffix)s.vcf.gz") % locals())  # create compressed version
        job.add("tabix IN:clinvar_alleles.%(fsuffix)s.vcf.gz" % locals(), output_filenames=["clinvar_alleles.%(fsuffix)s.vcf.gz.tbi" % locals()])
        job.add("cp IN:clinvar_alleles.%(fsuffix)s.vcf.gz IN:clinvar_alleles.%(fsuffix)s.vcf.gz.tbi %(output_dir)s/" % locals(), output_filenames=[
            "%(output_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz" % locals(),
            "%(output_dir)s/clinvar_alleles.%(fsuffix)s.vcf.gz.tbi" % locals()])

        # create uncompressed example files that contain the 1st 750 lines of the compressed tsvs so people can easily see typical values online on github
        job.add("gunzip -c IN:clinvar_alleles.%(fsuffix)s.vcf.gz | head -n 750 > OUT:%(output_dir)s/clinvar_alleles_example_750_rows.%(fsuffix)s.vcf" % locals())
        job.add("gunzip -c IN:clinvar_alleles.%(fsuffix)s.tsv.gz | head -n 750 > OUT:%(output_dir)s/clinvar_alleles_example_750_rows.%(fsuffix)s.tsv" % locals())
        job.add("gunzip -c IN:clinvar_allele_trait_pairs.%(fsuffix)s.tsv.gz | head -n 750 > OUT:%(output_dir)s/clinvar_allele_trait_pairs_example_750_rows.%(fsuffix)s.tsv" % locals())

        # create tsv table with extra fields from ExAC: filter, ac_adj, an_adj, popmax_ac, popmax_an, popmax
        if exac_sites_vcf and genome_build == "GRCh37":
            normalized_vcf = os.path.basename(exac_sites_vcf).split('.vcf')[0] + ".normalized.vcf.gz" % locals()
            job.add(("vt decompose -s IN:%(exac_sites_vcf)s | "
                     "vt normalize -r IN:%(reference_genome)s - | "
                     "bgzip -c > OUT:%(normalized_vcf)s") % locals())
            job.add("tabix IN:%(normalized_vcf)s" % locals(), output_filenames=["%(normalized_vcf)s.tbi" % locals()])
            job.add(("python -u IN:add_exac_fields.py -i IN:clinvar_alleles.%(fsuffix)s.tsv -e IN:%(normalized_vcf)s | "
                        "bgzip -c > OUT:clinvar_alleles_with_exac.%(fsuffix)s.tsv.gz") % locals())
            job.add("tabix -S 1 -s 1 -b 2 -e 2 IN:clinvar_alleles_with_exac.%(fsuffix)s.tsv.gz" % locals(), output_filenames=["clinvar_alleles_with_exac.%(fsuffix)s.tsv.gz.tbi" % locals()])
            job.add("cp IN:clinvar_alleles_with_exac.%(fsuffix)s.tsv.gz IN:clinvar_alleles_with_exac.%(fsuffix)s.tsv.gz.tbi %(output_dir)s/" % locals(), output_filenames=[
                "%(output_dir)s/clinvar_alleles_with_exac.%(fsuffix)s.tsv.gz" % locals(),
                "%(output_dir)s/clinvar_alleles_with_exac.%(fsuffix)s.tsv.gz.tbi" % locals()])

            job.add("gunzip -c IN:clinvar_alleles_with_exac.%(fsuffix)s.tsv.gz | head -n 750 > OUT:%(output_dir)s/clinvar_alleles_with_exac_example_750_rows.%(fsuffix)s.tsv" % locals())

        # create a stats file that summarizes some of the columns of clinvar_alleles.tsv.gz file
        # Columns: 1: chrom, 2: pos, 3: ref, 4: alt, 5: measureset_type, 6: measureset_id, 7: rcv, 8: allele_id,
        # 9: symbol, 10: hgvs_c, 11: hgvs_p, 12: molecular_consequence, 13: clinical_significance, 14: pathogenic, 15: benign, 16: conflicted, 17: review_status,
        # 18: gold_stars, 19: all_submitters, 20: all_traits, 21: all_pmids, 22: inheritance_modes,
        # 23: age_of_onset, 24: prevalence, 25: disease_mechanism, 26: origin, 27: xrefs

        python_oneliner_to_format_stats = "import sys;print(', '.join(['%s: %s' % (i+1, v) for l in sys.stdin for i,v in enumerate(l.split())]))"
        job.add("""echo \
        Columns: $(gunzip -c clinvar_alleles.%(fsuffix)s.tsv.gz | \
            head -n 1 | \
            python -c "%(python_oneliner_to_format_stats)s") > \
            OUT:clinvar_alleles_stats.%(fsuffix)s.txt &&
        echo ================ >> OUT:clinvar_alleles_stats.%(fsuffix)s.txt &&
        echo Total Rows: $(gunzip -c clinvar_alleles.%(fsuffix)s.tsv.gz | tail -n +2 | wc -l) >> OUT:clinvar_alleles_stats.%(fsuffix)s.txt &&
        for i in 13 14 15 16 17 18 22 23 24 25 26; do
            echo ================ >> OUT:clinvar_alleles_stats.%(fsuffix)s.txt ;
            gunzip -c IN:clinvar_alleles.%(fsuffix)s.tsv.gz | head -n 1 | cut -f $i >> OUT:clinvar_alleles_stats.%(fsuffix)s.txt ;
            gunzip -c IN:clinvar_alleles.%(fsuffix)s.tsv.gz | tail -n +2 | cut -f $i | tr ';' '\n' | sort | uniq -c | sort -r -n >> OUT:clinvar_alleles_stats.%(fsuffix)s.txt ;
        done
        """ % locals(), input_filenames=["clinvar_alleles.%(fsuffix)s.tsv.gz" % locals(), "master.py"])

        job.add("cp IN:clinvar_alleles_stats.%(fsuffix)s.txt OUT:%(output_dir)s/clinvar_alleles_stats.%(fsuffix)s.txt" % locals())

# run the above commands
jr.run(job)
