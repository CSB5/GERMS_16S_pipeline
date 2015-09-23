#!/usr/bin/env python
"""This is a remake of Vinutha's original 16S pipeline, which was
designed for Illumina shotgun Sequencing of 16S rRNA Amplicon
Sequences (Ong et al., 2013, PMID 23579286). At its core it's running
EMIRGE (Miller et al., 2011, PMID 21595876) or EMIRGE amplicon for
reconstructing the sequences and BLAST + Greengenes for classification.
"""


import os
import sys
import argparse
import logging
import copy
import json
import shutil
import subprocess
import getpass

try:
    # python 3
    from itertools import zip_longest
except ImportError:
    # python 2
    from itertools import izip_longest as zip_longest


__author__ = "Andreas Wilm"
__version__ = "2.0.0a"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "The MIT License (MIT)"


# global logger
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
    format='%(levelname)s [%(asctime)s]: %(message)s')


# The wrapper script created here calling snakemake
SNAKEMAKE_CLUSTER_WRAPPER = "snake.sh"
# root name of snakemake's make file
SNAKEMAKE_FILE = "snake.make"
# path to snakemake's makefile template
SNAKEMAKE_TEMPLATE = os.path.join(
    os.path.dirname(sys.argv[0]), SNAKEMAKE_FILE)
# root name of config file loaded in SNAKEMAKE_FILE


CONFIG_FILE = "conf.json"
# template variables to be written to CONFIG_FILE
# FIXME might be better to have a template as well and add user variables only
CONF = dict()
# databases and files
CONF['EMIRGE_BASEDIR'] = "/mnt/software/stow/emirge-v0.60-15-g0ddae1c-wilma/bin/"
#CONF['SSU_FA'] = '/mnt/genomeDB/misc/softwareDB/emirge/SSU_candidate_db.fasta'
CONF['SSU_FA'] = '/mnt/projects/wilma/16s/S16_V2/ssu/SSU_candidate_db.p-knowlesi-spikein.fasta'
CONF['SSU_DB'] = CONF['SSU_FA'].replace('.fasta', '')
CONF['SPIKEIN-NAME'] = "Plasmodium.knowlesi.profilin"
CONF['GG_REF'] = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus/rep_set/99_otus.fasta'
CONF['GG_REF'] = '/mnt/projects/wilma/16s/S16_V2/greengenes/99_otus.fasta'
CONF['GG_TAX'] = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus/taxonomy/99_otu_taxonomy.txt'
# programs
CONF['FAMAS'] = "/mnt/software/stow/famas-0.0.7/bin/famas"
CONF['PREFILTER'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "S16_V2_prefilter.py"))
CONF['PRIMER_TRIMMER'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "primer_trimmer.py"))
CONF['CLASSIFY_HITS'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "classify_hits.py"))
CONF['IDENT_TO_BAM'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "ident_to_bam.py"))
CONF['GRAPHMAP'] = '/mnt/software/stow/graphmap-0.2.2-dev-604a386/bin/graphmap'
CONF['BWA'] = '/mnt/software/stow/bwa-0.7.12/bin/bwa'
#CONF['BLASTN'] = '/mnt/software/stow/ncbi-blast-2.2.28+/bin/blastn'
CONF['CONVERT_TABLE'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "convert_table.py"))
CONF['DEBUG'] = False


def main():
    """The main function
    """

    for f in CONF.keys():
        if f in ['SSU_DB', 'SPIKEIN-NAME', 'DEBUG']:
            continue
        if not os.path.exists(CONF[f]):
            LOG.fatal("Missing file: {}".format(CONF[f]))
            sys.exit(1)
    assert os.path.exists(SNAKEMAKE_TEMPLATE)


    parser = argparse.ArgumentParser(description='16S pipeline: version 2')
    parser.add_argument('-1', "--fq1", required=True, nargs="+",
                        help="Paired-end FastQ file #1 (gzip supported). Multiple (split) input files allowed")
    parser.add_argument('-2', "--fq2", required=True, nargs="+",
                        help="Paired-end FastQ file #2 (gzip supported). Multiple (split) input files allowed")
    parser.add_argument('-o', "--outdir", required=True,
                        help='Output directory (may not exist, unless using --continue)')
    parser.add_argument('-i', '--ins-len', type=int, required=True,
                        help='Mean insert size')
    parser.add_argument('-m', '--max-read-len', type=int, required=True,
                        help='Max. read length i.e. sequencing library read size')

    default = 40
    parser.add_argument('--ins-stdev', type=int, default=default,
                        help='Insert size standard deviation (default {})'.format(default))
    default = 8
    parser.add_argument('-c', '--num-cores', type=int, default=default,
                        help='Number of cores to use (default = {})'.format(default))
    parser.add_argument('--verbose', action="store_true",
                        help='Be verbose')
    parser.add_argument('--debug', action="store_true",
                        help='Enable debugging output')
    parser.add_argument('--continue', action="store_true", dest="continue_run",
                        help='Continue interrupted run')
    parser.add_argument('--no-run', action="store_true",
                        help="Prepare output directory and files but don't actually submit job")
    # FIXME implement the following
    #parser.add_argument('--no-amplicon',
    #                    help='Use default EMIRGE, not the amplicon version')
    #parser.add_argument('--max-num-reads', default=500000,
    #                    help='Downsample to this number of reads')
    #parser.add_argument('--no-fastq-sort', action="store_true",
    #                    help="Don't sort FastQ but use given order in forward and reverse list of reads (careful now!)")
    #parser.add_argument('--phred64', action="store_true",
    #                    help='Assume Illumina 1.3-1.7 quality encoding')
    args = parser.parse_args()


    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    if args.continue_run:
        raise NotImplementedError()

    if os.path.exists(args.outdir):
        if not args.continue_run:
            LOG.fatal("Output directory must not exist: {}".format(args.outdir))
            sys.exit(1)
    else:
        if args.continue_run:
            LOG.fatal("Can't continue job that wasn't started yet")
            sys.exit(1)

    if args.fq1 == args.fq2:
        LOG.fatal("Paired-End FastQ files have identical names")
        sys.exit(1)
    for fq1, fq2 in zip_longest(args.fq1, args.fq2):
        # only i|zip_longest uses None if one is missing
        if fq1 is None or fq2 is None:
            LOG.fatal("Unequal number of FastQ files")
            sys.exit(1)
        # enforce gzipped fastq
        # check they all exist
        for f in [fq1, fq2]:
            if not os.path.exists(f):
                LOG.fatal("FastQ file {} does not exist".format(f))
                sys.exit(1)
            elif not f.endswith(".gz"):
                LOG.fatal("Non-gzipped FastQ files not supported")
                sys.exit(1)

    os.mkdir(args.outdir)
    os.mkdir(os.path.join(args.outdir, "results"))

    samples = []
    fqs1 = [os.path.abspath(f) for f in args.fq1]
    fqs2 = [os.path.abspath(f) for f in args.fq2]
    #if not args.no_fastq_sort:
    fqs1 = sorted(fqs1)
    fqs2 = sorted(fqs2)
    for i, (fq1, fq2) in enumerate(zip(fqs1, fqs2)):
        os.symlink(os.path.abspath(fq1), os.path.join(args.outdir, "{}_R1.fastq.gz".format(i+1)))
        os.symlink(os.path.abspath(fq2), os.path.join(args.outdir, "{}_R2.fastq.gz".format(i+1)))
        samples.append("{}_".format(i+1))


    config_file = os.path.join(args.outdir, CONFIG_FILE)
    conf = copy.deepcopy(CONF)
    conf['SAMPLES'] = samples
    conf['USE_AMPLICON'] = True
    # ints needs to saved as string in json to be loaded proprerly from within snakemake
    conf['INS_SIZE'] = str(args.ins_len)
    conf['INS_STDEV'] = str(args.ins_stdev)
    conf['MAX_READ_LEN'] = str(args.max_read_len)
    with open(config_file, 'w') as fh:
        json.dump(conf, fh, indent=4)

    shutil.copyfile(SNAKEMAKE_TEMPLATE, 
                    os.path.join(args.outdir, SNAKEMAKE_FILE))

    snakemake_cluster_wrapper = os.path.join(args.outdir, SNAKEMAKE_CLUSTER_WRAPPER)
    mail_option = "-m bes -M {}@gis.a-star.edu.sg".format(getpass.getuser())
    with open(snakemake_cluster_wrapper, 'w') as fh:
        fh.write('# snakemake requires python3\n')
        fh.write('source activate py3k;\n')
        fh.write('cd {};\n'.format(os.path.abspath(args.outdir)))
        fh.write('# qsub for snakemake itself\n')
        fh.write('qsub="qsub -pe OpenMP 1 -l mem_free=1G -l h_rt=48:00:00 {} -j y -V -b y -cwd";\n'.format(mail_option))
        fh.write('# -j in cluster mode is the maximum number of spawned jobs\n')
        fh.write('$qsub -N snakemake -o snakemake.qsub.log \'snakemake -j 8 -c "qsub -pe OpenMP {{threads}} -l mem_free=12G -l h_rt=24:00:00 -j y -V -b y -cwd" -s {} --configfile {}\';\n'.format(
                SNAKEMAKE_FILE, CONFIG_FILE))
        # FIXME max runtime should be defined per target in SNAKEMAKE_FILE

    cmd = ['bash', snakemake_cluster_wrapper]
    if args.no_run:
        print("Not actually submitting job")
        print("When ready run: {}".format(' '.join(cmd)))
    else:
        print("Running: {}".format(cmd))
        subprocess.check_output(cmd)
                

if __name__ == '__main__':
    main()

