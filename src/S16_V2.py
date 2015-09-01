#!/usr/bin/env python
"""This is a remake of Vinutha's original 16S pipeline, which was
designed for Illumina shotgun Sequencing of 16S rRNA Amplicon
Sequences (Ong et al., 2013, PMID 23579286). At its core it's running
EMIRGE (Miller et al., 2011, PMID 21595876) for reconstructing the
sequences and BLAST + Greengenes for classification

This pipelines works as follows:
- If necesssary, concat split fastq files, remove Q2 bases from 3', check paired end read name consistency
- Run FastQC
- Run EMIRGE or EMIRGE amplicon
- Trim primers
- Infer taxonomy from best hit against Greengenes
- Compile tables

TODO:
- BLAST vs Graphmap
- Use a pipeline framework like bpipe and serves
  only as a stop-gap.
- Add decont step
- Add downsampling option
"""


import os
import sys
import argparse
from itertools import groupby
import logging
import subprocess
import io
import gzip
import datetime
import shutil
import glob
from collections import OrderedDict
from collections import namedtuple
import csv

try:
    # python 3
    from itertools import zip_longest
except ImportError:
    # python 2
    from itertools import izip_longest as zip_longest


#import pysam
#
#PYSAM_VERSION = [int(x) for x in pysam.__version__.split(".")]
#assert PYSAM_VERSION >= [0, 8, 0]

IDENT_TAG = 'Xi'


__author__ = "Andreas Wilm"
__version__ = "2.0.0a"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "The MIT License (MIT)"



# global logger
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
    format='%(levelname)s [%(asctime)s]: %(message)s')



CONF = dict()
CONF['famas'] = "/mnt/software/stow/famas-0.0.7/bin/famas"
CONF['fastqc'] = "/mnt/software/stow/FastQC-0.10.1/bin/fastqc"
#CONF['emirge'] = "/mnt/software/stow/emirge-v0.60-15-g0ddae1c-anaconda/bin/emirge.py"
CONF['emirge'] = "/mnt/software/stow/emirge-v0.60-15-g0ddae1c-wilma/bin/emirge.py"
CONF['emirge-amplicon'] = os.path.join(os.path.dirname(CONF['emirge']), 'emirge_amplicon.py')
CONF['emirge-rename'] = os.path.join(os.path.dirname(CONF['emirge']),'emirge_rename_fasta.py')
#CONF['emirge-fasta'] = '/mnt/software/unstowable/16S_pipeline/SSU_candidate_db.fasta'
CONF['emirge-fasta'] = '/mnt/genomeDB/misc/softwareDB/emirge/SSU_candidate_db.fasta'
CONF['emirge-fasta-bowtie'] = CONF['emirge-fasta'].replace('.fasta', '')
CONF['primer-trimmer'] = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "primer_trimmer.py"))
#CONF['greengenes'] = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5.fasta'
#CONF['greengenes'] = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus/rep_set/99_otus.fasta'
CONF['greengenes'] = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "../greengenes/99_otus.fasta"))
CONF['greengenes-taxonomy'] = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus/taxonomy/99_otu_taxonomy.txt'
CONF['ident_to_bam'] = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "ident_to_bam.py"))
CONF['blastn'] = '/mnt/software/stow/ncbi-blast-2.2.28+/bin/blastn'
CONF['graphmap'] = '/mnt/software/stow/graphmap-0.2.2-dev-604a386/bin/graphmap'
CONF['best-blast-hit'] = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "best_hit_per_query_from_multi_blast.sh"))



BestHit = namedtuple('BestHit', ['seqid', 'pwid', 'tax'])


def timestamp():
    return datetime.datetime.now().isoformat()


def touch_completed(f):
    """FIXME
    """
    assert not os.path.exists(f)
    with open(f, 'w') as fh:
        fh.write(timestamp())


def fasta_iter(fasta_stream):
    """Brend Pedersen:  https://www.biostars.org/p/710/

    Given a fasta file. yield tuples of header, sequence
    """

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fasta_stream, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq



def get_fastq_nreads_maxlen(fastq):
    """FIXME
    """

    if fastq.endswith(".gz"):
        # gzip speedup with io.BufferedReader
        # see http://aripollak.com/pythongzipbenchmarks/
        # and http://www.reddit.com/r/Python/comments/2olhrf/fast_gzip_in_python/
        fh = io.BufferedReader(gzip.open(fastq))
    else:
        fh = open(fastq)

    max_read_len = num_lines = 0
    for num_lines, line in enumerate(fh):
        if num_lines % 4 == 0:
            assert line[0] == '@', (fastq, num_lines, line)
        elif num_lines % 4 == 1:
            l = len(line.strip())
            if l > max_read_len:
                max_read_len = l
    fh.close()

    assert (num_lines+1) % 4 == 0
    return (num_lines+1)/4, max_read_len


class JobFailedException(Exception):
    pass


def run_cmd(cmd, log_stdout=sys.stdout, log_stderr=sys.stderr):
    """Wrapper to Popen. Will raise an error if command execution
    failed with non-zero exit code
    """

    msg = "Running {}".format(' '.join(cmd))
    log_stdout.write(msg)
    LOG.info(msg)

    p = subprocess.Popen(cmd,# shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()

    for line in stdout.splitlines():
        log_stdout.write("stdout: %s\n" % line)
    for line in stderr.splitlines():
        log_stderr.write("stderr: %s\n" % line)

    # FIXME will only print stdout and stderr once done and also keep
    # all the output in memory right? Better to use shell and redirect
    # instead?

    if p.returncode != 0:
        raise JobFailedException(
            "Following command failed with exit status"
            " %s: '%s'. stderr was '%s'" % (
                p.returncode, ' '.join(cmd), stderr))

    return (stdout.splitlines(), stderr.splitlines())


def get_max_readlen_from_pair(fq1, fq2):
    """FIXME
    """
    num_reads_1, max_read_len_1 = get_fastq_nreads_maxlen(fq1)
    num_reads_2, max_read_len_2 = get_fastq_nreads_maxlen(fq2)
    assert num_reads_1 == num_reads_2
    max_read_len = max(max_read_len_1, max_read_len_2)
    return max_read_len



def read_blast_csv(blast_csv):
    """
    Expected fields in file:
    # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    """

    with open(blast_csv, 'rb') as fh:
        csvreader = csv.reader(fh, delimiter=" ")
        for row in csvreader:
            if row[0].startswith("#"):
                continue
            assert len(row)==6
    raise NotImplementedError


def parse_best_blast_hit(best_hit_file):
    """using our own format
    """

    query_to_hit = dict()
    with open(best_hit_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            query, hit, pwid = line.split(" ")[:3]
            pwid = abs(float(pwid))
            assert not query_to_hit.has_key(query)
            query_to_hit[query] = BestHit(seqid=hit, pwid=pwid, tax=None)
    return query_to_hit


def parse_graphmap(bam, ident_tag=IDENT_TAG):
    """FIXME
    """

    query_to_hit = dict()
    #cmd = ['samtools', 'view', bam]
    #proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)#, stderr=subprocess.STDOUT)
    #works in python 3.0+
    #for line in proc.stdout:
    #for line in iter(proc.stdout.readline,''):
    #    ls = line.rstrip().split("\t")
    #    query_to_hit[ls[0]] = ls[2]
    #return query_to_hit
    samfh = pysam.Samfile(bam)
    for r in samfh:
        pwid = round(r.get_tag(ident_tag), 1)
        ref = samfh.getrname(r.tid)
        query_to_hit[r.query_name] = BestHit(seqid=ref, pwid=pwid, tax=None)
    samfh.close()
    return query_to_hit


def read_greengenes_taxonomy(gg_tax_file):
    """returns greengenes taxonomy with OTU as key (string) and tax
    implemented as OrderedDict() as value

    example entry:
    4470549 k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__; s__
    """

    gg_tax = dict()
    with open(gg_tax_file) as fh:
        for line in fh:
            otu, tax_str = line.rstrip().split("\t")
            assert otu.isdigit()
            # remove semicolon and split as key__value
            assert not gg_tax.has_key(otu)
            gg_tax[otu] = OrderedDict([x.split("__") for x in tax_str.split("; ")])
            assert len(gg_tax[otu]) == 7
    return gg_tax


def main():
    """The main function
    """

    for f in CONF.keys():
        if f in ['emirge-fasta-bowtie']:
            continue
        if not os.path.exists(CONF[f]):
            LOG.fatal("Missing file: {}".format(CONF[f]))
            sys.exit(1)

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
                        help='Max. read length')

    default = 40
    parser.add_argument('--ins-stdev', type=int, default=default,
                        help='Insert size standard deviation (default {})'.format(default))
    default = 8
    parser.add_argument('-c', '--num-cores', type=int, default=default,
                        help='Number of cores to use (default = {})'.format(default))
    parser.add_argument('--no-amplicon',
                        help='Use default EMIRGE, not the amplicon version')
    parser.add_argument('--verbose', action="store_true",
                        help='Be verbose')
    parser.add_argument('--no-fastq-sort', action="store_true",
                        help="Don't sort FastQ but use given order in forward and reverse list of reads (careful now!)")
    parser.add_argument('--phred64', action="store_true",
                        help='Assume Illumina 1.3-1.7 quality encoding')
    parser.add_argument('--debug', action="store_true",
                        help='Enable debugging output')
    parser.add_argument('--continue', action="store_true", dest="continue_run",
                        help='Continue interrupted run')
    #parser.add_argument('--max-num-reads', default=500000,
    #                    help='Downsample to this number of reads')
    args = parser.parse_args()


    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    if os.path.exists(args.outdir):
        if not args.continue_run:
            LOG.fatal("Output directory must not exist: {}".format(args.outdir))
            sys.exit(1)
    else:
        if args.continue_run:
            LOG.fatal("Can't continue job that wasn't started yet")
            sys.exit(1)
        os.mkdir(args.outdir)

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


    # relative to args.outdir
    emirge_outdir = "emirge"
    emirge_out_fa = "emirge_out.fa"
    emirge_primer_trimmed_fa = emirge_out_fa.replace(".fa", "primer_trimmed.fa")
    gg_out =  "blast-greengenes.csv"
    gg_best = "blast-greengenes-best.txt"
    raw_table = "raw-table.csv"


    snakemake_file = os.path.join(args.outdir, "snake.make")
    # FIXME what if exists?
    snakemake_fh = open(snakemake_file, 'w')
    snakemake_fh.write('shell.prefix("set -o pipefail; ")\n')
    snakemake_fh.write('\n')


    samples = []
    fqs1 = [os.path.abspath(f) for f in args.fq1]
    fqs2 = [os.path.abspath(f) for f in args.fq2]
    if not args.no_fastq_sort:
        fqs1 = sorted(fqs1)
        fqs2 = sorted(fqs2)
    for i, (fq1, fq2) in enumerate(zip(fqs1, fqs2)):
        os.symlink(os.path.abspath(fq1), os.path.join(args.outdir, "{}_R1.fastq.gz".format(i+1)))
        os.symlink(os.path.abspath(fq2), os.path.join(args.outdir, "{}_R2.fastq.gz".format(i+1)))
        samples.append("{}_".format(i+1))


    snakemake_fh.write("SAMPLES = ['{}']\n".format("' ,'".join(samples)))
    snakemake_fh.write("\n")

    snakemake_fh.write("# must be first rule\n")
    snakemake_fh.write("\nrule final:\n")
    #snakemake_fh.write("  input: 'raw_table.csv'\n")
    snakemake_fh.write("  input: '{}'\n".format(gg_out))
    snakemake_fh.write("  message: 'This is the end. My only friend, the end'\n")
    snakemake_fh.write("\n")


    snakemake_fh.write("\nrule filter_fastq:\n")
    snakemake_fh.write("  input: fq1='{sample}R1.fastq.gz', fq2='{sample}R2.fastq.gz'\n")
    snakemake_fh.write("  output: fq1=temp('{sample}R1.flt.fastq.gz'), fq2=temp('{sample}R2.flt.fastq.gz')\n")
    if args.phred64:
        extra_arg = "--phred64"
    else:
        extra_arg = ""
    #snakemake_fh.write("# zcat into anonymous named pipes to stream into famas for filtering and merging\n")
    #snakemake_fh.write("  shell: '{} -i <(zcat {{input.fqs1}}) -j <(zcat {{input.fqs2}}) -o {{output.fq1}} -p {{output.fq2}} -q 3 -l 60 {} --append'\n".format(CONF['famas'], extra_arg))
    snakemake_fh.write("  shell: '{} -i {{input.fq1}} -j {{input.fq2}} -o {{output.fq1}} -p {{output.fq2}} -q 3 -l 60 {} --append'\n".format(CONF['famas'], extra_arg))
    snakemake_fh.write("\n")



    snakemake_fh.write("\nrule concat_fastq:\n")
    flt_fq1 = "R1.flt.fastq.gz"
    flt_fq2 = "R2.flt.fastq"# 2nd emirge input file has to be nonzipped
    snakemake_fh.write("  input: fqs1=expand('{sample}R1.flt.fastq.gz', sample=SAMPLES), fqs2=expand('{sample}R2.flt.fastq.gz', sample=SAMPLES)\n")
    snakemake_fh.write("  output: fq1=temp('{}'), fq2=temp('{}')\n".format(flt_fq1, flt_fq2))
    #snakemake_fh.write("  message: ''\n")
    snakemake_fh.write("  shell: 'zcat {input.fqs1} | gzip > {output.fq1} && zcat {input.fqs2} > {output.fq2}'\n")
    snakemake_fh.write("\n")

    # emirge
    #
    # this assumes we can't continue aborted emirge runs
    #if os.path.exists(emirge_outdir):
    #    shutil.rmtree(emirge_outdir)
    snakemake_fh.write("\nrule emirge:\n")
    snakemake_fh.write("  input: fq1=rules.concat_fastq.output.fq1, fq2=rules.concat_fastq.output.fq2\n")
    snakemake_fh.write("  threads: {}\n".format(args.num_cores))
    snakemake_fh.write("  output: 'emirge_succeeded'\n")# fake
    if args.no_amplicon:
        cmd = [CONF['emirge']]
    else:
        cmd = [CONF['emirge-amplicon']]
    cmd.extend(["-l", str(args.max_read_len), "-i", str(args.ins_len), "-s", str(args.ins_stdev)])
    if not args.phred64:
        cmd.extend(["--phred33"])
    sys.stderr.write("DEBUG only one iteration\n"); cmd.extend(["-n 1"])
    cmd.extend(["-a", str(args.num_cores), emirge_outdir])
    cmd.extend(["-f", CONF['emirge-fasta'], "-b", CONF['emirge-fasta-bowtie']])
    snakemake_fh.write("  shell:\n")
    snakemake_fh.write("    'test -d {} && rm -rf {}; {} -1 {{input.fq1}} -2 {{input.fq2}} && touch {{output}}'\n".format(
                emirge_outdir, emirge_outdir, ' '.join(cmd)))
    snakemake_fh.write("\n")


    snakemake_fh.write("\nrule emirge_rename:\n")
    snakemake_fh.write("  input: rules.emirge.output\n")
    snakemake_fh.write("  output: '{}'\n".format(emirge_out_fa))
    snakemake_fh.write('  run:\n')
    snakemake_fh.write('    import glob\n')
    snakemake_fh.write('    import subprocess\n')
    snakemake_fh.write('    last_iter = sorted(glob.glob(os.path.join("{}", "iter.*")))[-1]\n'.format(emirge_outdir))
    snakemake_fh.write('    cmd = ["{}", last_iter]\n'.format(CONF["emirge-rename"]));
    snakemake_fh.write('    with open(output, "w") as fh:\n')
    snakemake_fh.write('        subprocess.call(cmd, stdout=fh)\n')
    snakemake_fh.write("\n")


    snakemake_fh.write("\nrule emirge_trim_primer:\n")
    snakemake_fh.write("  input: rules.emirge_rename.output\n")
    snakemake_fh.write("  output: '{}'\n".format(emirge_primer_trimmed_fa))
    cmd = [CONF['primer-trimmer'], "-i", "{input}", "-o", "{output}",
           "--minlen", "200", "--maxlen", "1400"]
    snakemake_fh.write("  shell: '{}'\n".format(' '.join(cmd)))
    snakemake_fh.write("\n")


    snakemake_fh.write("\nrule emirge_vs_gg:\n")
    snakemake_fh.write("  input: rules.emirge_trim_primer.output, ref='{}'\n".format(CONF['greengenes']))
    snakemake_fh.write("  output: '{}'\n".format(gg_out))
    snakemake_fh.write("  threads: {}\n".format(args.num_cores))
    cmd = "{} -x illumina -t {{threads}} -r {{input.ref}} -d {{input}} | {} - {{output}} {{input.ref}}".format(
        CONF['graphmap'], CONF['ident_to_bam'])
    snakemake_fh.write("  shell: '{}'\n".format(cmd))
    sys.exit(1)



    sys.exit(0)

    LOG.fatal("Unfinished implementation")
    sys.exit(1)

"""
  out: class/abundance
  for all and
  lowid==novel, if id below p-80,f-90,g-95,s-97
  lowabd==rare, if abundance threshold 0.00005
"""


if __name__ == '__main__':
    LOG.warn("Test fastq split")
    LOG.warn("Make UGE ready")
    LOG.warn("Dump config")
    LOG.warn("Create log")
    main()

