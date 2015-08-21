#!/usr/bin/env python
"""This is a remake of Vinutha's original 16S pipeline, which was
designed for Illumina shotgun Sequencing of 16S rRNA Amplicon
Sequences (Ong et al., 2013, PMID 23579286). At its core it's running
EMIRGE (Miller et al., 2011, PMID 21595876) for reconstructing the
sequences and BLAST + Greengenes for classification

This pipelines works as follows:
- Concat files if split, remove Q2 reads from 3' (as in original
  EMIRGE paper), check paired end read name consistency
  (NotImplemented: optionally downsample)
- Run FastQC
- Run EMIRGE or EMIRGE amplicon
- Trim primers
- Blast against Greengenes
- Infer best Blasthit + abundance and classification

The previous pipeline worked as follows:
1. Check order of PE reads
2. Downsample
3. FastQC
4. Filter 3' Q2
5. Emirge
6. Remame 
7. Trim primers
8. Blast against Greengenes
9. Infer best Blasthit + abundance and classification
10. rename/reorg.cleanup


TODO: 

- This should really use a pipeline framework like bpipe and serves
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
import csv

try:
    # python 3
    from itertools import zip_longest
except ImportError:
    # python 2
    from itertools import izip_longest as zip_longest

#--- project specific imports
#
#/

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
CONF['emirge'] = "/mnt/software/stow/emirge-v0.60-15-g0ddae1c-anaconda/bin/emirge.py"
CONF['emirge-amplicon'] = "/mnt/software/stow/emirge-v0.60-15-g0ddae1c-anaconda/bin/emirge_amplicon.py"
CONF['emirge-rename'] = '/mnt/software/bin/emirge_rename_fasta.py'
#CONF['emirge-fasta'] = '/mnt/software/unstowable/16S_pipeline/SSU_candidate_db.fasta'
CONF['emirge-fasta'] = '/mnt/genomeDB/misc/softwareDB/emirge/SSU_candidate_db.fasta'
CONF['emirge-fasta-bowtie'] = CONF['emirge-fasta'].replace('.fasta', '')
CONF['primer-trimmer'] = os.path.join(os.path.dirname(sys.argv[0]), "primer_trimmer.py")
#CONF['greengenes'] = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5.fasta'
CONF['greengenes'] = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus/rep_set/99_otus.fasta'
CONF['greengenes-taxonomy'] = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus/taxonomy/99_otu_taxonomy.txt'
CONF['blastn'] = '/mnt/software/stow/ncbi-blast-2.2.28+/bin/blastn'
CONF['best-blast-hit'] = os.path.join(os.path.dirname(sys.argv[0]), "best_hit_per_query_from_multi_blast.sh")




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
            query, hit = line.split(" ")[:2]
            assert not query_to_hit.has_key(query)
            query_to_hit[query] = hit
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
                        help="Don't sort FastQ but use given order in forward and reverse list of reads")
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


    log_file = os.path.join(args.outdir, "{}.log".format(
            os.path.basename(sys.argv[0])))
    log_fh = open(log_file, 'a')
    # log_fh = sys.stdout
    log_fh.write("Starting {}: {}\n".format(timestamp(), sys.argv))


    # FIXME move below to CONF
    concat_filtered_fq1 = os.path.join(args.outdir, "concat_filtered_1.fastq.gz")
    concat_filtered_fq2 = os.path.join(args.outdir, "concat_filtered_2.fastq.gz")
    emirge_outdir = os.path.join(args.outdir, "emirge")
    emirge_out_fa = os.path.join(args.outdir, "emirge_out.fa")
    emirge_abundance_out = os.path.join(args.outdir, "abundance.txt")
    emirge_primer_trimmed_fa = emirge_out_fa.replace(".fa", "primer_trimmed.fa")
    blast_gg_out = os.path.join(args.outdir, "blast-greengenes.csv")
    blast_gg_best = os.path.join(args.outdir, "blast-greengenes-best.txt")


    # FIXME make below functions or create makefile


    # concat fastq, filter, merge with famas
    #
    stage = "concat-filter-and-merge-fasta"
    compl_file = os.path.join(args.outdir, ".{}.completed".format(stage))
    if os.path.exists(compl_file):
        LOG.info("Skipping stage {}".format(stage))

    else:
        if args.no_fastq_sort:
            fqs1 = sorted(args.fq1)
            fqs2 = sorted(args.fq2)
        else:
            fqs1 = args.fq1
            fqs2 = args.fq2
    
        for fq1, fq2 in zip_longest(fqs1, fqs2):
            cmd = [CONF['famas'], "-i", fq1, "-j", fq2]
            if args.phred64:
                cmd.extend(["--phred64"])
            cmd.extend(["-o", concat_filtered_fq1, "-p", concat_filtered_fq2])
            cmd.extend(["-q", "3", "-l", "60", "--append"])
            (stdout, stderr) = run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)

        touch_completed(compl_file)

            
    # plain fastqc
    #
    stage = "fastqc"
    compl_file = os.path.join(args.outdir, ".{}.completed".format(stage))
    if os.path.exists(compl_file):
        LOG.info("Skipping stage {}".format(stage))

    else:
        for f in [concat_filtered_fq1, concat_filtered_fq2]:
            cmd = [CONF['fastqc'], "--nogroup", "--quiet",
                   "--threads", str(args.num_cores), f]
            (stdout, stderr) = run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)

        touch_completed(compl_file)

    
    # emirge
    #
    stage = "emirge"
    compl_file = os.path.join(args.outdir, ".{}.completed".format(stage))
    if os.path.exists(compl_file):
        LOG.info("Skipping stage {}".format(stage))

    else:
        # this assumes we can't continue aborted emirge runs
        if os.path.exists(emirge_outdir):
            shutil.rmtree(emirge_outdir)

        max_read_len = get_max_readlen_from_pair(
            concat_filtered_fq1, concat_filtered_fq2)

        # 2nd fastq may not be gzipped
        concat_filtered_fq2_unzipped = concat_filtered_fq2.replace(".gz", "")
        with open(concat_filtered_fq2_unzipped, 'w') as fh:
            subprocess.call(["gzip", "-dc", concat_filtered_fq2], stdout=fh)

        if args.no_amplicon:
            emirge = CONF['emirge']        
        else:
            emirge = CONF['emirge-amplicon']

        cmd = [emirge, "-1", concat_filtered_fq1, "-2", concat_filtered_fq2_unzipped]
        cmd.extend(["-f", CONF['emirge-fasta'], "-b", CONF['emirge-fasta-bowtie']])
        cmd.extend(["-l", str(max_read_len), "-i", str(args.ins_len), "-s", str(args.ins_stdev)])
        if not args.phred64:
            cmd.extend(["--phred33"])
        cmd.extend(["-a", str(args.num_cores), emirge_outdir])
        LOG.warn("One iteration only")
        cmd.extend(["-n", "1"])
        (stdout, stderr) = run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)
    
        # remove temporarily unzipped 2nd fastq file
        os.unlink(concat_filtered_fq2_unzipped)

        # clean BAMs from emirge
        for root, directories, filenames in os.walk(emirge_outdir):
            for filename in filenames:
                if filename.endswith(".bam"):
                    os.unlink(os.path.join(root, filename))

        touch_completed(compl_file)


    # rename emirge output from last iteration 
    #
    stage = "emirge-rename"
    compl_file = os.path.join(args.outdir, ".{}.completed".format(stage))
    if os.path.exists(compl_file):
        LOG.info("Skipping stage {}".format(stage))

    else:
        # FIXME iteration hardcoded
        last_iter = sorted(glob.glob(os.path.join(emirge_outdir, "iter.*")))[-1]
        LOG.info("Using EMIRGE fasta from {}".format(last_iter))
        cmd = [CONF['emirge-rename'], last_iter]
        with open(emirge_out_fa, 'w') as fh:
            subprocess.call(cmd, stdout=fh)

        touch_completed(compl_file)


    # trim primers
    #
    # preserves abundances, but might kick out some sequences due to
    # length restrictions
    #
    stage = "primer-trimming"
    compl_file = os.path.join(args.outdir, ".{}.completed".format(stage))
    if os.path.exists(compl_file):
        LOG.info("Skipping stage {}".format(stage))

    else:
        # FIXME hardcoded params
        cmd = [CONF['primer-trimmer'], "-i", emirge_out_fa, "-o", emirge_primer_trimmed_fa,
               "--minlen", "200", "--maxlen", "1400"]
        (stdout, stderr) = run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)

        touch_completed(compl_file)


    # FIXME: we could in theory also use the EMIRGE/ARB markup for direct taxonomy assignment
    # Never benchmarked though and not sure how assignment happens after split within EMIRGE.
    # Other option is to use the GG database directly (but split problem applies here as well)

    # blast against greengenes
    #
    # FIXME this is slow. why not use graphmap. also assigns simply
    # best hit, so no need to infer it like for blast
    # 
    stage = "blast-against-greengenes"
    compl_file = os.path.join(args.outdir, ".{}.completed".format(stage))
    if os.path.exists(compl_file):
        LOG.info("Skipping stage {}".format(stage))

    else:
        cmd = [CONF['blastn'], "-db", CONF['greengenes'], "-query", emirge_primer_trimmed_fa,
               "-outfmt", "7", "-out", blast_gg_out, "-num_threads", str(args.num_cores)]
        (stdout, stderr) = run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)

        touch_completed(compl_file)

        
    # best blast hit
    #
    # FIXME uses unreadable bash script
    # can't we just take the first hit?
    #
    stage = "best-blast-hit"
    compl_file = os.path.join(args.outdir, ".{}.completed".format(stage))
    if os.path.exists(compl_file):
        LOG.info("Skipping stage {}".format(stage))

    else:
        cmd = [CONF['best-blast-hit'], blast_gg_out]
        with open(blast_gg_best, 'w') as fh:
            subprocess.call(cmd, stdout=fh)

        touch_completed(compl_file)


    #stage  = "abundance-inference"
    #compl_file = os.path.join(args.outdir, ".{}.completed".format(stage))
    #if os.path.exists(compl_file):
    #    LOG.info("Skipping stage {}".format(stage))
    #
    #else:
    #    assert not os.path.exists(emirge_abundance_out)
    #    with open(emirge_out_fa) as fh_in, open(emirge_abundance_out, 'w') as fh_out:
    #        for line in fh_in:
    #            if not line.startswith(">"):
    #                continue
    #            ls = line.strip().split()
    #            fh_out.write("{}\t{}\n".format(ls[0], ls[3].replace("NormPrior=", "")))
    #            
    #    touch_completed(compl_file)

    
    # FIXME
    # graphmap  -x  illumina -t 2 -r /mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus/rep_set/99_otus.fasta  -d schmock/emirge_outprimer_trimmed.fa | samtools view -bS - | samtools sort - schmock/emirge_outprimer_trimmed_gg_13_5_99_otus
    # all start with 50bp deletion? samtools view emirge_outprimer_trimmed_gg_13_5_99_otus.bam
        
    LOG.info("Reading taxonomy from {}".format(CONF['greengenes-taxonomy']))
    gg_tax = read_greengenes_taxonomy(CONF['greengenes-taxonomy'])
    query_to_hit = parse_best_blast_hit(blast_gg_best)
    LOG.info("Assigning taxonomy".format())
    with open(emirge_primer_trimmed_fa) as fh:
        for (s_id, s_seq) in fasta_iter(fh):
            # id's look as follows:
            # >47|AJ704791.1.1593 Prior=0.045024 Length=744 NormPrior=0.045018
            s_id_split = s_id.split()
            abundance = float(s_id_split[-1].replace("NormPrior=", ""))
            s_id = s_id_split[0]
            best_hit = query_to_hit.get(s_id, None)
            # since blast was run against gg otus we can read from gg tax
            if best_hit is not None:
                tax = gg_tax.get(best_hit, None)
            else:
                tax = None
            print "{}\t{}\t{}\t{}".format(s_id, best_hit, abundance, tax)


    if log_fh != sys.stdout:
        log_fh.close()
    LOG.fatal("Unfinished implementation")
    sys.exit(1)
"""
    
    OUT_CLASSIFN = OUTPUT + "/classification.txt"
    CMD_GET_CLASSIFICATION = "python %s %s %s %s > %s"% (GET_CLASSIFICATION, OUT_BEST_HIT, CURRENT_GG, OUT_ABUN, OUT_CLASSIFN)
    run_command(CMD_GET_CLASSIFICATION)
    
    #OUT_LOW_HITS = OUTPUT + "/Low_identity_hits.txt"
    #OUT_FINAL_TABLE = OUTPUT + "/Final_classification_table.txt"
    print "Step 10: Creating final tables"
    
    #CMD_GET_TABLE = "python %s %s %s %s"% (GET_TABLE, OUT_CLASSIFN, OUTPUT, OUTPUT + "/Results")
    cmd_mv = "mv %s %s"% (OUT_CLASSIFN, OUTPUT + "/Analysis/")
    run_command(cmd_mv)
    CMD_GET_TABLE = "python %s %s %s"% (GET_TABLE, OUTPUT, OUTPUT + "/Results")
    #print CMD_GET_TABLE
    run_command(CMD_GET_TABLE)
    print "Done"
    """
    

if __name__ == '__main__':
    LOG.warn("Test fastq split")
    LOG.warn("Make UGE ready")
    LOG.warn("Dump config")
    LOG.warn("Create log")
    main()
    
