#!/mnt/software/unstowable/anaconda/envs/pysam-0.8.2/bin/python
"""16S pipeline prefilter. Determines ratio between spike-in,
unspecific and 16S product based on given mapping against
16S(+spikein) database
"""


import sys
import os
import argparse
import io
import gzip

import pysam
PYSAM_VERSION = [int(x) for x in pysam.__version__.split(".")]
assert PYSAM_VERSION >= [0, 8, 0]


# DNA base complements
COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}


def write_read(fastq, read):
    """
    Write read to open FASTQ file.

    From https://github.com/martijnvermaat/bio-playground/blob/master/bam-to-fastq/bam_to_fastq.py
    """
    info = {'index': int(not read.is_read1) + 1,
            'name':  read.qname}
    if read.is_reverse:
        info.update({'quality':  read.qual[::-1],
                     'sequence': reverse_complement(read.seq)})
    else:
        info.update({'quality':  read.qual,
                     'sequence': read.seq})
    fastq.write('@{name}/{index}\n{sequence}\n+\n{quality}\n'.format(**info))


def reverse_complement(sequence):
    """
    Return reverse complement of DNA sequence.

    From https://github.com/martijnvermaat/bio-playground/blob/master/bam-to-fastq/bam_to_fastq.py
    """
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])


def pairs_from_name_srtd_bam(samfh):
    """Get pairs from name sorted BAM
    """
    for r1 in samfh:
        while r1.is_secondary:# or r1.is_supplementary:
            r1 = samfh.next()
        r2 = samfh.next()
        while r2.is_secondary:# or r2.is_supplementary:
            r2 = samfh.next()
        assert r1.query_name == r2.query_name, (r1.query_name, r2.query_name)
        yield (r1, r2)


def  main():
    """FIXME
    """

    parser = argparse.ArgumentParser(description='16S pipeline: version 2')
    parser.add_argument('-i', "--bam", required=True,
                        help="BAM input")
    parser.add_argument('-1', "--fq1", required=True,
                        help="Paired-end FastQ file #1 (gzip supported)")
    parser.add_argument('-2', "--fq2", required=True,
                        help="Paired-end FastQ file #2 (gzip supported)")
    parser.add_argument('-s', "--spikein", required=True,
                        help="Spike-in name, e.g. Plasmodium.knowlesi.profilin")

    args = parser.parse_args()

    #if args.verbose:
    #    LOG.setLevel(logging.INFO)
    #if args.debug:
    #    LOG.setLevel(logging.DEBUG)


    # open BAM/SAM(stdin)
    #
    # if stdin assuming SAM, not BAM
    if args.bam == "-":
        rmode = "r"
    else:
        assert os.path.exists(args.bam)
        rmode = "rb"
    samfh = pysam.AlignmentFile(args.bam, rmode)


    # open output fq
    #
    assert not os.path.exists(args.fq1)
    assert not os.path.exists(args.fq2)
    # gzip speedup with io.BufferedReader/io.BufferedWriter
    # see http://aripollak.com/pythongzipbenchmarks/
    # and http://www.reddit.com/r/Python/comments/2olhrf/fast_gzip_in_python/
    if args.fq1.endswith(".gz"):
        fq1fh = io.BufferedWriter(gzip.open(args.fq1, 'wb'))
    else:
        fq1fh = open(args.fq1, 'w')
    if args.fq2.endswith(".gz"):
        fq2fh = io.BufferedWriter(gzip.open(args.fq2, 'wb'))
    else:
        fq2fh = open(args.fq2, 'w')


    counts = dict()
    for (r1, r2) in pairs_from_name_srtd_bam(samfh):
        if r1.is_unmapped and r2.is_unmapped:
            k = 'Unspecific'
        elif args.spikein in [samfh.getrname(r1.tid), samfh.getrname(r2.tid)]:
            k = 'Spike-in'
        else:
            k = '16S'
            write_read(fq1fh, r1)
            write_read(fq2fh, r2)
        counts[k] = counts.get(k, 0) + 1

    counts_total = sum(counts.values())
    if counts_total == 0:
        sys.stderr.write("No read pairs encountered!\n")
        sys.exit(1)
    print "#type\trel.abundance"
    for k, v in counts.items():
        print "{}\t{}".format(k, v/float(counts_total))

    fq1fh.close()
    fq2fh.close()


if __name__ == "__main__":
    main()
