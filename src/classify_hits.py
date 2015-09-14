#!/mnt/software/unstowable/anaconda/envs/pysam-0.8.2/bin/python
"""Assign taxonomy based on greengenes hits
"""


import os
import sys
import argparse
from itertools import groupby
import logging
from collections import OrderedDict
from collections import namedtuple

import pysam
PYSAM_VERSION = [int(x) for x in pysam.__version__.split(".")]
# FIXME replace pysam dependency with https://github.com/mdshw5/simplesam


__author__ = "Andreas Wilm"
__version__ = "2.0.0a"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "The MIT License (MIT)"


# global logger
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
    format='%(levelname)s [%(asctime)s]: %(message)s')

IDENT_TAG = 'Xi'

BestHit = namedtuple('BestHit', ['seqid', 'pwid', 'tax'])


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


def parse_greengenes_taxonomy(gg_tax_file):
    """returns greengenes taxonomy with OTU as key (string) and tax
    implemented as OrderedDict() as value

    example entry:
    4470549 k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__; s__
    """

    gg_tax = dict()
    with open(gg_tax_file) as fh:
        for line in fh:
            otu, tax_str = line.rstrip().split("\t")
            assert otu.isdigit()# doesn't work with spike-in
            # remove semicolon and split as key__value
            assert not gg_tax.has_key(otu)
            gg_tax[otu] = OrderedDict([x.split("__") for x in tax_str.split("; ")])
            assert len(gg_tax[otu]) == 7

            # Taxa with square brackets "are suggested, but not verified, taxonomies."
            # https://groups.google.com/forum/#!topic/qiime-forum/WnhQiB5q9Hc
            # We treat them as verified here
            for k, v in gg_tax[otu].items():
                gg_tax[otu][k] = v.translate(None, "[]")
    return gg_tax


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


def parse_bam(bam, ident_tag=IDENT_TAG):
    """parses bam to which IDENT_TAG was added

    - FIXME replace pysam dependency with https://github.com/mdshw5/simplesam
    - FIXME test with pre 0.8 versions
    """

    query_to_hit = dict()
    samfh = pysam.Samfile(bam)
    for r in samfh:
        pwid = round(r.get_tag(ident_tag), 1)
        ref = samfh.getrname(r.tid)
        assert not query_to_hit.has_key(r.query_name)
        query_to_hit[r.query_name] = BestHit(seqid=ref, pwid=pwid, tax=None)
    samfh.close()
    return query_to_hit


def main():
    """The main function
    """

    parser = argparse.ArgumentParser(
        description=os.path.basename(sys.argv[0]) + ": " + __doc__)
    parser.add_argument('-q', "--query-file", required=True,
                        help="Fasta file containing the EMIRGE reconstruct sequences used to match greengenes")
    parser.add_argument('-i', "--hit-file", required=True,
                        help="Greengenes hits produced with BLAST (best hit only as csv) or BAM with {} tag".format(IDENT_TAG))
    parser.add_argument('-t', "--green-tax", required=True,
                        help="Greengenes taxonomy file")
    parser.add_argument('-o', "--out-table", required=True,
                        help='Output table')
    parser.add_argument('-y', "--hit-type", choices=['csv', 'bam'],
                        help='Format of hits file (automatically guessed if not provided)')

    parser.add_argument('--verbose', action="store_true",
                        help='Enable verbose output')
    parser.add_argument('--debug', action="store_true",
                        help='Enable debugging output')
    args = parser.parse_args()


    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)


    for f in [args.query_file, args.green_tax, args.hit_file]:
        if not os.path.exists(f):
            LOG.fatal("Non-existant file {}".format(f))
            sys.exit(1)

    for f in [args.out_table]:
        if os.path.exists(f):
            LOG.fatal("Cowardly refusing to overwrite {}".format(f))
            sys.exit(1)


    LOG.info("Load taxonomy from  {}".format(args.green_tax))
    gg_tax = parse_greengenes_taxonomy(args.green_tax)
    if args.hit_type == 'blast':
        p = parse_best_blast_hit
    elif args.hit_type is None and args.hit_file.endswith(".csv"):
        p = parse_best_blast_hit
    elif args.hit_type == 'bam':
        p = parse_bam
    elif args.hit_type is None and args.hit_file.endswith(".bam"):
        p = parse_bam
    else:
        raise ValueError(args.hit_type, args.hit_file)
    query_to_hit = p(args.hit_file)



    LOG.info("Assigning taxonomy to hits in {}".format(args.query_file))
    with open(args.query_file) as fh, open(args.out_table, 'w') as fhout:
        for (s_id, s_seq) in fasta_iter(fh):
            # emirge id's look as follows:
            # >47|AJ704791.1.1593 Prior=0.045024 Length=744 NormPrior=0.045018
            s_id_split = s_id.split()
            abundance = float(s_id_split[-1].replace("NormPrior=", ""))
            s_id = s_id_split[0]
            besthit = query_to_hit.get(s_id, None)
            # since blast was run against gg otus we can read from gg tax
            if besthit is not None:
                tax = gg_tax.get(besthit.seqid, None)
                besthit = besthit._replace(tax=tax)
                fhout.write("{}\t{}\t{}\t{}\t{}\n".format(
                    s_id, abundance, besthit.seqid, besthit.pwid,
                    '\t'.join(["{}:{}".format(k, v) for (k, v) in besthit.tax.items()])))
            else:
                fhout.write("{}\tNA\n".format(s_id))


if __name__ == '__main__':
    main()
