#!/usr/bin/env python
"""Primer chopping for the 16S project.

Uses difflib for locating the subsequence best (min 0.6 identity; see
http://docs.python.org/library/difflib.html) matching the predefined
primers. reported if within a certain sequence range and obviously if
the primer orientation is correct. Ambiguity characters are treated as
mismatches.
"""

#--- standard library imports
#
import logging
import difflib
import sys
import os
from optparse import OptionParser
import gzip

#--- third-party imports
#
from Bio import SeqIO

#--- project specific imports
#
# /

__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

PRIMER = dict()
PRIMER['fw'] = 'ACTCCTACGGGAGGCAGC'
PRIMER['rv'] = "GTCGTCAGCTCGTGYYG"
#PRIMER['rv'] = 'GTCGTCAGCTCGTG'
DEFAULT_MIN_LEN = 500
DEFAULT_MAX_LEN = 1000


def guess_seqformat(fseq):
    """
    FIXME
    """
    default = 'fasta'

    # Table for guessing the alignment format from the file extension. 
    # See http://www.biopython.org/wiki/SeqIO
    #
    # Only define the ones I actually came accors here:
    ext_to_fmt_table = dict(
        aln = 'clustal',
        embl = 'embl',
        fasta = 'fasta',
        fa = 'fasta',
        genbank = 'genbank',
        gb = 'genbank',
        phylip = 'phylip',
        phy = 'phylip',
        ph = 'phylip',
        pir = 'pir',
        stockholm = 'stockholm',
        st = 'stockholm',
        stk = 'stockholm')

    try:
        fext = os.path.splitext(fseq)[1]
        fext = fext[1:].lower()
        fmt =  ext_to_fmt_table[fext]
    except KeyError:
        return default
    return fmt


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: " + __doc__ + "\n" \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="debugging")
    parser.add_option("-i", "--input",
                      dest="fseq_in", # type="string|int|float"
                      help="Sequence input file (gzip supported)")
    parser.add_option("-o", "--output",
                      dest="fseq_out", # type="string|int|float"
                      help="Sequence output file")
    parser.add_option("", "--minlen",
                      dest="minlen",
                      type="int",
                      default=DEFAULT_MIN_LEN,
                      help="Minimum allowed sequence length after chopping (default: %d)" % DEFAULT_MIN_LEN)
    parser.add_option("", "--maxlen",
                      dest="maxlen",
                      type="int",
                      default=DEFAULT_MAX_LEN,
                      help="Maximum allowed sequence length after chopping (default: %d)" % DEFAULT_MAX_LEN)
    return parser




def main():
    """
    The main function
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)
        
    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        
    if not opts.fseq_in:
        parser.error("Input file argument missing.")
        sys.exit(1)
    if not os.path.exists(opts.fseq_in):
        sys.stderr.write(
            "Input file '%s' does not exist.\n" % opts.fseq_in)
        sys.exit(1)

    if not opts.fseq_out:
        parser.error("Output file argument missing.")
        sys.exit(1)
    if os.path.exists(opts.fseq_out):
        sys.stderr.write(
            "Cowardly refusing to overwrite existing file '%s'.\n" % opts.fseq_out)
        sys.exit(1)

    if opts.maxlen < opts.minlen:
        sys.stderr.write(
            "Maximum length %d is smaller than minimum length %d.\n" % (opts.minlen, opts.maxlen))
        sys.exit(1)
        

    LOG.info("Trimming sequences in '%s' and writing to '%s'. Allowed length range = %d-%d." % (
            opts.fseq_in, opts.fseq_out, opts.minlen, opts.maxlen))
    if opts.fseq_in[-3:] == ".gz":
        fh_in = gzip.open(opts.fseq_in, 'r')
    else:
        fh_in = open(opts.fseq_in, 'r')

    fh_out = open(opts.fseq_out, 'w')
    seq_fmt = guess_seqformat(opts.fseq_in)
    for seqrec in SeqIO.parse(fh_in, seq_fmt):
        pos = dict()
        for ori in ['fw', 'rv']:
            primer = PRIMER[ori]

            # creates words for all possible windows. huge speed
            # impact. needs optimisation
            all_seq_words = [str(seqrec.seq[i:i+len(primer)]) for i in range(len(seqrec.seq))]
            pos[ori] = seqrec.seq.find(primer)
            if pos[ori] == -1:
                closest_matches = difflib.get_close_matches(primer, all_seq_words)
                if len(closest_matches) == 0:
                    LOG.warn("Skipping %s because not even a fuzzy %s match found" % (seqrec.id, ori))
                    continue            
                pos[ori] = seqrec.seq.find(closest_matches[0])
    
        if pos['fw'] > pos['rv']:
            LOG.info("Skipping %s because of primer matches with wrong orientation" % (seqrec.id))
            continue
    
        outseq = str(seqrec.seq[pos['fw']:pos['rv']+len(PRIMER['rv'])])
        if len(outseq) < opts.minlen or len(outseq) > opts.maxlen:
            LOG.info("Skipping seq '%s' because trimmed sequence length (%d) outside of allowed range (min %d; max %d)" % (
                    seqrec.id, len(outseq), opts.minlen, opts.maxlen))
            continue        

        # lazy! fasta only. could use biopython instead
        # print description which keeps emirges prior==abundance
        fh_out.write(">%s\n%s\n" % (seqrec.description, outseq))
    fh_out.close()
    fh_in.close()

if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
