#!/mnt/software/unstowable/anaconda/envs/pysam-0.8.2/bin/python
#!/usr/bin/env python
"""Adds pairwise identity tag (matches/readlength) to BAM file

TODO:
- add tests

- use argparse
- allow user defined tag as well as overwriting of tag
"""


import sys
import os



import pysam
PYSAM_VERSION = [int(x) for x in pysam.__version__.split(".")]
assert PYSAM_VERSION >= [0, 8, 0]


IDENT_TAG = 'Xi'


def calc_ident(r, ref):
    """Computes pairwise identity of read r
    """

    if r.is_unmapped:
        return -1.0
    op_counts = {'indel': 0, 'match': 0, 'mismatch': 0}
    for (qpos, rpos) in r.aligned_pairs:
        if qpos is None or rpos is None:
            op = 'indel'
        else:
            assert rpos < len(ref)
            if ref[rpos].upper() == r.seq[qpos].upper():
                op = 'match'
            else:
                op = 'mismatch'
        op_counts[op] += + 1
    ident = op_counts['match']/float(sum(op_counts.values()))
    return ident


def main(sam_fh_in, sam_fh_out, fasta_fh):
    """FIXME
    """
    last_tid = -1
    for r in sam_fh_in:
        if r.tid != last_tid:
            ref = fasta_fh.fetch(sam_fh_in.getrname(r.tid))
            last_tid = r.tid
        ident = calc_ident(r, ref)

        # NOTE pysam's set_tag inferred value_type 'd' which is
        # undefined according to
        # https://samtools.github.io/hts-specs/SAMv1.pdf
        r.set_tag(IDENT_TAG, round(ident*100, 1), value_type='f', replace=REPLACE_TAG)
        sam_fh_out.write(r)


if __name__ == "__main__":
    REPLACE_TAG = False

    try:
        sam_in, sam_out, ref_fa = sys.argv[1:]
    except IndexError:
        sys.stderr.write("FATAL: Need input and output BAM (stdout supported)  as well as the (indexed) reference as (only) arguments\n")
        sys.exit(1)
    assert not os.path.exists(sam_out)

    fasta_fh = pysam.Fastafile(ref_fa)
    sam_fh_in = pysam.Samfile(sam_in)# mode automatically inferred
    out_mode = 'w'
    if sam_out.endswith(".bam"):
        out_mode += "b"
    sam_fh_out = pysam.Samfile(sam_out, out_mode, template=sam_fh_in)
    main(sam_fh_in, sam_fh_out, fasta_fh)
