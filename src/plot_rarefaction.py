#!/mnt/software/unstowable/anaconda/bin/python
"""Plots a pseudo-rarefaction curve which is based on given abundances and detection limit

Works for Python2 and 3
"""

import os
import sys
import argparse

# non interactive
import matplotlib as mpl; mpl.use('Agg')

from matplotlib import pyplot as plt


DEFAULT_DETECT_LIM = 10**-4


def main(raw_table, detect_lim, outplot, sample_name=None):
    """main function"""
    abundances = []
    with open(raw_table) as fh:
        for line in fh:
            abd = float(line.split()[1])
            abundances.append(abd)
    #print("DEBUG min abundance={}".format(min(abundances)))

    plotx = []
    ploty = []

    step = 2
    for rate in range(step, 100, step):
        #for x in abundances:
        #    print("x={} x*rate={} downsampled={} >{}={}".format(x, x*rate, x*rate/100.0, detect_lim, x*rate/100.0 >= detect_lim))
        num_otus_detected = len([x for x in abundances if x*rate/100.0 >= detect_lim])
        #print("DEBUG num_otus_detected={} at rate {}".format(num_otus_detected, rate))
        plotx.append(rate)
        ploty.append(num_otus_detected)

    #print("DEBUG detect_lim={}".format(detect_lim))
    #print("DEBUG abundances={}".format(abundances))

    plt.plot(plotx, ploty, '-', color='red', linewidth=3.0)
    plt.axhline(len(abundances), color='blue', ls="--", linewidth=3.0)
    plt.xlabel("Downsampling Rate [%]")
    plt.ylabel("# Detected OTUs")
    title = "Rarefaction Curve"
    if sample_name:
        title += " for {}".format(sample_name)
    plt.title(title)
    plt.ylim(0, 1.05*max(ploty))
    #plt.xscale('log', base=10)
    plt.savefig(outplot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=os.path.basename(sys.argv[0]) + ": " + __doc__)

    parser.add_argument("-o", "--outplot", required=True,
                        help="Output plot")
    parser.add_argument("-t", "--table", required=True,
                        help="Table with abundances")
    default = DEFAULT_DETECT_LIM
    parser.add_argument("-l", "--det-lim", default=default,
                        help="Detection limit (default {})".format(default))
    parser.add_argument("-n", "--sample-name",
                        help="Sample name (added to plot)")
    parser.add_argument("-f", "--overwrite",
                        action="store_true",
                        help="Force overwrite")
    args = parser.parse_args()
    if os.path.exists(args.outplot) and not args.overwrite:
        sys.stderr.write("FATAL: Cowardly refusing to overwrite {}\n".format(args.outplot))
        sys.exit(1)
    if not os.path.exists(args.table):
        sys.stderr.write("FATAL: Missing input file {}\n".format(args.table))
        sys.exit(1)

    main(args.table, args.det_lim, args.outplot, sample_name=args.sample_name)
