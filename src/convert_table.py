#!/mnt/software/unstowable/anaconda/bin/python
"""FIXME
"""

# FIXME hardcoded interpreter path
# can't find matplotlib otherwise?!

from collections import OrderedDict
import logging
import sys

# make sure we don't need X for plotting
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt


__author__ = "Andreas Wilm"
__version__ = "2.0.0a"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "The MIT License (MIT)"


# global logger
LOG = logging.getLogger("")
logging.basicConfig(
    level=logging.WARN,
    format='%(levelname)s [%(asctime)s]: %(message)s')


# Yarza et al. (2014): http://www.ncbi.nlm.nih.gov/pubmed/25118885
PWID_THRESHOLD = OrderedDict([('p', 75), ('c', 78.5), ('o', 82),
                              ('f', 86.5), ('g', 94.5), ('s', 97)])

TAXA_ABBREV_MAP = {'k':'Kingdom', 'p':'Phylum', 'c':'Class',
                   'o':'Order', 'f':'Family', 'g':'Genus', 's':'Species'}


def main(raw_table, outprefix):
    """main function
    """

    piechart = OrderedDict()
    for k, v in PWID_THRESHOLD.items():
        piechart[k] = dict()


    with open(raw_table) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            ls = line.rstrip().split("\t")
            (_emirge_name, abundance, _gg_name, gg_pwid), tax = ls[:4], ls[4:]
            abundance, gg_pwid = float(abundance), float(gg_pwid)
            tax = OrderedDict([x.split(":") for x in tax])

            # clean any entries below given threshold, i.e. if we have a hit at 90%
            # remove genus and species
            for k, v in PWID_THRESHOLD.items():
                if gg_pwid < v:
                    tax[k] = ''

            # add genus to species if both are set
            if tax['g'] != '' and tax['s'] != '':
                tax['s'] = tax['g'] + " " + tax['s']

            # add to abundance piechart
            for k, v in PWID_THRESHOLD.items():
                piechart[k][tax[k]] = piechart[k].get(tax[k], 0) + abundance


    # abundances don't add up to 1 because we removed reconstructed reads
    # whose length was outside defined range. save to assume they are
    # unknown. also rename to 'NA'
    for k, v in piechart.items():
        extra_unknown = 1.0 - sum(v.values())
        v['NA'] = v.get('', 0) + extra_unknown
        if '' in v:
            del v['']

    for k, v in piechart.items():
        assert abs(sum(v.values())-1) <= 0.00000000001

        outplot = "{}{}-piechart.pdf".format(outprefix, TAXA_ABBREV_MAP[k])
        table = "{}{}-table.csv".format(outprefix, TAXA_ABBREV_MAP[k])

        with open(table, 'w') as fh:
            sizes = []
            labels = []
            for name, abd in sorted(v.items(), key=lambda x: x[1], reverse=True):
                sizes.append(abd*100)
                labels.append(name)
                fh.write("{}\t{:f}\n".format(name, abd))

        # http://matplotlib.org/examples/pie_and_polar_charts/pie_demo_features.html
        plt.pie(sizes, labels=labels, autopct='%1.4f%%', shadow=True, startangle=90)
        # Set aspect ratio to be equal so that pie is drawn as a circle.
        plt.axis('equal')
        plt.savefig(outplot)
        plt.clf()


if __name__ == '__main__':
    assert len(sys.argv) == 3
    raw_table = sys.argv[1]
    outprefix = sys.argv[2]
    main(raw_table, outprefix)
