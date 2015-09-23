This is a remake of the original 16S pipeline, designed for Illumina
shotgun sequencing of 16S rRNA amplicon sequences
([Ong et al., 2013, PMID 23579286](http://www.ncbi.nlm.nih.gov/pubmed/23579286)).

## Running the pipeline

The pipeline is based on
[snakemake](http://www.ncbi.nlm.nih.gov/pubmed/22908215). The main
program (`S16_V2.py`) will write a config file (`conf.json`) and
snakemake file (`snake.make`) in the given output directory. These are
then used to call snakemake via `qsub` using the also created wrapper
script `snake.sh`. The system will send an email to you upon
completion (be it successful or not).

For help see `S16_V2.py --help`.

Only upon successful completion the output directory will contain an
empty file called `COMPLETE`.

Results (abundance tables and piecharts) can then be found in
`results` subdirectory (see `report.html` there).


## Steps involved

- Reads are first preprocessed (merging of multiple files, trimming
etc.) with [famas](https://github.com/andreas-wilm/famas).
- Preprocessed reads are classified into 'unspecific', 'spike-in' and '16S', based on a
[BWA-MEM](http://arxiv.org/abs/1303.3997) mapping against the original
EMIRGE/ARB SSU database containing the spike-in (see `ratios.txt` for
corresponding ratios).
- Only 16S sequences are used for reconstructing the full-length
sequences with EMIRGE
([Miller et al., 2011, PMID 21595876](http://www.ncbi.nlm.nih.gov/pubmed/21595876))
or EMIRGE amplicon
([Miller et al., 2013; PMID 23405248](http://www.ncbi.nlm.nih.gov/pubmed/23405248))
- Primers are clipped from the reconstructed sequences
- Clipped sequences are  mapped with
[Graphmap](http://biorxiv.org/content/early/2015/06/10/020719) against
a preclustered version (99% OTU) of the Greengenes database for
classification.
- Abundance-tables and piecharts are  created in the `results`
subdirectory of the output directory. Pairwise identity thresholds
for the different taxonomi ranks are implemented as determined by
[Yarza et al. (2014; PMID 25118885)](http://www.ncbi.nlm.nih.gov/pubmed/25118885).


## Tip

If you want to first see what the pipeline would do, call `S16_V2.py` with `--no-run`.
Then check the created files (see above). To get a graphical representation of the workflow, run (from the output directory):

    snakemake -s snake.make --configfile conf.json --dag --forceall | dot -Tpdf > dag.pdf

and have a look at `dag.pdf`. Once you're satisfied just run `bash snake.sh`.

