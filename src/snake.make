import shutil
from snakemake.utils import report

# simulate a bash login shell, see https://bitbucket.org/johanneskoester/snakemake/wiki/FAQ
shell.executable("/bin/bash")
# "unofficial bash strict mode" http://www.redsymbol.net/articles/unofficial-bash-strict-mode/

shell.prefix("source ~/.bashrc; set -euo pipefail;")


localrules: emirge_rename, emirge_trim_primer, classify, rarefaction, convert_tables, final


# emirge needs bowtie and usearch
assert shutil.which("usearch")# 6.0.307 ok
assert shutil.which("bowtie")# 0.12.8 ok


if config['USE_AMPLICON']:
   EMIRGE = config['EMIRGE_BASEDIR'] + '/' + 'emirge_amplicon.py'
else:
   EMIRGE = config['EMIRGE_BASEDIR'] + '/' + 'emirge.py'
EMIRGE_RENAME = config['EMIRGE_BASEDIR'] + '/' + 'emirge_rename_fasta.py'


# must be first rule

rule final:
  input:   'convert_tables.succeeded', 'results/report.html', 'results/rarefaction.pdf'
  message: 'This is the end. My only friend, the end'
  output: 'COMPLETE'
  shell:  'touch {output}'


rule trim_fastq:
  input:  fq1='{sample}R1.fastq.gz', fq2='{sample}R2.fastq.gz'
  output: fq1=temp('{sample}R1.trimmed.fastq.gz'), fq2=temp('{sample}R2.trimmed.fastq.gz')
  params: min3qual='3', minlen='60'
  shell:  '{config[FAMAS]} -i {input.fq1} -j {input.fq2} -o {output.fq1} -p {output.fq2} -q {params.min3qual} -l {params.minlen}'


rule concat_fastq:
  input:  fqs1=expand('{sample}R1.trimmed.fastq.gz', sample=config['SAMPLES']),
          fqs2=expand('{sample}R2.trimmed.fastq.gz', sample=config['SAMPLES'])
  output: fq1=temp('R1.trimmed.fastq.gz'), fq2=temp('R2.trimmed.fastq.gz')
  shell:  'zcat {input.fqs1} | gzip > {output.fq1} && zcat {input.fqs2} | gzip > {output.fq2}'
  #output: fq1=temp('R1.trimmed.fastq.gz'), fq2=temp('R2.trimmed.fastq.gz')
  #shell:  'zcat {input.fqs1} | gzip > {output.fq1} && zcat {input.fqs2} | gzip > {output.fq2}'


rule pre_filter:
  input:    fq1=rules.concat_fastq.output.fq1, fq2=rules.concat_fastq.output.fq2
  output:   fq1=temp('R1.trimmed.flt.fastq.gz'), fq2=temp('R2.trimmed.flt.fastq'), ratios='results/ratios.txt'
  params:   ssu_fa=config['SSU_FA'], spikein_name=config['SPIKEIN-NAME']
  threads:  8
  shell:    '{config[BWA]} mem -t {threads} -k 30 -M {params.ssu_fa} {input.fq1} {input.fq2} | {config[PREFILTER]} -i - -1 {output.fq1} -2 {output.fq2} -s {params.spikein_name} > {output.ratios}'


rule emirge:
  input:     fq1=rules.pre_filter.output.fq1, fq2=rules.pre_filter.output.fq2
  threads:   8
  output:    'emirge.succeeded'
  benchmark: "emirge-benchmark.json"
  params :    max_read_len=config['MAX_READ_LEN'], ins_size=config['INS_SIZE'], ins_stdev=config['INS_STDEV'], ssu_fa=config['SSU_FA'], ssu_db=config['SSU_DB'],
              iter_arg="-n 1" if config['DEBUG'] else "-n 10"
  #params:    max_read_len=config['MAX_READ_LEN'], ssu_fa=config['SSU_FA']
  # existing emirge directory can't be reused so delete if existing
  shell:     'test -d emirge && rm -rf emirge; {EMIRGE} -l {params.max_read_len} -i {params.ins_size} -s {params.ins_stdev} --phred33 {params.iter_arg} -a {threads} emirge -f {params.ssu_fa} -b {params.ssu_db} -1 {input.fq1} -2 {input.fq2} && touch {output}'


rule emirge_rename:
  input:  rules.emirge.output
  output: 'emirge_out.fa'
  run:
    import glob
    import subprocess
    import fnmatch
    import os

    last_iter = sorted(glob.glob(os.path.join("emirge", "iter.*")))[-1]
    cmd = [EMIRGE_RENAME, last_iter]
    print("Using {}".format(last_iter))
    with open(output[0], "w") as fh:
        subprocess.call(cmd, stdout=fh)
    # clean BAMs from emirge. FIXME make sep target?
    if not config['DEBUG']:
        for root, dirnames, filenames in os.walk('emirge'):
            for filename in fnmatch.filter(filenames, '*.bam'):
                os.unlink(os.path.join(root, filename))


rule emirge_trim_primer:
  input:  rules.emirge_rename.output
  output: 'emirge_outprimer_trimmed.fa'
  params: minlen='200', maxlen='1400'
  shell:  '{config[PRIMER_TRIMMER]} -f -i {input} -o {output} --minlen {params.minlen} --maxlen {params.maxlen}'


rule emirge_vs_gg:
  input:
        fa=rules.emirge_trim_primer.output,
	ref=config['GG_REF']
  #output: 'greengenes-hits-blast.csv'
  output:  'greengenes-hits-graphmap.bam'
  threads: 8
  shell:  '{config[GRAPHMAP]} -t {threads} -r {input.ref} -d {input.fa} | {config[IDENT_TO_BAM]}  - {output} {input.ref}'


rule classify:
  input: query=rules.emirge_trim_primer.output, hits=rules.emirge_vs_gg.output, gg_tax=config['GG_TAX']
  output: 'results/raw-table.txt'
  message: 'classifying hits'
  shell: '{config[CLASSIFY_HITS]} -f -q {input.query} -i {input.hits} -t {input.gg_tax} -o {output}'


rule rarefaction:
  input: 'results/raw-table.txt'
  output: 'results/rarefaction.pdf'
  message: 'Plotting rarefaction curve'
  shell: '{config[PLOT_RAREFACTION]} -o {output} -t {input} -f'

rule convert_tables:
  input:  rules.classify.output
  output: 'convert_tables.succeeded'
  message: 'converting raw table'
  shell:  '{config[CONVERT_TABLE]} {input} results/ && touch {output}'


rule report:
  input:  ratios=rules.pre_filter.output.ratios,
          rawtable=rules.classify.output[0],
  output: html="results/report.html"
  run:    report("""
          ===================
          16S Pipeline Report
          ===================

          Ratios of unspecific products, spike-in
          ({config[SPIKEIN-NAME]}) and 16S are listed in the
          ratios_ file.

          All raw results can be found in file rawtable_

          Abundances are scaled to 100% and are only based on
          valid 16S sequences. Abundances for each taxonomic rank are listed in the files called
	      rank-table.csv and (quick and dirty piecharts) in rank-piechart.pdf
          """, output.html, metadata="Andreas WILM", **input)
