# FIXME
# - produce final table
# - ensure usearch and bowtie are in PATH as well
# - all vars before first rule should go into config https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-configuration
# - add more benchmark rules: https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-benchmark-rules
# - add cluster config rules per target e.g. params: runtime="4h" and use in cluster arg
# qsub="qsub -pe OpenMP 1 -l mem_free=1G -l h_rt=48:00:00  -j y -V -b y -cwd"
# $qsub -N snakemake -o snakemake.qsub.log  'snakemake  --forceall -j 8 -c "qsub -pe OpenMP {threads} -l mem_free=12G -l h_rt=06:00:00  -j y -V -b y -cwd" -s snake.make --configfile conf.json'


# simulate a bash login shell, see https://bitbucket.org/johanneskoester/snakemake/wiki/FAQ
shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; set -o pipefail; ")

if config['USE_AMPLICON']:
   EMIRGE = config['EMIRGE_BASEDIR'] + '/' + 'emirge_amplicon.py'
else:
   EMIRGE = config['EMIRGE_BASEDIR'] + '/' + 'emirge.py'
EMIRGE_RENAME = config['EMIRGE_BASEDIR'] + '/' + 'emirge_rename_fasta.py'


# must be first rule

rule final:
  #input:  'blast-greengenes.csv'
  input:   'raw-table.txt'
  message: 'This is the end. My only friend, the end'


rule filter_fastq:
  input:  fq1='{sample}R1.fastq.gz', fq2='{sample}R2.fastq.gz'
  output: fq1=temp('{sample}R1.flt.fastq.gz'), fq2=temp('{sample}R2.flt.fastq.gz')
  params: min3qual='3', minlen='60'
  shell:  '{config[FAMAS]} -i {input.fq1} -j {input.fq2} -o {output.fq1} -p {output.fq2} -q {params.min3qual} -l {params.minlen}'


rule concat_fastq:
  input:  fqs1=expand('{sample}R1.flt.fastq.gz', sample=config['SAMPLES']),
          fqs2=expand('{sample}R2.flt.fastq.gz', sample=config['SAMPLES'])
  output: fq1=temp('R1.flt.fastq.gz'), fq2=temp('R2.flt.fastq')
  shell:  'zcat {input.fqs1} | gzip > {output.fq1} && zcat {input.fqs2} > {output.fq2}'


rule emirge:
  input:     fq1=rules.concat_fastq.output.fq1, fq2=rules.concat_fastq.output.fq2
  # FIXME make max_num_threads
  threads:   8
  output:    'emirge_succeeded'
  message:   'WARN: one iteration only!'
  benchmark: "emirge-benchmark.json"
  params: max_read_len=config['MAX_READ_LEN'], ins_size=config['INS_SIZE'], ins_stdev=config['INS_STDEV'], ssu_fa=config['SSU_FA'], ssu_db=config['SSU_DB']
  # existing emirge directory can't be reused so delete if existing
  shell:     'test -d emirge && rm -rf emirge; {EMIRGE} -l {params.max_read_len} -i {params.ins_size} -s {params.ins_stdev} --phred33 -n 1 -a {threads} emirge -f {params.ssu_fa} -b {params.ssu_db} -1 {input.fq1} -2 {input.fq2} && touch {output}'


rule emirge_rename:
  input: rules.emirge.output
  output: 'emirge_out.fa'
  run:
    import glob
    import subprocess
    last_iter = sorted(glob.glob(os.path.join("emirge", "iter.*")))[-1]
    cmd = [EMIRGE_RENAME, last_iter]
    with open(output[0], "w") as fh:
        subprocess.call(cmd, stdout=fh)


rule emirge_trim_primer:
  input: rules.emirge_rename.output
  output: 'emirge_outprimer_trimmed.fa'
  params: minlen='200', maxlen='1400'
  shell: '{config[PRIMER_TRIMMER]} -i {input} -o {output} --minlen {params.minlen} --maxlen {params.maxlen}'


rule emirge_vs_gg:
  input:  rules.emirge_trim_primer.output, ref=config['GG_REF']
  #output: 'greengenes-hits-blast.csv'
  output: 'greengenes-hits-graphmap.bam'
  # FIXME make max_num_threads
  threads: 8
  shell:  '{config[GRAPHMAP]} -x illumina -t {threads} -r {input.ref} -d {input} | {config[IDENT_TO_BAM]}  - {output} {input.ref}'


rule classify:
  input: query=rules.emirge_trim_primer.output, hits=rules.emirge_vs_gg.output, gg_tax=config['GG_TAX']
  output: 'raw-table.txt'
  shell: '{config[CLASSIFY_HITS]} -q {input.query} -i {input.hits} -t {input.gg_tax} -o {output}'
