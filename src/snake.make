# FIXME
# - add benchmark rules: https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-benchmark-rules
# - add cluster config rules per target e.g. params: runtime="4h" and use in cluster arg
# - make sure usearch and bowtie are in PATH as well
# - all vars before first rule should go into config https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-configuration

# simulate a bash login shell, see https://bitbucket.org/johanneskoester/snakemake/wiki/FAQ
shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; set -o pipefail; ")

USE_AMPLICON = True

SAMPLES = ['1_' ,'2_' ,'3_' ,'4_' ,'5_' ,'6_']

GG_REF = '/mnt/projects/wilma/16s/S16_V2/greengenes/99_otus.fasta'
GG_TAX = '/mnt/genomeDB/misc/greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus/taxonomy/99_otu_taxonomy.txt'

SSU_FA = '/mnt/genomeDB/misc/softwareDB/emirge/SSU_candidate_db.fasta'
SSU_DB = '/mnt/genomeDB/misc/softwareDB/emirge/SSU_candidate_db'

FAMAS = '/mnt/software/stow/famas-0.0.7/bin/famas'
EMIRGE_BASEDIR = '/mnt/software/stow/emirge-v0.60-15-g0ddae1c-wilma/bin/'
EMIRGE_RENAME = EMIRGE_BASEDIR + '/' + 'emirge_rename_fasta.py'
if USE_AMPLICON:
   EMIRGE = EMIRGE_BASEDIR + '/' + 'emirge_amplicon.py'
else:
   EMIRGE = EMIRGE_BASEDIR + '/' + 'emirge.py'

MAX_READ_LEN = 100
INS_SIZE = 200
INS_STDEV = 40

PRIMER_TRIMMER = '/mnt/projects/wilma/16s/S16_V2/src/primer_trimmer.py'
GRAPHMAP = '/mnt/software/stow/graphmap-0.2.2-dev-604a386/bin/graphmap'
IDENT_TO_BAM = '/mnt/projects/wilma/16s/S16_V2/src/ident_to_bam.py'
CLASSIFY_HITS = '/mnt/projects/wilma/16s/S16_V2/src/classify_hits.py'


# must be first rule

rule final:
  #input: 'blast-greengenes.csv'
  input: 'raw-table.txt'
  message: 'This is the end. My only friend, the end'


rule filter_fastq:
  input: fq1='{sample}R1.fastq.gz', fq2='{sample}R2.fastq.gz'
  output: fq1=temp('{sample}R1.flt.fastq.gz'), fq2=temp('{sample}R2.flt.fastq.gz')
  shell: '{FAMAS} -i {input.fq1} -j {input.fq2} -o {output.fq1} -p {output.fq2} -q 3 -l 60  --append'


rule concat_fastq:
  input: fqs1=expand('{sample}R1.flt.fastq.gz', sample=SAMPLES), fqs2=expand('{sample}R2.flt.fastq.gz', sample=SAMPLES)
  output: fq1=temp('R1.flt.fastq.gz'), fq2=temp('R2.flt.fastq')
  shell: 'zcat {input.fqs1} | gzip > {output.fq1} && zcat {input.fqs2} > {output.fq2}'


rule emirge:
  input: fq1=rules.concat_fastq.output.fq1, fq2=rules.concat_fastq.output.fq2
  threads: 8
  output: 'emirge_succeeded'
  message: 'WARN: one iteration only!'
  # existing emirge directory can't be reused so delete if existing
  shell:
    'test -d emirge && rm -rf emirge; {EMIRGE} -l {MAX_READ_LEN} -i {INS_SIZE} -s {INS_STDEV} --phred33 -n 1 -a {threads} emirge -f {SSU_FA} -b {SSU_DB} -1 {input.fq1} -2 {input.fq2} && touch {output}'


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
  shell: '{PRIMER_TRIMMER} -i {input} -o {output} --minlen 200 --maxlen 1400'


rule emirge_vs_gg:
  input: rules.emirge_trim_primer.output, ref=GG_REF
  #output: 'greengenes-hits-blast.csv'
  output: 'greengenes-hits-graphmap.bam'
  threads: 8
  shell: '{GRAPHMAP} -x illumina -t {threads} -r {input.ref} -d {input} | {IDENT_TO_BAM}  - {output} {input.ref}'


rule classify:
  input: query=rules.emirge_trim_primer.output, hits=rules.emirge_vs_gg.output, gg_tax=GG_TAX
  output: 'raw-table.txt'
  shell: '{CLASSIFY_HITS} -q {input.query} -i {input.hits} -t {input.gg_tax} -o {output}'
