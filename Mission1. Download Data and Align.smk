import os

SAMPLES = """clip-35L33G
polya-siLin28a
polya-siLuc
polya-untreated
rpf-siLin28a
rpf-siLuc""".split()

# Genome Index
GENOME_INDEX = "reference-data/genome/mouse/mm39/genome-index/"
STAR_CONTAM = "~/genome_index/star_mouse_contam"
STAR_MM39_TX = os.path.join(GENOME_INDEX, "star-mm39-gencode-transcript")
# STAR_MM39_54 = os.path.join(GENOME_INDEX, "star-mm39-gencode-54nt")
# STAR_MM39_78 = os.path.join(GENOME_INDEX, "star-mm39-gencode-78nt")
# STAR_MM39_36 = os.path.join(GENOME_INDEX, "star-mm39-gencode-36nt")
"""
STAR_GENOME = {
    'clip-35L33G': STAR_MM39_78,
    'polya-untreated': STAR_MM39_54,
    'rpf-siLuc': STAR_MM39_36,
    'rpf-siLin28a': STAR_MM39_36,
    'polya-siLuc': STAR_MM39_54,
    'polya-siLin28a': STAR_MM39_54
}
"""
# Genome References
# MM39_GTF = "reference-data/genome/mouse/mm39/gencode/gencode.vM28.annotation.gtf"
MM39_TX_FASTA = "reference-data/genome/mouse/mm39/gencode/gencode.vM28.transcripts.fa"
# MM39_PRI_FASTA = "reference-data/genome/mouse/mm39/gencode/GRCm39.primary_assembly.genome.fa"

# Tools
# TOOL_PATH="~/tools"

# Cutadapt
modban_illumina = 'CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTG'
illumina_SRA15 = 'ATCTCGTATGCCGTCTTCTGCTTG'
illumina_TruSeq = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC'

ADAPTER_SEQ = {
    'clip-35L33G': modban_illumina,
    'polya-untreated': modban_illumina,
    'rpf-siLuc': illumina_SRA15,
    'rpf-siLin28a': illumina_SRA15,
    'polya-siLuc': illumina_SRA15,
    'polya-siLin28a': illumina_SRA15
}

QUALITY=25
MINLEN=20
TIMES=2

TARGETS=[
#    expand("filtered/{sample}.filtered.fq.gz", sample=SAMPLES),
#    expand("alignments/{sample}.contam.bam", sample=SAMPLES),
    expand("alignments/{sample}.bam", sample=SAMPLES)
#    expand(f"{PROJECT_PATH}" + "/alignments/{sample}.cram.crai", sample=SAMPLES)
#    expand(f"{PROJECT_PATH}" + "/annotations/{sample}.bedintersect.gz", sample=SAMPLES)
    ]

localrules: all

rule all:
    input: TARGETS

rule adapter_trimming:
# And discarding reads with row quality end
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "filtered/{sample}.filtered.fq.gz"
    threads: 5
    run:
        adapter_seq = ADAPTER_SEQ[wildcards.sample]
        shell("""
            cutadapt \
                -a {adapter_seq} -n {TIMES} \
                -q {QUALITY},{QUALITY} \
                --minimum-length {MINLEN} \
                -o {output} \
                -j {threads} \
                {input}
              """)

rule map_to_contam:
# Illumina PRIMERS & ADAPTERS, rRNA, etc.
    input:
        "filtered/{sample}.filtered.fq.gz"
    output:
        "alignments/{sample}.contam.bam"
    threads: 10
    run:
        temp_dir = str(output).replace('.contam.bam', '.contam.temp_dir')
        if wildcards.sample.startswith('clip'):
            ratio = 0.1
        elif wildcards.sample.startswith('rpf'):
            ratio = 0.05
        elif wildcards.sample.startswith('polya'):
            ratio = 0.05
        else:
            raise NameError(wildcards.sample)
        shell("""
            STAR --readFilesIn {input} \
              --genomeDir {STAR_CONTAM} \
              --runThreadN {threads} \
              --readFilesCommand zcat \
              --outStd SAM \
              --outSAMunmapped Within \
              --outFilterMismatchNoverLmax {ratio} \
              --outTmpDir {temp_dir} \
              | samtools sort -n -@ 6 - -o {output}
              """)

rule collect_unmapped_to_contam:
    input:
        'alignments/{sample}.contam.bam'
    output:
        'alignments/{sample}.unmapped.fq'
    threads: 2
    shell: "samtools view -h -f 4 {input} | bedtools bamtofastq -i - -fq {output}"

rule map_unmapped_to_genome:
    input: 'alignments/{sample}.unmapped.fq'
    output: 'alignments/{sample}.bam'
    threads: 10
    run:
        temp_dir = str(output).replace('.bam', '.temp_dir')
        indx = STAR_MM39_TX
        if wildcards.sample.startswith('clip'):
            ratio = 0.1
        elif wildcards.sample.startswith('rpf'):
            ratio = 0.05
        elif wildcards.sample.startswith('polya'):
            ratio = 0.05
        else:
            raise NameError(wildcards.sample)
        shell("""
            STAR --readFilesIn {input} \
              --genomeDir {indx} \
              --runThreadN {threads} \
              --outStd SAM \
              --outSAMunmapped Within \
              --outFilterMismatchNoverLmax {ratio} \
              --outTmpDir {temp_dir} \
              | samtools sort -n -@ 6 - -o {output}
              """)
