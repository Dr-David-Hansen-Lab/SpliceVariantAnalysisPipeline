# Snakefile for SpliceVariantAnalysisPipeline
# This pipeline performs RNA-seq isoform analysis using Salmon, HISAT2, StringTie, and R-based tools.
# Edit config/config.yaml and config/samples.csv to customize for your dataset.


import os
import pandas as pd
import yaml

# Load config early so it is available to any Python executed at parse time
configfile: "config/config.yaml"

# Directory constants from config to avoid unintended wildcards
TRIMMED_DIR = config["trimmed_dir"]
SALMON_QUANT_DIR = config["salmon_quant_dir"]
ALIGNMENTS_DIR = config["alignments_dir"]
LOGS_DIR = config["logs_dir"]
HISAT2_INDEX_DIR = config["hisat2_index"]
SALMON_INDEX_DIR = config["salmon_index"]
ANALYSIS_DIR = config.get("analysis_dir", "output/analysis")

# Load sample sheet once and reuse mappings
def _read_samples_df():
	# Robust CSV parsing: skip comments, trim spaces, normalize headers
	df = pd.read_csv(config["samples_csv"], comment="#", skipinitialspace=True)
	df.columns = [str(c).strip() for c in df.columns]
	required = ["sample", "fastq_1", "fastq_2"]
	missing = [c for c in required if c not in df.columns]
	if missing:
		raise ValueError(f"samples.csv missing required columns: {missing}")
	df["sample"] = df["sample"].astype(str).str.strip()
	df["fastq_1"] = df["fastq_1"].astype(str).str.strip()
	df["fastq_2"] = df["fastq_2"].astype(str).str.strip()
	return df

SAMPLES_DF = _read_samples_df()
SAMPLES = SAMPLES_DF["sample"].tolist()
FASTQ1 = SAMPLES_DF.set_index("sample")["fastq_1"].to_dict()
FASTQ2 = SAMPLES_DF.set_index("sample")["fastq_2"].to_dict()

rule all:
	input:
		"data/reference/genome.fa",
		"data/reference/transcripts.fa",
		"data/reference/reference.gtf",
		expand(f"{TRIMMED_DIR}/{{sample}}_R1.trimmed.fastq.gz", sample=SAMPLES),
		expand(f"{TRIMMED_DIR}/{{sample}}_R2.trimmed.fastq.gz", sample=SAMPLES),
		expand(f"{SALMON_QUANT_DIR}/{{sample}}/quant.sf", sample=SAMPLES),
		expand(f"{ALIGNMENTS_DIR}/{{sample}}.sorted.bam", sample=SAMPLES),
		f"{ANALYSIS_DIR}/merged.stringtie.gtf",
		f"{ANALYSIS_DIR}/tx2gene.tsv",
		expand(f"{ANALYSIS_DIR}/isoforms/{{sample}}_proportions.csv", sample=SAMPLES),
		expand(f"{ANALYSIS_DIR}/isoforms/{{sample}}_dominant.csv", sample=SAMPLES),
		expand(f"{ANALYSIS_DIR}/figures/{{gene}}_proportions.png", gene=config["r_analysis"]["genes_of_interest"]),
		expand(f"{ANALYSIS_DIR}/gene_tables/{{gene}}_isoforms.csv", gene=config["r_analysis"]["genes_of_interest"]),
		f"{ANALYSIS_DIR}/report/isoform_report.md",
		f"{LOGS_DIR}/pipeline_complete.txt"

rule download_reference:
	output:
		genome="data/reference/genome.fa",
		transcripts="data/reference/transcripts.fa",
		gtf="data/reference/reference.gtf"
	shell:
		"bash scripts/download_reference_grch38.sh"

rule fastqc:
	input:
		genome="data/reference/genome.fa",
			r1=lambda wc: FASTQ1[wc.sample],
			r2=lambda wc: FASTQ2[wc.sample]
	output:
		f"{LOGS_DIR}/{{sample}}_fastqc.log"
	shell:
		"fastqc {input.r1} {input.r2} > {output} 2>&1"

rule trim_reads:
	input:
		r1=lambda wc: FASTQ1[wc.sample],
		r2=lambda wc: FASTQ2[wc.sample]
	output:
		r1_trimmed=f"{TRIMMED_DIR}/{{sample}}_R1.trimmed.fastq.gz",
		r2_trimmed=f"{TRIMMED_DIR}/{{sample}}_R2.trimmed.fastq.gz"
	log:
		f"{TRIMMED_DIR}/{{sample}}_trim.log"
	params:
		trimmer=config["qc"]["trimmer"]
	shell:
		"fastp -i {input.r1} -I {input.r2} -o {output.r1_trimmed} -O {output.r2_trimmed} > {log} 2>&1"

rule salmon_index:
	input:
		transcripts=config["reference"]["transcripts_fa"]
	output:
		directory(SALMON_INDEX_DIR)
	params:
		kmer=config["salmon"]["kmer_size"]
	shell:
		"salmon index -t {input.transcripts} -i {output} -k {params.kmer}"

rule salmon_quant:
	input:
		index=rules.salmon_index.output,
		r1=f"{TRIMMED_DIR}/{{sample}}_R1.trimmed.fastq.gz",
		r2=f"{TRIMMED_DIR}/{{sample}}_R2.trimmed.fastq.gz"
	output:
		quant=f"{SALMON_QUANT_DIR}/{{sample}}/quant.sf",
		log=f"{SALMON_QUANT_DIR}/{{sample}}/{{sample}}_salmon.log"
	params:
		threads=config["salmon"]["threads"],
		outdir=lambda wc: f"{SALMON_QUANT_DIR}/{wc.sample}"
	shell:
		"salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -p {params.threads} -o {params.outdir} > {output.log} 2>&1"

rule hisat2_index:
	input:
		genome=config["reference"]["genome_fa"]
	output:
		touch(f"{HISAT2_INDEX_DIR}/.done")
	params:
		index_dir=HISAT2_INDEX_DIR,
		threads=config["hisat2"]["threads"]
	shell:
		"hisat2-build -p {params.threads} {input.genome} {params.index_dir}/GRCh38 && touch {output}"

rule hisat2_align:
	input:
		index=rules.hisat2_index.output,
		r1=f"{TRIMMED_DIR}/{{sample}}_R1.trimmed.fastq.gz",
		r2=f"{TRIMMED_DIR}/{{sample}}_R2.trimmed.fastq.gz"
	output:
		bam=f"{ALIGNMENTS_DIR}/{{sample}}.sorted.bam",
		log=f"{ALIGNMENTS_DIR}/{{sample}}_hisat2.log"
	params:
		threads=config["hisat2"]["threads"],
		index_dir=HISAT2_INDEX_DIR
	shell:
		"hisat2 -p {params.threads} --dta -x {params.index_dir}/GRCh38 -1 {input.r1} -2 {input.r2} | samtools sort -@ {params.threads} -o {output.bam} > {output.log} 2>&1 && samtools index {output.bam}"

rule stringtie_assemble:
	input:
		bam=f"{ALIGNMENTS_DIR}/{{sample}}.sorted.bam"
	output:
		gtf=f"{ANALYSIS_DIR}/{{sample}}.stringtie.gtf",
		log=f"{ANALYSIS_DIR}/{{sample}}_stringtie.log"
	params:
		gtf=config["reference"]["gtf"]
	shell:
		"mkdir -p $(dirname {output.gtf}) && stringtie {input.bam} -G {params.gtf} -o {output.gtf} > {output.log} 2>&1"

rule stringtie_merge:
	input:
		gtfs=expand(f"{ANALYSIS_DIR}/{{sample}}.stringtie.gtf", sample=SAMPLES)
	output:
		merged=f"{ANALYSIS_DIR}/merged.stringtie.gtf",
		log=f"{ANALYSIS_DIR}/merged_stringtie_merge.log"
	params:
		gtf=config["reference"]["gtf"]
	shell:
		"mkdir -p $(dirname {output.merged}) && stringtie --merge -G {params.gtf} -o {output.merged} {input.gtfs} > {output.log} 2>&1"

rule gffcompare:
	input:
		merged=f"{ANALYSIS_DIR}/merged.stringtie.gtf"
	output:
		cmp=f"{ANALYSIS_DIR}/merged.merged.stringtie.gtf.tmap",
		log=f"{ANALYSIS_DIR}/merged_gffcompare.log"
	params:
		gtf=config["reference"]["gtf"]
	shell:
		"mkdir -p {ANALYSIS_DIR} && gffcompare -r {params.gtf} -o {ANALYSIS_DIR}/merged {input.merged} > {output.log} 2>&1"

rule pipeline_complete:
	input:
		f"{ANALYSIS_DIR}/merged.merged.stringtie.gtf.tmap"
	output:
		touch(f"{LOGS_DIR}/pipeline_complete.txt")
	shell:
		"mkdir -p $(dirname {output}) && echo 'Pipeline completed successfully.' > {output}"

# Build transcript->gene mapping from reference GTF
rule tx2gene_map:
	input:
		gtf=config["reference"]["gtf"]
	output:
		f"{ANALYSIS_DIR}/tx2gene.tsv"
	shell:
		"mkdir -p $(dirname {output}) && python analysis/scripts/tx2gene.py --gtf {input.gtf} --out {output}"

# Compute isoform proportions and dominant isoforms per sample from Salmon quant + tx2gene map
rule isoform_tables:
	input:
		quant=f"{SALMON_QUANT_DIR}/{{sample}}/quant.sf",
		tx2gene=f"{ANALYSIS_DIR}/tx2gene.tsv"
	output:
		props=f"{ANALYSIS_DIR}/isoforms/{{sample}}_proportions.csv",
		dom=f"{ANALYSIS_DIR}/isoforms/{{sample}}_dominant.csv"
	shell:
		"python analysis/scripts/compute_isoform_tables.py --quant {input.quant} --tx2gene {input.tx2gene} --props {output.props} --dominant {output.dom}"

# Generate figures for genes of interest across samples
rule isoform_figures:
	input:
		props=expand(f"{ANALYSIS_DIR}/isoforms/{{sample}}_proportions.csv", sample=SAMPLES)
	output:
		expand(f"{ANALYSIS_DIR}/figures/{{gene}}_proportions.png", gene=config["r_analysis"]["genes_of_interest"])
	params: 
		genes_str=(" ".join(config["r_analysis"]["genes_of_interest"]) if isinstance(config["r_analysis"]["genes_of_interest"], list) else str(config["r_analysis"]["genes_of_interest"])),
		metric=str(config.get('plots', {}).get('figure_metric', 'proportion'))
	shell:
		"python analysis/scripts/plot_isoform_figures.py --props_dir {ANALYSIS_DIR}/isoforms --genes {params.genes_str} --out_dir {ANALYSIS_DIR}/figures --metric {params.metric}"

# Export gene-specific isoform tables across samples
rule gene_isoform_tables:
	input:
		props=expand(f"{ANALYSIS_DIR}/isoforms/{{sample}}_proportions.csv", sample=SAMPLES)
	output:
		expand(f"{ANALYSIS_DIR}/gene_tables/{{gene}}_isoforms.csv", gene=config["r_analysis"]["genes_of_interest"])
	params:
		genes_str=(" ".join(config["r_analysis"]["genes_of_interest"]) if isinstance(config["r_analysis"]["genes_of_interest"], list) else str(config["r_analysis"]["genes_of_interest"]))
	shell:
		"python analysis/scripts/export_gene_isoforms.py --props_dir {ANALYSIS_DIR}/isoforms --genes {params.genes_str} --out_dir {ANALYSIS_DIR}/gene_tables"

# Create a simple Markdown report linking figures
rule isoform_report:
	input:
		figs=expand(f"{ANALYSIS_DIR}/figures/{{gene}}_proportions.png", gene=config["r_analysis"]["genes_of_interest"])
	output:
		f"{ANALYSIS_DIR}/report/isoform_report.md"
	shell:
		"python analysis/scripts/write_report.py --figs {input.figs} --out {output}"

